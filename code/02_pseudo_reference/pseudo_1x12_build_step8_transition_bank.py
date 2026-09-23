#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pseudo_1x12_build_step8_transition_bank.py
==========================================

FAST Step-8 upstream builder for the ClonoDynamics pseudo-1x12 ensemble.

This utility DOES NOT rerun:
    - Step 2 aggregate-only
    - Step 4
    - Step 5

Instead it reuses the already fitted 66 pair-level Step-2 tables and builds
only the quantities required by the pseudo fluctuation technical-null analysis.

Why this is exact for Step 8
----------------------------
Current production Step 4 propagates the pair-level latent consensus and
reconstructs observed log-frequency directly from count/depth without a
pseudocount.

Current production Step 5:
    - inner-joins endpoint states by aaSeqCDR3;
    - defines xmid_latent as the endpoint midpoint;
    - defines xstar_latent by inverse-variance weighting of endpoint latent
      log-frequency when both endpoint SDs are valid, else midpoint fallback;
    - defines observed displacement as endpoint log-frequency difference;
    - defines common4 as both replicates positive at both endpoints.

Those are exactly the operations implemented here.

Critical optimization
---------------------
There are 2,000 randomized configurations, each containing 15 pseudo-time
intervals, but only a much smaller set of UNIQUE pair-vs-pair technical
comparisons.

This utility therefore builds:

1. 66 compact pair-state caches, once;
2. each unique pair-vs-pair common4 transition block, once;
3. a configuration -> interval -> transition-block map;
4. configuration-level support counts by pseudo-lag.

No configuration-local Step-2/4/5 tables are produced.

Output
------
<pseudo-root>/step8_transition_bank_v1/
    00_build_config.json
    00_pair_state_manifest.csv
    00_transition_block_manifest.csv
    00_configuration_interval_map.csv
    00_configuration_support_by_dt.csv

    pair_states/
        P001_S0.parquet
        ...

    transition_blocks/
        B000001.parquet
        ...

Each transition block contains common4 rows only:

    xmid_latent
    xstar_latent
    xmid_observed
    dx_observed_rep1
    dx_observed_rep2

The configuration interval map contains the physical pseudo-time labels and a
`reverse_sign` flag. Reversing endpoint order changes the sign of both
displacements but leaves covariance, variance and all conditioning coordinates
unchanged.

Typical
-------
python code/utils_4/pseudo_1x12_build_step8_transition_bank.py \
  ./results_4/pseudo \
  --expected-configs 2000

Testing:
  add --max-configs 3

Use --restart-bank only if the 66-pair source tables or this bank contract have
changed.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v1-direct-pairbank-step8-transition-bank-2026-09-17"
BANK_DIRNAME = "step8_transition_bank_v1"
PAIRBANK_DIRNAME = "pair_bank_multirepresentation"

PAIR_STATE_REQUIRED = {
    "aaSeqCDR3",
    "x_latent_median",
    "count_rep1",
    "count_rep2",
    "depth_rep1",
    "depth_rep2",
}

PAIR_STATE_OPTIONAL = {"x_latent_sd"}

PAIR_STATE_OUTPUT_COLUMNS = [
    "aaSeqCDR3",
    "x_latent",
    "x_latent_sd",
    "x_observed_rep1",
    "x_observed_rep2",
    "read_positive_rep1",
    "read_positive_rep2",
]

BLOCK_OUTPUT_COLUMNS = [
    "xmid_latent",
    "xstar_latent",
    "xmid_observed",
    "dx_observed_rep1",
    "dx_observed_rep2",
]


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------

def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if np.isnan(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk_size), b""):
            h.update(block)
    return h.hexdigest()


def stable_hash(payload: Mapping[str, object]) -> str:
    raw = json.dumps(
        json_safe(dict(payload)),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def parquet_columns(path: Path) -> List[str]:
    try:
        import pyarrow.parquet as pq
        return list(pq.ParquetFile(path).schema.names)
    except Exception:
        return list(pd.read_parquet(path).columns)


def pair_key(a: str, b: str) -> Tuple[str, str]:
    return tuple(sorted((str(a), str(b))))  # type: ignore[return-value]


def variant_token(pair_index: int, swapped: bool) -> str:
    return f"P{int(pair_index):03d}_S{int(bool(swapped))}"


def canonical_block_key(
    token_a: str,
    token_b: str,
) -> Tuple[str, str, bool]:
    """
    Return canonical left/right token plus whether actual A->B interval is
    reversed relative to the canonical stored block direction.
    """
    if token_a <= token_b:
        return token_a, token_b, False
    return token_b, token_a, True


def normalize_bool(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    if pd.api.types.is_numeric_dtype(series):
        return (
            pd.to_numeric(series, errors="coerce")
            .fillna(0).ne(0).to_numpy(dtype=bool)
        )
    return (
        series.astype(str).str.strip().str.lower()
        .isin(["true", "t", "1", "yes", "y"])
        .to_numpy(dtype=bool)
    )


# -----------------------------------------------------------------------------
# Manifests
# -----------------------------------------------------------------------------

def load_configuration_manifest(
    pseudo_root: Path,
    expected_configs: int,
) -> pd.DataFrame:
    path = pseudo_root / "pseudo_configurations.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)

    d = pd.read_csv(path, sep="\t")
    required = {
        "configuration_id",
        "pseudo_time",
        "replicate",
        "source_id",
    }
    missing = required - set(d.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required columns {sorted(missing)}"
        )

    d["configuration_id"] = d["configuration_id"].astype(str)
    d["pseudo_time"] = pd.to_numeric(
        d["pseudo_time"], errors="raise"
    ).astype(int)
    d["replicate"] = pd.to_numeric(
        d["replicate"], errors="raise"
    ).astype(int)
    d["source_id"] = d["source_id"].astype(str)

    ids = list(dict.fromkeys(d["configuration_id"].tolist()))
    if int(expected_configs) > 0:
        expected = [
            f"C{i:06d}"
            for i in range(1, int(expected_configs) + 1)
        ]
        if ids != expected:
            raise ValueError(
                "Configuration IDs do not match expected contiguous "
                f"C000001..C{int(expected_configs):06d}"
            )

    for cid, g in d.groupby("configuration_id", sort=False):
        if len(g) != 12:
            raise ValueError(
                f"{cid}: expected 12 assignment rows, found {len(g)}"
            )
        for t, tg in g.groupby("pseudo_time"):
            if set(tg["replicate"].astype(int)) != {1, 2}:
                raise ValueError(
                    f"{cid} T{t}: expected replicate labels 1 and 2"
                )
        if g["source_id"].nunique() != 12:
            raise ValueError(
                f"{cid}: source measurement reused or omitted"
            )

    return d


def load_pair_manifest(
    pseudo_root: Path,
) -> pd.DataFrame:
    path = (
        pseudo_root
        / PAIRBANK_DIRNAME
        / "pair_bank_manifest.tsv"
    )
    if not path.is_file():
        raise FileNotFoundError(path)

    d = pd.read_csv(path, sep="\t")
    required = {
        "pair_index",
        "source_rep1",
        "source_rep2",
    }
    missing = required - set(d.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required columns {sorted(missing)}"
        )

    d["pair_index"] = pd.to_numeric(
        d["pair_index"], errors="raise"
    ).astype(int)
    d["source_rep1"] = d["source_rep1"].astype(str)
    d["source_rep2"] = d["source_rep2"].astype(str)

    if len(d) != 66:
        raise ValueError(
            f"Expected 66 pair-bank entries, found {len(d)}"
        )
    return d


def pair_table_path(
    pseudo_root: Path,
    pair_index: int,
) -> Path:
    base = (
        pseudo_root
        / PAIRBANK_DIRNAME
        / "2-clonotype_state_inference"
        / "per_clone_latent_logP"
    )
    for suffix in (".parquet", ".pq"):
        p = base / f"latent_per_clone_logP_1_{int(pair_index)}{suffix}"
        if p.is_file():
            return p
    raise FileNotFoundError(
        f"Pair table not found for pair_index={pair_index} under {base}"
    )


# -----------------------------------------------------------------------------
# Pair-state compact cache
# -----------------------------------------------------------------------------

def build_pair_state_cache(
    source_path: Path,
    output_path: Path,
    pair_index: int,
    restart: bool,
) -> Dict[str, object]:
    if output_path.is_file() and not restart:
        d = pd.read_parquet(output_path)
        missing = sorted(
            set(PAIR_STATE_OUTPUT_COLUMNS) - set(d.columns)
        )
        if missing:
            raise ValueError(
                f"{output_path}: stale/incomplete pair-state cache; "
                f"missing {missing}. Use --restart-bank."
            )
        if d["aaSeqCDR3"].astype(str).duplicated().any():
            raise ValueError(
                f"{output_path}: duplicate aaSeqCDR3 values"
            )
        return {
            "pair_index": int(pair_index),
            "source_pair_table": str(source_path.resolve()),
            "pair_state": str(output_path.resolve()),
            "n_states": int(len(d)),
            "n_two_positive": int(
                (
                    normalize_bool(d["read_positive_rep1"])
                    & normalize_bool(d["read_positive_rep2"])
                ).sum()
            ),
            "execution": "resumed",
        }

    cols = set(parquet_columns(source_path))
    missing = sorted(PAIR_STATE_REQUIRED - cols)
    if missing:
        raise ValueError(
            f"{source_path}: missing required pair-state columns {missing}"
        )

    read_cols = sorted(
        PAIR_STATE_REQUIRED
        | (PAIR_STATE_OPTIONAL & cols)
    )
    d = pd.read_parquet(source_path, columns=read_cols)

    if d.empty:
        raise ValueError(f"{source_path}: empty pair table")

    d["aaSeqCDR3"] = d["aaSeqCDR3"].astype(str)
    if d["aaSeqCDR3"].duplicated().any():
        raise ValueError(
            f"{source_path}: duplicate aaSeqCDR3 values"
        )

    count1 = pd.to_numeric(d["count_rep1"], errors="coerce").to_numpy(float)
    count2 = pd.to_numeric(d["count_rep2"], errors="coerce").to_numpy(float)
    depth1 = pd.to_numeric(d["depth_rep1"], errors="coerce").to_numpy(float)
    depth2 = pd.to_numeric(d["depth_rep2"], errors="coerce").to_numpy(float)
    x_latent = pd.to_numeric(
        d["x_latent_median"], errors="coerce"
    ).to_numpy(float)

    if (
        np.any(~np.isfinite(count1))
        or np.any(~np.isfinite(count2))
        or np.any(count1 < 0)
        or np.any(count2 < 0)
    ):
        raise ValueError(
            f"{source_path}: invalid replicate counts"
        )
    if (
        np.any(~np.isfinite(depth1))
        or np.any(~np.isfinite(depth2))
        or np.any(depth1 <= 0)
        or np.any(depth2 <= 0)
    ):
        raise ValueError(
            f"{source_path}: invalid replicate depths"
        )
    if np.any(~np.isfinite(x_latent)):
        raise ValueError(
            f"{source_path}: non-finite x_latent_median"
        )

    pos1 = count1 > 0
    pos2 = count2 > 0

    xobs1 = np.full(len(d), np.nan, dtype=float)
    xobs2 = np.full(len(d), np.nan, dtype=float)
    xobs1[pos1] = np.log(count1[pos1] / depth1[pos1])
    xobs2[pos2] = np.log(count2[pos2] / depth2[pos2])

    if "x_latent_sd" in d.columns:
        xsd = pd.to_numeric(
            d["x_latent_sd"], errors="coerce"
        ).to_numpy(float)
    else:
        xsd = np.full(len(d), np.nan, dtype=float)

    out = pd.DataFrame(
        {
            "aaSeqCDR3": d["aaSeqCDR3"].to_numpy(),
            "x_latent": x_latent,
            "x_latent_sd": xsd,
            "x_observed_rep1": xobs1,
            "x_observed_rep2": xobs2,
            "read_positive_rep1": pos1,
            "read_positive_rep2": pos2,
        }
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(
        output_path,
        index=False,
        compression="zstd",
    )

    return {
        "pair_index": int(pair_index),
        "source_pair_table": str(source_path.resolve()),
        "pair_state": str(output_path.resolve()),
        "n_states": int(len(out)),
        "n_two_positive": int((pos1 & pos2).sum()),
        "execution": "computed",
    }


def load_pair_variant(
    path: Path,
    swapped: bool,
) -> pd.DataFrame:
    d = pd.read_parquet(path)
    if not swapped:
        return d

    out = d.copy()
    out[
        ["x_observed_rep1", "x_observed_rep2"]
    ] = out[
        ["x_observed_rep2", "x_observed_rep1"]
    ].to_numpy()
    out[
        ["read_positive_rep1", "read_positive_rep2"]
    ] = out[
        ["read_positive_rep2", "read_positive_rep1"]
    ].to_numpy()
    return out


# -----------------------------------------------------------------------------
# Configuration interval map
# -----------------------------------------------------------------------------

def build_configuration_interval_map(
    configs: pd.DataFrame,
    pair_manifest: pd.DataFrame,
    selected_ids: Sequence[str],
) -> pd.DataFrame:
    lookup = {
        pair_key(row.source_rep1, row.source_rep2):
            int(row.pair_index)
        for row in pair_manifest.itertuples(index=False)
    }

    rows: List[Dict[str, object]] = []

    for cid in selected_ids:
        g = configs[
            configs["configuration_id"].astype(str).eq(str(cid))
        ].copy()

        time_variant: Dict[int, Dict[str, object]] = {}

        for t in range(1, 7):
            tg = (
                g[g["pseudo_time"].eq(t)]
                .sort_values("replicate")
                .reset_index(drop=True)
            )
            if len(tg) != 2:
                raise ValueError(
                    f"{cid} T{t}: expected two replicate assignments"
                )

            a = str(tg.loc[0, "source_id"])
            b = str(tg.loc[1, "source_id"])
            canonical = pair_key(a, b)
            if canonical not in lookup:
                raise KeyError(
                    f"{cid} T{t}: pair {canonical} absent from pair bank"
                )

            pair_index = lookup[canonical]
            swapped = (a, b) != canonical
            token = variant_token(pair_index, swapped)

            time_variant[t] = {
                "pair_index": pair_index,
                "swapped": bool(swapped),
                "token": token,
                "source_rep1": a,
                "source_rep2": b,
            }

        for t0 in range(1, 6):
            for t1 in range(t0 + 1, 7):
                a = time_variant[t0]
                b = time_variant[t1]
                left, right, reverse = canonical_block_key(
                    str(a["token"]),
                    str(b["token"]),
                )

                rows.append(
                    {
                        "configuration_id": cid,
                        "t0": int(t0),
                        "t1": int(t1),
                        "dt": int(t1 - t0),
                        "token_t0": str(a["token"]),
                        "token_t1": str(b["token"]),
                        "block_left_token": left,
                        "block_right_token": right,
                        "reverse_sign": bool(reverse),
                    }
                )

    out = pd.DataFrame(rows)
    if len(out) != 15 * len(selected_ids):
        raise RuntimeError(
            "Configuration interval map has unexpected row count"
        )
    return out


# -----------------------------------------------------------------------------
# Unique transition blocks
# -----------------------------------------------------------------------------

def parse_variant_token(token: str) -> Tuple[int, bool]:
    # P001_S0
    if not token.startswith("P") or "_S" not in token:
        raise ValueError(f"Invalid pair-variant token: {token}")
    left, right = token.split("_S", 1)
    return int(left[1:]), bool(int(right))


def transition_block_id(
    left_token: str,
    right_token: str,
) -> str:
    raw = f"{left_token}|{right_token}"
    digest = hashlib.blake2b(
        raw.encode("utf-8"),
        digest_size=8,
    ).hexdigest()
    return f"B_{digest}"


def build_transition_block(
    left_token: str,
    right_token: str,
    pair_state_dir: Path,
    output_path: Path,
    restart: bool,
) -> Dict[str, object]:
    meta_path = output_path.with_suffix(".meta.json")
    if output_path.is_file() and not restart:
        d = pd.read_parquet(output_path)
        missing = sorted(
            set(BLOCK_OUTPUT_COLUMNS) - set(d.columns)
        )
        if missing:
            raise ValueError(
                f"{output_path}: incomplete transition block; "
                f"missing {missing}. Use --restart-bank."
            )
        if not meta_path.is_file():
            raise ValueError(
                f"{output_path}: transition block metadata missing. "
                "Use --restart-bank once to rebuild the bank under the current contract."
            )
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if (
            str(meta.get("left_token")) != left_token
            or str(meta.get("right_token")) != right_token
        ):
            raise ValueError(
                f"{meta_path}: transition-block key mismatch. Use --restart-bank."
            )
        if int(meta.get("n_common4", -1)) != int(len(d)):
            raise ValueError(
                f"{meta_path}: n_common4 does not match Parquet row count."
            )
        meta["block_path"] = str(output_path.resolve())
        meta["execution"] = "resumed"
        return meta

    left_index, left_swap = parse_variant_token(left_token)
    right_index, right_swap = parse_variant_token(right_token)

    left_path = pair_state_dir / f"P{left_index:03d}_S0.parquet"
    right_path = pair_state_dir / f"P{right_index:03d}_S0.parquet"

    left = load_pair_variant(left_path, left_swap)
    right = load_pair_variant(right_path, right_swap)

    # Step-5 transition universe: inner join by aaSeqCDR3.
    l = left.rename(
        columns={
            c: f"{c}_L"
            for c in left.columns
            if c != "aaSeqCDR3"
        }
    )
    r = right.rename(
        columns={
            c: f"{c}_R"
            for c in right.columns
            if c != "aaSeqCDR3"
        }
    )

    joined = l.merge(
        r,
        on="aaSeqCDR3",
        how="inner",
        validate="one_to_one",
        sort=False,
    )
    n_transitions = int(len(joined))

    p1_l = normalize_bool(joined["read_positive_rep1_L"])
    p2_l = normalize_bool(joined["read_positive_rep2_L"])
    p1_r = normalize_bool(joined["read_positive_rep1_R"])
    p2_r = normalize_bool(joined["read_positive_rep2_R"])

    common4 = p1_l & p2_l & p1_r & p2_r
    c = joined.loc[common4].copy()

    if c.empty:
        raise RuntimeError(
            f"No common4 rows for {left_token} vs {right_token}"
        )

    x0 = pd.to_numeric(
        c["x_latent_L"], errors="coerce"
    ).to_numpy(float)
    x1 = pd.to_numeric(
        c["x_latent_R"], errors="coerce"
    ).to_numpy(float)
    sd0 = pd.to_numeric(
        c["x_latent_sd_L"], errors="coerce"
    ).to_numpy(float)
    sd1 = pd.to_numeric(
        c["x_latent_sd_R"], errors="coerce"
    ).to_numpy(float)

    xA0 = pd.to_numeric(
        c["x_observed_rep1_L"], errors="coerce"
    ).to_numpy(float)
    xA1 = pd.to_numeric(
        c["x_observed_rep1_R"], errors="coerce"
    ).to_numpy(float)
    xB0 = pd.to_numeric(
        c["x_observed_rep2_L"], errors="coerce"
    ).to_numpy(float)
    xB1 = pd.to_numeric(
        c["x_observed_rep2_R"], errors="coerce"
    ).to_numpy(float)

    # common4 guarantees these observed endpoint values are defined.
    if not (
        np.isfinite(xA0).all()
        and np.isfinite(xA1).all()
        and np.isfinite(xB0).all()
        and np.isfinite(xB1).all()
    ):
        raise RuntimeError(
            f"{left_token} vs {right_token}: common4 contains non-finite "
            "observed log-frequencies"
        )

    xmid = 0.5 * (x0 + x1)

    valid_sd = (
        np.isfinite(sd0)
        & np.isfinite(sd1)
        & (sd0 > 0)
        & (sd1 > 0)
    )
    xstar = xmid.copy()
    if valid_sd.any():
        w0 = 1.0 / (sd0[valid_sd] ** 2 + 1e-12)
        w1 = 1.0 / (sd1[valid_sd] ** 2 + 1e-12)
        xstar[valid_sd] = (
            w0 * x0[valid_sd]
            + w1 * x1[valid_sd]
        ) / (w0 + w1)

    xmid_observed = 0.25 * (xA0 + xA1 + xB0 + xB1)

    out = pd.DataFrame(
        {
            "xmid_latent": xmid,
            "xstar_latent": xstar,
            "xmid_observed": xmid_observed,
            "dx_observed_rep1": xA1 - xA0,
            "dx_observed_rep2": xB1 - xB0,
        }
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(
        output_path,
        index=False,
        compression="zstd",
    )

    meta = {
        "left_token": left_token,
        "right_token": right_token,
        "block_path": str(output_path.resolve()),
        "n_transitions": n_transitions,
        "n_common4": int(len(out)),
        "fraction_common4": float(
            len(out) / n_transitions
            if n_transitions else np.nan
        ),
        "xstar_inverse_variance_rows": int(valid_sd.sum()),
        "xstar_midpoint_fallback_rows": int((~valid_sd).sum()),
    }
    meta_path.write_text(
        json.dumps(json_safe(meta), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {**meta, "execution": "computed"}


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build a fast reusable Step-8 transition bank directly from the "
            "already fitted 66-pair pseudo technical bank."
        ),
    )
    p.add_argument("pseudo_root", type=Path)
    p.add_argument("--expected-configs", type=int, default=2000)
    p.add_argument(
        "--max-configs",
        type=int,
        default=None,
        help=(
            "Testing only. The transition bank is built only for unique "
            "blocks required by the first N configurations."
        ),
    )
    p.add_argument(
        "--restart-bank",
        action="store_true",
        help="Rebuild pair-state and transition-block caches.",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()

    if args.expected_configs < 0:
        raise ValueError("--expected-configs must be >=0")
    if args.max_configs is not None and args.max_configs < 1:
        raise ValueError("--max-configs must be >=1")

    root = args.pseudo_root.expanduser().resolve(strict=True)
    bank_root = root / BANK_DIRNAME
    pair_state_dir = bank_root / "pair_states"
    block_dir = bank_root / "transition_blocks"

    if args.restart_bank and bank_root.exists():
        shutil.rmtree(bank_root)

    pair_state_dir.mkdir(parents=True, exist_ok=True)
    block_dir.mkdir(parents=True, exist_ok=True)

    configs = load_configuration_manifest(
        root,
        int(args.expected_configs),
    )
    pair_manifest = load_pair_manifest(root)

    all_ids = list(
        dict.fromkeys(configs["configuration_id"].astype(str).tolist())
    )
    selected_ids = (
        all_ids
        if args.max_configs is None
        else all_ids[: int(args.max_configs)]
    )

    print("[INFO] configurations available :", len(all_ids))
    print("[INFO] configurations selected  :", len(selected_ids))
    print("[INFO] pair fits reused          : 66")
    print("[INFO] Step 2 aggregate-only     : SKIPPED")
    print("[INFO] Step 4                    : SKIPPED")
    print("[INFO] Step 5                    : SKIPPED")

    # Pair-state compact cache.
    pair_state_rows = []
    for row in pair_manifest.sort_values("pair_index").itertuples(index=False):
        pair_index = int(row.pair_index)
        source = pair_table_path(root, pair_index)
        out = pair_state_dir / f"P{pair_index:03d}_S0.parquet"
        audit = build_pair_state_cache(
            source,
            out,
            pair_index,
            restart=bool(args.restart_bank),
        )
        pair_state_rows.append(audit)

    pair_state_manifest = pd.DataFrame(pair_state_rows)
    pair_state_manifest.to_csv(
        bank_root / "00_pair_state_manifest.csv",
        index=False,
    )

    # Configuration interval -> unique block mapping.
    interval_map = build_configuration_interval_map(
        configs,
        pair_manifest,
        selected_ids,
    )

    unique_keys = (
        interval_map[
            ["block_left_token", "block_right_token"]
        ]
        .drop_duplicates()
        .sort_values(
            ["block_left_token", "block_right_token"]
        )
        .reset_index(drop=True)
    )

    print("[INFO] physical intervals requested :", len(interval_map))
    print("[INFO] unique transition blocks     :", len(unique_keys))

    block_rows = []
    block_lookup: Dict[Tuple[str, str], Dict[str, object]] = {}

    for i, row in unique_keys.iterrows():
        left = str(row["block_left_token"])
        right = str(row["block_right_token"])
        bid = transition_block_id(left, right)
        out = block_dir / f"{bid}.parquet"

        audit = build_transition_block(
            left,
            right,
            pair_state_dir,
            out,
            restart=bool(args.restart_bank),
        )
        audit["block_id"] = bid
        block_rows.append(audit)
        block_lookup[(left, right)] = audit

        if (i + 1) % 50 == 0 or i == len(unique_keys) - 1:
            print(
                f"[block {i+1}/{len(unique_keys)}] {bid} "
                f"{left} vs {right}"
            )

    block_manifest = pd.DataFrame(block_rows)
    block_manifest.to_csv(
        bank_root / "00_transition_block_manifest.csv",
        index=False,
    )

    # Enrich interval map with block IDs and support counts.
    block_meta = block_manifest[
        [
            "left_token",
            "right_token",
            "block_id",
            "block_path",
            "n_transitions",
            "n_common4",
            "fraction_common4",
        ]
    ].rename(
        columns={
            "left_token": "block_left_token",
            "right_token": "block_right_token",
        }
    )

    interval_map = interval_map.merge(
        block_meta,
        on=["block_left_token", "block_right_token"],
        how="left",
        validate="many_to_one",
    )

    if interval_map["block_id"].isna().any():
        raise RuntimeError(
            "Some configuration intervals could not be mapped to a transition block"
        )

    interval_map.to_csv(
        bank_root / "00_configuration_interval_map.csv",
        index=False,
    )

    support = (
        interval_map
        .groupby(["configuration_id", "dt"], as_index=False)
        .agg(
            n_transitions=("n_transitions", "sum"),
            n_common4=("n_common4", "sum"),
            n_intervals=("block_id", "size"),
        )
    )
    support["fraction_common4"] = (
        support["n_common4"] / support["n_transitions"]
    )
    support.to_csv(
        bank_root / "00_configuration_support_by_dt.csv",
        index=False,
    )

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "pseudo_root": str(root),
        "bank_root": str(bank_root),
        "source_biological_subjects": 1,
        "source_measurements": 12,
        "pair_fits_reused": 66,
        "step2_aggregate_only_rerun": False,
        "step4_rerun": False,
        "step5_rerun": False,
        "n_configurations_available": int(len(all_ids)),
        "n_configurations_selected": int(len(selected_ids)),
        "n_physical_intervals_selected": int(len(interval_map)),
        "n_unique_transition_blocks": int(len(block_manifest)),
        "observed_log_frequency_definition":
            "log(count/depth) for positive counts; no pseudocount",
        "transition_universe":
            "inner join of pair-state tables by aaSeqCDR3",
        "common4_definition":
            "both technical replicates positive at both endpoints",
        "xmid_latent_definition":
            "0.5*(x_latent_t0+x_latent_t1)",
        "xstar_latent_definition":
            "inverse-variance endpoint weighted when both SD valid; midpoint fallback",
        "xmid_observed_definition":
            "0.25*(xA0+xA1+xB0+xB1) on common4",
        "configuration_id_is_biological_subject": False,
    }

    (bank_root / "00_build_config.json").write_text(
        json.dumps(
            json_safe(config),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print("\n[DONE] fast Step-8 transition bank")
    print("Selected configurations :", len(selected_ids))
    print("Unique transition blocks:", len(block_manifest))
    print("Bank root               :", bank_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
