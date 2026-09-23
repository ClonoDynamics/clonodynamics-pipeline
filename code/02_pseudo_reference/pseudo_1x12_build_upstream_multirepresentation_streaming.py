#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pseudo_1x12_build_upstream_multirepresentation.py
=================================================

Build the production ClonoDynamics pseudo-1x12 upstream ensemble under the
replicate-resolved multirepresentation framework.

DESIGN
------
The pseudo reference contains ONE biological source sample from which PBMCs
were isolated together and then split into 12 independently processed aliquots.
The resulting 12 TCRbeta RepSeq measurements are treated as exchangeable
technical measurements of one biological state.

The design manifest is generated upstream by `pseudo_1x12_generate.py`:

    1 biological source sample
    12 technical measurements
        -> randomized assignment
        -> 6 pseudo-timepoints x 2 technical replicates
        -> N randomization configurations (production: 2,000)

`configuration_id` is a randomization index. It is never a biological subject
and never contributes to biological N.

WHY A 66-PAIR BANK IS VALID
---------------------------
With 12 source measurements there are C(12,2)=66 unordered technical pairs.
The Step-2 paired Negative-Binomial latent fit is pair-specific and does not
depend on the arbitrary pseudo-time label. Therefore each unordered pair is
fitted only once.

For each pseudo configuration this script then:

    1. selects the six required pair fits from the 66-pair bank;
    2. restores the randomized A/B orientation encoded in the configuration;
    3. remaps the six pair states to pseudo-times 1..6;
    4. writes a configuration-local six-pair fit manifest, remapped to
       pseudo-times 1..6;
    5. reruns Step-2 AGGREGATION ONLY on those six pair states;
    6. thereby recomputes configuration-specific empirical p-values and
       operational observability;
    6. runs the current Step 4 multirepresentation trajectory assembly;
    7. runs the current Step 5 generic longitudinal transition assembly.

No pair-level latent fit is repeated across configurations.

The configuration-local fit manifest is essential for recent Step-2 versions
that protect aggregate-only mode by restricting input to pair IDs marked as
successful in `<fit-cache-dir>/latent_params_by_pair.csv`. Pointing that guard
at the 66-pair bank while supplying only six remapped pseudo-time files would
incorrectly make Step 2 expect all 66 cached files in every configuration.
This builder therefore supplies a six-pair configuration-local cache manifest.

IMPORTANT OBSERVABILITY RULE
----------------------------
The fitted pair posterior can be cached, but the empirical p-value reference
pool / operational observation annotation cannot be copied from the 66-pair
bank. Those quantities depend on the six pseudo-time states contained in the
specific configuration and are therefore recomputed by Step 2 in
`--aggregate-only` mode for every configuration.

CURRENT PRODUCTION CONTRACT
---------------------------
Required code in --code-dir:

    2-multirepresentation_clonotype_state_inference.py
    4-multirepresentation_trajectory_assembly.py
    5-longitudinal_transition_assembly.py
    noiseK_latent.py

Per-configuration output:

    <pseudo-root>/ensemble_results_v3/C000001/
        2-clonotype_state_inference/
            per_clone_latent_subject.parquet
            per_clone_latent_logP/
        4-multirepresentation_trajectory_assembly/
            multirepresentation_trajectories_long.parquet
        5-longitudinal_transition_assembly/
            longitudinal_transitions.parquet

The Step-5 transition table is the direct input to the new Step 7 and Step 8
pseudo technical-null analyses.

STALE-CACHE PROTECTION
----------------------
The 66-pair bank is guarded by a stable pair-bank contract signature. Builder
changes that affect only configuration assembly (such as the six-pair
aggregate-only cache fix in v2) do not invalidate already completed pair fits.

The signature contains:
    - hashes of the current Step-2 and noise-model scripts;
    - source-manifest hash;
    - alpha;
    - pair-bank filename pattern.

If an existing pair bank does not match the current signature, execution stops
and requires --restart-pairbank. This prevents an old latent-only Step-2 cache
from being silently reused after the multirepresentation framework update.

TYPICAL PRODUCTION RUN
----------------------
python3 pseudo_1x12_build_upstream_multirepresentation.py \
    /path/to/pseudo_1x12 \
    --code-dir /path/to/clonodynamics-pipeline/code/dynamics \
    --expected-configs 2000 \
    --alpha 0.05 \
    --n-jobs 0 \
    --n-jobs-profile balanced

Testing:
    add --max-configs 3

The generator manifests themselves are not modified by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


SCRIPT_VERSION = "v3-step7-compact-cache-2026-09-15"
PAIRBANK_CONTRACT_VERSION = "v1-multirepresentation-ensemble-v3-2026-09-10"

STEP2_NAME = "2-multirepresentation_clonotype_state_inference.py"
STEP4_NAME = "4-multirepresentation_trajectory_assembly.py"
STEP5_NAME = "5-longitudinal_transition_assembly.py"
NOISE_NAME = "noiseK_latent.py"

ENSEMBLE_DIRNAME = "ensemble_results_v3"
PAIRBANK_DIRNAME = "pair_bank_multirepresentation"

CONFIG_RE = re.compile(r"^C\d{6}$")

PAIR_PATTERN = (
    r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])"
    r"(?:\.(?:tsv|csv|parquet|pq|feather|arrow))?$"
)

STEP7_REQUIRED_TRANSITION_COLUMNS = {
    "configuration_id",
    "subject",
    "aaSeqCDR3",
    "t0",
    "t1",
    "dt",
    "x0_latent",
    "xmid_latent",
    "xstar_latent",
    "dx_latent",
    "x_observed_rep1_t0",
    "x_observed_rep2_t0",
    "dx_observed_rep1",
    "dx_observed_rep2",
    "forward_ab_eligible",
    "forward_ba_eligible",
    "common4",
    "p_value_t0",
    "p_value_t1",
    "observable_reference_t0",
    "observable_reference_t1",
    "obs_class_reference",
}

# Minimal cache consumed by Step 7. It keeps only the primary lag and the
# columns actually used by the pseudo forward technical-null analysis.
STEP7_COMPACT_COLUMNS = [
    "configuration_id",
    "subject",
    "dt",
    "x0_latent",
    "xmid_latent",
    "xstar_latent",
    "dx_latent",
    "x_observed_rep1_t0",
    "x_observed_rep2_t0",
    "dx_observed_rep1",
    "dx_observed_rep2",
    "forward_ab_eligible",
    "forward_ba_eligible",
    "common4",
]
STEP7_COMPACT_REL = Path("step7_cache") / "step7_dt1.parquet"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            h.update(block)
    return h.hexdigest()


def stable_hash(payload: Dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def run(cmd: Sequence[object], log: Path, cwd: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    print("\n$ " + " ".join(map(str, cmd)))
    print("[log]", log)
    with log.open("w", encoding="utf-8", buffering=1) as handle:
        process = subprocess.Popen(
            [str(x) for x in cmd],
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            handle.write(line)
        rc = int(process.wait())
    if rc != 0:
        raise RuntimeError(f"Command failed ({rc}). See:\n  {log}")


def safe_symlink(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    # Use an absolute target because the original measurement directory is not
    # necessarily inside the pseudo project.
    destination.symlink_to(source.resolve())


def load_source_manifest(root: Path) -> pd.DataFrame:
    path = root / "source_measurements.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Run pseudo_1x12_generate.py first."
        )
    frame = pd.read_csv(path, sep="\t")
    required = {"source_id", "source_path"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if len(frame) != 12:
        raise ValueError(f"Expected 12 source measurements, found {len(frame)}")
    if frame["source_id"].astype(str).duplicated().any():
        raise ValueError("source_measurements.tsv contains duplicate source_id values")
    for raw in frame["source_path"]:
        p = Path(str(raw)).expanduser()
        if not p.is_file():
            raise FileNotFoundError(p)
    return frame


def load_configuration_manifest(
    root: Path,
    expected_configs: int,
) -> pd.DataFrame:
    path = root / "pseudo_configurations.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Run pseudo_1x12_generate.py first."
        )

    frame = pd.read_csv(path, sep="\t")
    required = {
        "configuration_id",
        "pseudo_time",
        "replicate",
        "source_id",
        "source_path",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    frame["configuration_id"] = frame["configuration_id"].astype(str)
    frame["pseudo_time"] = pd.to_numeric(frame["pseudo_time"], errors="raise").astype(int)
    frame["replicate"] = pd.to_numeric(frame["replicate"], errors="raise").astype(int)

    ids = list(dict.fromkeys(frame["configuration_id"].tolist()))
    if expected_configs > 0:
        expected = [f"C{i:06d}" for i in range(1, expected_configs + 1)]
        if ids != expected:
            raise ValueError(
                "Configuration IDs do not match the expected contiguous "
                f"C000001..C{expected_configs:06d} production set."
            )

    for cid, group in frame.groupby("configuration_id", sort=False):
        if len(group) != 12:
            raise ValueError(f"{cid}: expected 12 assignment rows, found {len(group)}")
        if sorted(group["pseudo_time"].tolist()) != [1,1,2,2,3,3,4,4,5,5,6,6]:
            raise ValueError(f"{cid}: invalid pseudo-time geometry")
        for t, tg in group.groupby("pseudo_time"):
            if set(tg["replicate"].tolist()) != {1, 2}:
                raise ValueError(f"{cid} T{t}: expected replicate labels 1 and 2")
        used = group["source_id"].astype(str).tolist()
        if len(set(used)) != 12:
            raise ValueError(f"{cid}: a source measurement was reused or omitted")

    # Exact assignment duplicates are not allowed by the generator.
    wide_signature = (
        frame.sort_values(["configuration_id", "pseudo_time", "replicate"])
        .groupby("configuration_id")["source_id"]
        .apply(lambda s: "|".join(s.astype(str)))
    )
    if wide_signature.duplicated().any():
        dup = wide_signature[wide_signature.duplicated(keep=False)]
        raise ValueError(
            "Exact duplicate pseudo configurations found: "
            + ", ".join(dup.index[:10].astype(str))
        )

    return frame


def pair_key(a: str, b: str) -> Tuple[str, str]:
    return tuple(sorted((str(a), str(b))))  # type: ignore[return-value]


def build_pair_bank_raw(root: Path, source: pd.DataFrame) -> pd.DataFrame:
    bank_root = root / PAIRBANK_DIRNAME
    raw_dir = bank_root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    by_id = {str(row.source_id): row for row in source.itertuples(index=False)}
    rows: List[Dict[str, object]] = []

    for pair_index, (a, b) in enumerate(
        itertools.combinations(sorted(by_id), 2),
        start=1,
    ):
        pa = Path(str(by_id[a].source_path)).expanduser().resolve(strict=True)
        pb = Path(str(by_id[b].source_path)).expanduser().resolve(strict=True)

        safe_symlink(pa, raw_dir / f"1_{pair_index}-1{pa.suffix.lower()}")
        safe_symlink(pb, raw_dir / f"1_{pair_index}-2{pb.suffix.lower()}")

        rows.append(
            {
                "pair_index": pair_index,
                "pair_id": f"P{pair_index:03d}",
                "source_rep1": a,
                "source_rep2": b,
                "source_path_rep1": str(pa),
                "source_path_rep2": str(pb),
            }
        )

    frame = pd.DataFrame(rows)
    if len(frame) != 66:
        raise RuntimeError(f"Expected 66 unique source pairs, obtained {len(frame)}")
    frame.to_csv(bank_root / "pair_bank_manifest.tsv", sep="\t", index=False)
    return frame


def pairbank_signature_payload(
    root: Path,
    step2: Path,
    noise: Path,
    alpha: float,
) -> Dict[str, object]:
    return {
        "builder_version": PAIRBANK_CONTRACT_VERSION,
        "step2_sha256": sha256_file(step2),
        "noise_sha256": sha256_file(noise),
        "source_manifest_sha256": sha256_file(root / "source_measurements.tsv"),
        "alpha": float(alpha),
        "pair_pattern": PAIR_PATTERN,
        "n_source_measurements": 12,
        "n_pair_fits": 66,
    }


def pairbank_complete(bank_results: Path) -> bool:
    params = bank_results / "latent_params_by_pair.csv"
    pair_dir = bank_results / "per_clone_latent_logP"
    if not params.is_file() or not pair_dir.is_dir():
        return False
    tables = list(pair_dir.glob("latent_per_clone_logP_1_*.parquet"))
    return len(tables) == 66


def ensure_pairbank(
    root: Path,
    code: Path,
    step2: Path,
    noise: Path,
    alpha: float,
    n_jobs: int,
    n_jobs_profile: str,
    restart: bool,
) -> Tuple[Path, Path, pd.DataFrame]:
    bank_root = root / PAIRBANK_DIRNAME
    bank_results = bank_root / "2-clonotype_state_inference"
    bank_pair_dir = bank_results / "per_clone_latent_logP"
    signature_path = bank_root / "00_pairbank_signature.json"

    source = load_source_manifest(root)
    pair_manifest = build_pair_bank_raw(root, source)

    expected_payload = pairbank_signature_payload(root, step2, noise, alpha)
    expected_signature = stable_hash(expected_payload)

    if restart and bank_results.exists():
        shutil.rmtree(bank_results)
    if restart and signature_path.exists():
        signature_path.unlink()

    if bank_results.exists():
        if not signature_path.is_file():
            raise RuntimeError(
                "Existing pair bank has no multirepresentation signature. "
                "Use --restart-pairbank."
            )
        previous = json.loads(signature_path.read_text(encoding="utf-8"))
        if previous.get("signature") != expected_signature:
            raise RuntimeError(
                "Existing pair bank was created with different source/code/alpha. "
                "Use --restart-pairbank to refit all 66 pairs."
            )

    if not pairbank_complete(bank_results):
        cmd = [
            sys.executable,
            "-u",
            str(step2),
            "--data-dir",
            str(bank_root / "raw"),
            "--results-dir",
            str(bank_results),
            "--file-sep",
            "\t",
            "--pattern",
            PAIR_PATTERN,
            "--alpha",
            repr(float(alpha)),
            "--null",
            "subject",
            "--tail",
            "low",
            "--n-jobs",
            str(int(n_jobs)),
            "--n-jobs-profile",
            str(n_jobs_profile),
            "--intermediate-format",
            "parquet",
            "--parquet-compression",
            "zstd",
        ]
        run(cmd, root / "logs" / "pairbank_step2.log", code)

        if not pairbank_complete(bank_results):
            raise RuntimeError("66-pair multirepresentation Step-2 bank is incomplete")

        signature_path.write_text(
            json.dumps(
                {
                    "signature": expected_signature,
                    "payload": expected_payload,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    else:
        print("[resume] current 66-pair multirepresentation bank is complete")

    return bank_results, bank_pair_dir, pair_manifest


def swap_rep_orientation(frame: pd.DataFrame) -> pd.DataFrame:
    """Swap every matched rep1/rep2 column pair and preserve all other fields."""
    out = frame.copy()
    done = set()
    for c1 in list(out.columns):
        if c1 in done or "rep1" not in c1:
            continue
        c2 = c1.replace("rep1", "rep2")
        if c2 == c1 or c2 not in out.columns:
            continue
        v1 = out[c1].copy()
        out[c1] = out[c2].to_numpy()
        out[c2] = v1.to_numpy()
        done.update({c1, c2})
    return out


def find_pair_table(pair_dir: Path, pair_index: int) -> Path:
    for suffix in (".parquet", ".pq"):
        path = pair_dir / f"latent_per_clone_logP_1_{pair_index}{suffix}"
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"No cached per-clone table for pair {pair_index} under {pair_dir}"
    )


def find_pair_posterior(pair_dir: Path, pair_index: int) -> Optional[Path]:
    path = pair_dir / f"latent_per_clone_logP_1_{pair_index}.posterior.npz"
    return path if path.is_file() else None



def build_configuration_fit_cache(
    bank_results: Path,
    step2_dir: Path,
    map_rows: Sequence[Dict[str, object]],
) -> Path:
    """
    Build the six-pair fit manifest expected by Step-2 aggregate-only mode.

    The 66-pair bank contains successful fits indexed as 1_1 ... 1_66.
    A pseudo configuration contains only six selected pair tables, remapped to
    pseudo-time IDs 1_1 ... 1_6. Recent Step-2 implementations use
    `<fit-cache-dir>/latent_params_by_pair.csv` as a stale-cache/success guard.
    Therefore the fit-cache manifest passed to aggregate-only must describe
    these SIX remapped pairs, not all 66 bank pairs.

    No model parameter is re-estimated here. Pair-level parameters are copied
    from the corresponding bank row and identifiers are remapped to pseudo-time.
    Rep1/rep2-labelled metadata are swapped when the configuration orientation
    is reversed.
    """
    bank_params_path = Path(bank_results) / "latent_params_by_pair.csv"
    if not bank_params_path.is_file():
        raise FileNotFoundError(
            f"Missing pair-bank parameter table: {bank_params_path}"
        )

    bank = pd.read_csv(bank_params_path)
    if bank.empty:
        raise ValueError(f"Empty pair-bank parameter table: {bank_params_path}")

    # Normalize identifiers used to locate the bank row.
    if "subject" in bank.columns:
        bank["subject"] = pd.to_numeric(bank["subject"], errors="coerce")
    if "time" in bank.columns:
        bank["time"] = pd.to_numeric(bank["time"], errors="coerce")

    config_rows: List[pd.Series] = []

    def swap_rep_metadata(row: pd.Series) -> pd.Series:
        row = row.copy()
        done = set()
        for c1 in list(row.index):
            if c1 in done or "rep1" not in str(c1):
                continue
            c2 = str(c1).replace("rep1", "rep2")
            if c2 == c1 or c2 not in row.index:
                continue
            v1 = row[c1]
            row[c1] = row[c2]
            row[c2] = v1
            done.update({str(c1), c2})
        return row

    for mapping in sorted(map_rows, key=lambda r: int(r["pseudo_time"])):
        pseudo_time = int(mapping["pseudo_time"])
        pair_index = int(mapping["pair_bank_index"])
        orientation_swapped = bool(int(mapping["orientation_swapped"]))

        candidates = bank.copy()
        if {"subject", "time"}.issubset(candidates.columns):
            hit = candidates[
                (pd.to_numeric(candidates["subject"], errors="coerce") == 1)
                & (pd.to_numeric(candidates["time"], errors="coerce") == pair_index)
            ]
        elif "pair_id" in candidates.columns:
            hit = candidates[
                candidates["pair_id"].astype(str).eq(f"1_{pair_index}")
            ]
        else:
            raise ValueError(
                "Pair-bank latent_params_by_pair.csv must contain either "
                "subject/time or pair_id."
            )

        if len(hit) != 1:
            raise ValueError(
                f"Expected exactly one parameter row for pair-bank index "
                f"{pair_index}; found {len(hit)}"
            )

        row = hit.iloc[0].copy()
        if orientation_swapped:
            row = swap_rep_metadata(row)

        # Remap bank pair identity -> configuration pseudo-time identity.
        if "subject" in row.index:
            row["subject"] = 1
        if "time" in row.index:
            row["time"] = pseudo_time
        if "pair_id" in row.index:
            row["pair_id"] = f"1_{pseudo_time}"

        # Preserve bank provenance explicitly.
        row["configuration_id"] = str(mapping["configuration_id"])
        row["pseudo_time"] = pseudo_time
        row["pair_bank_index"] = pair_index
        row["pair_bank_pair_id"] = f"1_{pair_index}"
        row["source_rep1_id"] = str(mapping["replicate1_source"])
        row["source_rep2_id"] = str(mapping["replicate2_source"])
        row["orientation_swapped_from_pair_bank"] = int(orientation_swapped)

        # Ensure the guard sees the remapped pair as successful.
        if "success" not in row.index:
            row["success"] = True

        config_rows.append(row)

    config_params = pd.DataFrame(config_rows)
    if len(config_params) != 6:
        raise RuntimeError(
            f"Configuration-local fit cache must contain 6 rows, got "
            f"{len(config_params)}"
        )

    expected_ids = {f"1_{t}" for t in range(1, 7)}
    if "pair_id" in config_params.columns:
        got_ids = set(config_params["pair_id"].astype(str))
        if got_ids != expected_ids:
            raise RuntimeError(
                f"Configuration-local fit-cache pair IDs are {sorted(got_ids)}, "
                f"expected {sorted(expected_ids)}"
            )

    # Only successful fits are valid in a configuration assembled from the
    # successful 66-pair bank.
    if "success" in config_params.columns:
        success = (
            config_params["success"]
            if pd.api.types.is_bool_dtype(config_params["success"])
            else config_params["success"].astype(str).str.strip().str.lower()
                 .isin(["true", "t", "1", "yes", "y"])
        )
        if not bool(success.all()):
            raise RuntimeError(
                "A configuration-local pair was not marked successful."
            )

    out = Path(step2_dir) / "latent_params_by_pair.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    config_params.to_csv(out, index=False)
    return out


def transition_schema(path: Path) -> List[str]:
    try:
        import pyarrow.parquet as pq
        return list(pq.ParquetFile(path).schema.names)
    except Exception:
        return list(pd.read_parquet(path).columns)


def validate_transition_output(path: Path, cid: str) -> Dict[str, object]:
    columns = set(transition_schema(path))
    missing = sorted(STEP7_REQUIRED_TRANSITION_COLUMNS - columns)
    if missing:
        raise ValueError(
            f"{cid}: Step-5 transition table is not Step-7 ready; missing {missing}"
        )

    cols = [
        "configuration_id",
        "subject",
        "t0",
        "t1",
        "dt",
        "forward_ab_eligible",
        "forward_ba_eligible",
        "common4",
    ]
    d = pd.read_parquet(path, columns=cols)

    configs = set(d["configuration_id"].astype(str).dropna().unique())
    if configs != {cid}:
        raise ValueError(f"{cid}: configuration_id values in output are {sorted(configs)}")

    subjects = pd.to_numeric(d["subject"], errors="coerce").dropna().astype(int).unique()
    if set(subjects) != {1}:
        raise ValueError(f"{cid}: expected pseudo biological subject 1, found {subjects}")

    pairs = d[["t0", "t1"]].drop_duplicates()
    if len(pairs) != 15:
        raise ValueError(f"{cid}: expected 15 pseudo-time pairs, found {len(pairs)}")

    dt_values = sorted(pd.to_numeric(d["dt"], errors="coerce").dropna().astype(int).unique())
    if dt_values != [1, 2, 3, 4, 5]:
        raise ValueError(f"{cid}: expected dt 1..5, found {dt_values}")

    def n_true(column: str) -> int:
        s = d[column]
        if pd.api.types.is_bool_dtype(s):
            return int(s.fillna(False).sum())
        return int(
            s.astype(str).str.strip().str.lower()
            .isin(["true", "t", "1", "yes", "y"])
            .sum()
        )

    counts = {
        "n_rows": int(len(d)),
        "n_forward_ab_eligible": n_true("forward_ab_eligible"),
        "n_forward_ba_eligible": n_true("forward_ba_eligible"),
        "n_common4": n_true("common4"),
    }
    if min(
        counts["n_forward_ab_eligible"],
        counts["n_forward_ba_eligible"],
        counts["n_common4"],
    ) <= 0:
        raise RuntimeError(f"{cid}: one or more primary support sets are empty")

    return counts


def build_step7_compact_cache(
    transition: Path,
    compact_path: Path,
    cid: str,
    dt: int,
) -> Dict[str, object]:
    """Write a compact dt-specific Parquet containing only Step-7 fields."""
    try:
        import pyarrow.parquet as pq
        table = pq.read_table(
            transition,
            columns=STEP7_COMPACT_COLUMNS,
            filters=[("dt", "=", int(dt))],
        )
        if table.num_rows <= 0:
            raise RuntimeError(f"{cid}: no dt={dt} rows available for compact Step-7 cache")
        compact_path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, compact_path, compression="zstd")
    except Exception:
        d = pd.read_parquet(transition, columns=STEP7_COMPACT_COLUMNS)
        d = d[pd.to_numeric(d["dt"], errors="coerce").eq(int(dt))].copy()
        if d.empty:
            raise RuntimeError(f"{cid}: no dt={dt} rows available for compact Step-7 cache")
        compact_path.parent.mkdir(parents=True, exist_ok=True)
        d.to_parquet(compact_path, index=False, compression="zstd")
    return validate_step7_compact_cache(compact_path, cid, dt)


def validate_step7_compact_cache(path: Path, cid: str, dt: int) -> Dict[str, object]:
    cols = set(transition_schema(path))
    missing = sorted(set(STEP7_COMPACT_COLUMNS) - cols)
    if missing:
        raise ValueError(f"{cid}: compact Step-7 cache missing columns {missing}")
    d = pd.read_parquet(
        path,
        columns=[
            "configuration_id", "subject", "dt",
            "forward_ab_eligible", "forward_ba_eligible", "common4",
        ],
    )
    if d.empty:
        raise ValueError(f"{cid}: compact Step-7 cache is empty")
    configs = set(d["configuration_id"].astype(str).dropna().unique())
    if configs != {cid}:
        raise ValueError(f"{cid}: compact cache configuration_id values are {sorted(configs)}")
    subjects = pd.to_numeric(d["subject"], errors="coerce").dropna().astype(int).unique()
    if set(subjects) != {1}:
        raise ValueError(f"{cid}: compact cache expected biological subject 1, found {subjects}")
    dts = sorted(pd.to_numeric(d["dt"], errors="coerce").dropna().astype(int).unique())
    if dts != [int(dt)]:
        raise ValueError(f"{cid}: compact cache expected dt={dt}, found {dts}")

    def n_true(column: str) -> int:
        x = d[column]
        if pd.api.types.is_bool_dtype(x):
            return int(x.fillna(False).sum())
        return int(x.astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes", "y"]).sum())

    out = {
        "n_rows": int(len(d)),
        "n_forward_ab_eligible": n_true("forward_ab_eligible"),
        "n_forward_ba_eligible": n_true("forward_ba_eligible"),
        "n_common4": n_true("common4"),
    }
    if min(out["n_forward_ab_eligible"], out["n_forward_ba_eligible"], out["n_common4"]) <= 0:
        raise RuntimeError(f"{cid}: compact Step-7 cache has an empty primary support set")
    return out


def cleanup_heavy_configuration_intermediates(cdir: Path) -> None:
    """Delete heavy configuration-local intermediates after compact cache validation."""
    for name in [
        "2-clonotype_state_inference",
        "4-multirepresentation_trajectory_assembly",
        "5-longitudinal_transition_assembly",
    ]:
        shutil.rmtree(cdir / name, ignore_errors=True)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build the multirepresentation Step-4/5 pseudo-1x12 ensemble from "
            "the 12-measurement design manifest and a 66-pair Step-2 bank."
        ),
    )
    p.add_argument("pseudo_root", type=Path)
    p.add_argument("--code-dir", type=Path, required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--expected-configs", type=int, default=2000)
    p.add_argument("--max-configs", type=int, default=None)
    p.add_argument("--n-jobs", type=int, default=0)
    p.add_argument(
        "--n-jobs-profile",
        choices=["memory", "balanced", "speed"],
        default="balanced",
    )
    p.add_argument("--restart-pairbank", action="store_true")
    p.add_argument("--restart-configs", action="store_true")
    p.add_argument(
        "--posterior-step5",
        action="store_true",
        help="Optional full Step-5 posterior displacement propagation.",
    )
    p.add_argument("--posterior-samples", type=int, default=1000)
    p.add_argument("--posterior-seed", type=int, default=123)
    p.add_argument(
        "--compact-step7-cache",
        action="store_true",
        help="Write a compact dt-specific Step-7 cache after each completed configuration.",
    )
    p.add_argument(
        "--compact-step7-dt", type=int, default=1,
        help="Temporal lag retained in the compact Step-7 cache.",
    )
    p.add_argument(
        "--cleanup-after-compact",
        action="store_true",
        help="Delete heavy Step-2/4/5 intermediates after compact cache validation.",
    )
    return p


def main() -> int:
    args = parser().parse_args()

    if args.expected_configs < 0:
        raise ValueError("--expected-configs must be >=0")
    if args.max_configs is not None and args.max_configs < 1:
        raise ValueError("--max-configs must be >=1")
    if not (0.0 < float(args.alpha) < 1.0):
        raise ValueError("--alpha must be between 0 and 1")
    if int(args.compact_step7_dt) < 1:
        raise ValueError("--compact-step7-dt must be >=1")
    if args.cleanup_after_compact and not args.compact_step7_cache:
        raise ValueError("--cleanup-after-compact requires --compact-step7-cache")

    root = args.pseudo_root.expanduser().resolve(strict=True)
    code = args.code_dir.expanduser().resolve(strict=True)

    step2 = code / STEP2_NAME
    step4 = code / STEP4_NAME
    step5 = code / STEP5_NAME
    noise = code / NOISE_NAME
    for path in (step2, step4, step5, noise):
        if not path.is_file():
            raise FileNotFoundError(path)

    source = load_source_manifest(root)
    configs = load_configuration_manifest(root, int(args.expected_configs))

    bank_results, bank_pair_dir, pair_manifest = ensure_pairbank(
        root=root,
        code=code,
        step2=step2,
        noise=noise,
        alpha=float(args.alpha),
        n_jobs=int(args.n_jobs),
        n_jobs_profile=str(args.n_jobs_profile),
        restart=bool(args.restart_pairbank),
    )

    lookup = {
        pair_key(row.source_rep1, row.source_rep2): int(row.pair_index)
        for row in pair_manifest.itertuples(index=False)
    }

    config_ids = list(dict.fromkeys(configs["configuration_id"].astype(str).tolist()))
    if args.max_configs is not None:
        config_ids = config_ids[: int(args.max_configs)]

    ensemble_root = root / ENSEMBLE_DIRNAME
    ensemble_root.mkdir(parents=True, exist_ok=True)

    status_rows: List[Dict[str, object]] = []

    for index, cid in enumerate(config_ids, start=1):
        print("\n" + "=" * 84)
        print(f"Configuration {cid} ({index}/{len(config_ids)})")
        print("=" * 84)

        cdir = ensemble_root / cid
        step2_dir = cdir / "2-clonotype_state_inference"
        pair_dir = step2_dir / "per_clone_latent_logP"
        aggregate = step2_dir / "per_clone_latent_subject.parquet"

        step4_dir = cdir / "4-multirepresentation_trajectory_assembly"
        trajectory = step4_dir / "multirepresentation_trajectories_long.parquet"

        step5_dir = cdir / "5-longitudinal_transition_assembly"
        transition = step5_dir / "longitudinal_transitions.parquet"
        compact_transition = cdir / STEP7_COMPACT_REL

        if args.restart_configs and cdir.exists():
            shutil.rmtree(cdir)

        if args.compact_step7_cache and compact_transition.is_file():
            audit = validate_step7_compact_cache(compact_transition, cid, int(args.compact_step7_dt))
            print("[resume] compact Step-7 cache complete")
            status_rows.append({
                "configuration_id": cid, "status": "complete",
                "step2": "", "trajectory": "",
                "transitions": str(compact_transition),
                "storage_mode": "step7_compact", **audit,
            })
            continue

        if transition.is_file() and trajectory.is_file() and aggregate.is_file():
            audit = validate_transition_output(transition, cid)
            print("[resume] configuration complete and Step-7 ready")
            if args.compact_step7_cache:
                audit = build_step7_compact_cache(transition, compact_transition, cid, int(args.compact_step7_dt))
                if args.cleanup_after_compact:
                    cleanup_heavy_configuration_intermediates(cdir)
                transition_for_status = compact_transition
                storage_mode = "step7_compact"
            else:
                transition_for_status = transition
                storage_mode = "full"
            status_rows.append({
                "configuration_id": cid, "status": "complete",
                "step2": str(aggregate) if not args.cleanup_after_compact else "",
                "trajectory": str(trajectory) if not args.cleanup_after_compact else "",
                "transitions": str(transition_for_status),
                "storage_mode": storage_mode, **audit,
            })
            continue

        rows = configs[configs["configuration_id"] == cid].copy()
        if len(rows) != 12:
            raise ValueError(f"{cid}: expected 12 configuration rows")

        pair_dir.mkdir(parents=True, exist_ok=True)
        for old in pair_dir.iterdir():
            if old.is_file() or old.is_symlink():
                old.unlink()

        map_rows: List[Dict[str, object]] = []

        for t in range(1, 7):
            time_rows = (
                rows[rows["pseudo_time"] == t]
                .sort_values("replicate")
                .reset_index(drop=True)
            )
            if len(time_rows) != 2 or set(time_rows["replicate"]) != {1, 2}:
                raise ValueError(f"{cid} T{t}: invalid replicate assignment")

            a = str(time_rows.loc[0, "source_id"])
            b = str(time_rows.loc[1, "source_id"])
            canonical = pair_key(a, b)
            pair_index = lookup[canonical]

            source_pair = find_pair_table(bank_pair_dir, pair_index)
            frame = pd.read_parquet(source_pair)
            orientation_swapped = (a, b) != canonical
            if orientation_swapped:
                frame = swap_rep_orientation(frame)

            frame["subject"] = 1
            frame["time"] = int(t)
            frame["configuration_id"] = cid
            frame["source_rep1_id"] = a
            frame["source_rep2_id"] = b
            frame["pair_bank_index"] = int(pair_index)
            if "pair_id" in frame.columns:
                frame["pair_id"] = f"1_{t}"

            out_pair = pair_dir / f"latent_per_clone_logP_1_{t}.parquet"
            frame.to_parquet(out_pair, index=False, compression="zstd")

            posterior_src = find_pair_posterior(bank_pair_dir, pair_index)
            posterior_dst = pair_dir / f"latent_per_clone_logP_1_{t}.posterior.npz"
            if args.posterior_step5:
                if posterior_src is None:
                    raise FileNotFoundError(
                        f"{cid} T{t}: pair-bank posterior NPZ not available"
                    )
                safe_symlink(posterior_src, posterior_dst)

            map_rows.append(
                {
                    "configuration_id": cid,
                    "pseudo_time": t,
                    "replicate1_source": a,
                    "replicate2_source": b,
                    "pair_bank_index": pair_index,
                    "orientation_swapped": int(orientation_swapped),
                    "pair_bank_file": str(source_pair),
                }
            )

        pd.DataFrame(map_rows).to_csv(
            step2_dir / "configuration_pair_map.tsv",
            sep="\t",
            index=False,
        )

        # Step-2 aggregate-only cache guard must describe the SIX remapped
        # pseudo-time pairs, not the complete 66-pair bank.
        local_fit_manifest = build_configuration_fit_cache(
            bank_results=bank_results,
            step2_dir=step2_dir,
            map_rows=map_rows,
        )
        print(
            "[config-fit-cache] wrote 6 remapped successful pairs:",
            local_fit_manifest,
        )

        # Recompute configuration-specific p_value / observable from these six
        # pseudo states; do not refit the pair models.
        aggregate_cmd = [
            sys.executable,
            "-u",
            str(step2),
            "--aggregate-only",
            "--results-dir",
            str(step2_dir),
            "--fit-cache-dir",
            str(step2_dir),
            "--perclone-dir",
            str(pair_dir),
            "--aggregate-results-dir",
            str(step2_dir),
            "--alpha",
            repr(float(args.alpha)),
            "--null",
            "subject",
            "--tail",
            "low",
            "--intermediate-format",
            "parquet",
            "--parquet-compression",
            "zstd",
            "--file-sep",
            "\t",
        ]
        run(
            aggregate_cmd,
            root / "logs" / cid / "step2_aggregate_only.log",
            code,
        )

        if not aggregate.is_file():
            raise RuntimeError(f"{cid}: Step 2 aggregate output missing")

        aggregate_frame = pd.read_parquet(aggregate)
        aggregate_frame["configuration_id"] = cid
        aggregate_frame.to_parquet(aggregate, index=False, compression="zstd")

        # Step 4: preserve latent consensus, observed A/B, and observation layer.
        step4_dir.mkdir(parents=True, exist_ok=True)
        cmd4 = [
            sys.executable,
            "-u",
            str(step4),
            "--per-clone",
            str(aggregate),
            "--outdir",
            str(step4_dir),
            "--reference-alpha",
            repr(float(args.alpha)),
            "--extra-group-cols",
            "configuration_id",
            "--output-format",
            "parquet",
            "--parquet-compression",
            "zstd",
            "--streaming",
        ]
        run(cmd4, root / "logs" / cid / "step4.log", code)

        if not trajectory.is_file():
            raise RuntimeError(f"{cid}: Step 4 output missing: {trajectory}")

        # Step 5: build the full transition universe; no TT restriction.
        step5_dir.mkdir(parents=True, exist_ok=True)
        cmd5 = [
            sys.executable,
            "-u",
            str(step5),
            "--trajectories",
            str(trajectory),
            "--out",
            str(transition),
            "--extra-group-cols",
            "configuration_id",
            "--min-dt",
            "1",
            "--max-dt",
            "5",
            "--parquet-compression",
            "zstd",
            "--parquet-row-group-size",
            "100000",
            "--work-dir",
            str(step5_dir / ".step5_work"),
            "--report-md",
            str(step5_dir / "transition_assembly_report.md"),
        ]

        if args.posterior_step5:
            cmd5 += [
                "--posterior-dir",
                str(pair_dir),
                "--allow-shared-posterior-files",
                "--n-posterior-samples",
                str(int(args.posterior_samples)),
                "--posterior-seed",
                str(int(args.posterior_seed)),
            ]

        run(cmd5, root / "logs" / cid / "step5.log", code)

        if not transition.is_file():
            raise RuntimeError(f"{cid}: Step 5 output missing: {transition}")

        audit = validate_transition_output(transition, cid)
        transition_for_status = transition
        storage_mode = "full"
        if args.compact_step7_cache:
            audit = build_step7_compact_cache(transition, compact_transition, cid, int(args.compact_step7_dt))
            transition_for_status = compact_transition
            storage_mode = "step7_compact"
            print(f"[compact] wrote Step-7 dt={int(args.compact_step7_dt)} cache: {compact_transition}")
            if args.cleanup_after_compact:
                cleanup_heavy_configuration_intermediates(cdir)
                print("[cleanup] removed heavy Step-2/4/5 configuration intermediates")
        status_rows.append({
            "configuration_id": cid, "status": "complete",
            "step2": str(aggregate) if not args.cleanup_after_compact else "",
            "trajectory": str(trajectory) if not args.cleanup_after_compact else "",
            "transitions": str(transition_for_status),
            "storage_mode": storage_mode, **audit,
        })

        pd.DataFrame(status_rows).to_csv(
            root / "ensemble_build_status_v3.tsv",
            sep="\t",
            index=False,
        )

    status = pd.DataFrame(status_rows)
    status.to_csv(root / "ensemble_build_status_v3.tsv", sep="\t", index=False)

    meta = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "ensemble_root": str(ensemble_root),
        "pairbank_root": str(root / PAIRBANK_DIRNAME),
        "source_biological_subjects": 1,
        "source_measurements": 12,
        "unique_pair_fits": 66,
        "configurations_processed": int(len(status)),
        "alpha_reference": float(args.alpha),
        "configuration_specific_p_value_and_observability": True,
        "aggregate_only_fit_cache_scope": "six remapped pseudo-time pairs per configuration",
        "pair_fit_parameters_reestimated_per_configuration": False,
        "configuration_id_is_biological_subject": False,
        "full_transition_universe_retained": True,
        "primary_forward_support_fields": [
            "forward_ab_eligible",
            "forward_ba_eligible",
        ],
        "primary_fluctuation_support_field": "common4",
        "posterior_step5": bool(args.posterior_step5),
        "step7_compact_cache": bool(args.compact_step7_cache),
        "step7_compact_dt": int(args.compact_step7_dt),
        "cleanup_after_compact": bool(args.cleanup_after_compact),
    }
    (root / "00_ensemble_build_v3.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("\n[DONE]")
    print("Biological source subjects      : 1")
    print("Source technical measurements   : 12")
    print("Unique Step-2 pair fits         : 66")
    print("Configurations processed        :", len(status))
    print("Configuration as biological N   : NO")
    print("Full transition universe        : YES")
    print("Ensemble root                   :", ensemble_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
