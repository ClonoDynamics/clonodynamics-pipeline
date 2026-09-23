#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-pseudo_fluctuation_conditioning_validation_pairbank.py
========================================================

FAST ClonoDynamics Step 8 for the pseudo-1x12 technical-null ensemble.

This analyzer consumes the reusable transition bank produced by:

    pseudo_1x12_build_step8_transition_bank.py

It does NOT consume 2,000 full Step-5 transition tables and does NOT require
configuration-local Step-2/4/5 reconstruction.

Statistical design
------------------
- 1 biological source subject.
- 12 independent technical RepSeq measurements.
- 2,000 randomized configurations.
- configuration_id is a randomization realization, never biological N.
- randomization Q2.5-Q97.5 intervals are not biological confidence intervals.
- no subject bootstrap.
- pseudo-lag is arbitrary technical ordering, not biological duration.

Primary support
---------------
common4 only:
    both technical replicates positive at both endpoints.

Primary fluctuation decomposition
---------------------------------
Within configuration / pseudo-lag / abundance bin:

    cross_cov
        = Cov(dx_observed_rep1, dx_observed_rep2)

    same_var_mean
        = 0.5 * [Var(dx_observed_rep1) + Var(dx_observed_rep2)]

    replicate_specific_excess
        = same_var_mean - cross_cov

cross_cov remains signed and is never clipped at zero.

Primary conditioning
--------------------
    xmid_latent

Sensitivities
-------------
    xstar_latent
    xmid_observed

All sensitivities reuse the exact primary abundance grid.

Complete-case policy
--------------------
- pooled covariance: finite paired observed displacements only;
- each conditioning coordinate: its own complete cases;
- direct conditioning-sensitivity contrasts: pairwise-matched complete cases.

Technical optimization
----------------------
The transition bank contains each unique pair-vs-pair common4 block once.
For a configuration, the analyzer selects its 15 blocks and applies the
appropriate endpoint-order sign to dxA/dxB.  Reversal changes the sign of both
displacements but not covariance, variance, xmid, xstar or xmid_observed.

Outputs intentionally match the historical Step-8 publication contract:

    01_configuration_dt_metrics.csv
    02_ensemble_dt_summary.csv
    03_configuration_binned_metrics_long.csv
    04_binned_randomization_envelope.csv
    05_configuration_temporal_slopes.csv
    06_temporal_slope_randomization_summary.csv
    07_interval_metrics_by_configuration.csv
    08_conditioning_sensitivity_by_configuration.csv
    09_support_by_configuration_dt.csv
    10_global_bin_edges.csv
    11_configuration_curve_descriptors.csv
    12_configuration_audit.csv

Typical
-------
python code/dynamics_4/8-pseudo_fluctuation_conditioning_validation.py \
  --bank-root ./results_4/pseudo/step8_transition_bank_v1 \
  --outdir ./results_4/8-pseudo_fluctuation_conditioning_validation \
  --expected-configs 2000 \
  --dt 1 2 3 4 5 \
  --n-bins 20 \
  --min-n 10 \
  --min-valid-bins 6 \
  --min-global-n 20 \
  --seed 123

Testing:
  add --max-configs 3
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v4-fast-pair-transition-bank-2026-09-17"
PRIMARY_CONDITIONING = "xmid_latent"
CONDITIONINGS = ("xmid_latent", "xstar_latent", "xmid_observed")

PRIMARY_METRICS = (
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
    "shared_fraction",
)

GLOBAL_SLOPE_METRICS = (
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
    "equal_interval_cross_cov",
    "equal_interval_same_var_mean",
    "equal_interval_replicate_specific_excess",
)

BLOCK_COLUMNS = [
    "xmid_latent",
    "xstar_latent",
    "xmid_observed",
    "dx_observed_rep1",
    "dx_observed_rep2",
]


# -----------------------------------------------------------------------------
# Utilities
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


def stable_hash(payload: Dict[str, object]) -> str:
    raw = json.dumps(
        json_safe(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def path_identity(path: Path) -> Dict[str, object]:
    p = path.expanduser().resolve(strict=True)
    st = p.stat()
    return {
        "path": str(p),
        "size_bytes": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns),
    }


def finite_numeric(values) -> np.ndarray:
    x = pd.to_numeric(
        pd.Series(values), errors="coerce"
    ).to_numpy(dtype=float)
    return x[np.isfinite(x)]


def safe_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(ok.sum()) < 2:
        return np.nan
    x = x[ok]
    y = y[ok]
    if np.ptp(x) <= 0:
        return np.nan
    return float(np.polyfit(x, y, 1)[0])


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(ok.sum()) < 2:
        return np.nan
    x = x[ok]
    y = y[ok]
    if np.std(x) <= 0 or np.std(y) <= 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(ok.sum()) < 3:
        return np.nan
    xr = pd.Series(x[ok]).rank(method="average").to_numpy(float)
    yr = pd.Series(y[ok]).rank(method="average").to_numpy(float)
    if np.std(xr) <= 0 or np.std(yr) <= 0:
        return np.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def sample_values(
    values: np.ndarray,
    n: int,
    seed: int,
) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size <= int(n):
        return x
    rng = np.random.default_rng(int(seed))
    idx = rng.choice(x.size, size=int(n), replace=False)
    return x[idx]


# -----------------------------------------------------------------------------
# Bank loading
# -----------------------------------------------------------------------------

def load_bank_manifests(
    bank_root: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    interval_path = bank_root / "00_configuration_interval_map.csv"
    block_path = bank_root / "00_transition_block_manifest.csv"
    support_path = bank_root / "00_configuration_support_by_dt.csv"

    for path in [interval_path, block_path, support_path]:
        if not path.is_file():
            raise FileNotFoundError(path)

    intervals = pd.read_csv(interval_path)
    blocks = pd.read_csv(block_path)
    support = pd.read_csv(support_path)

    interval_required = {
        "configuration_id",
        "t0",
        "t1",
        "dt",
        "block_id",
        "reverse_sign",
    }
    block_required = {
        "block_id",
        "block_path",
        "n_transitions",
        "n_common4",
    }
    support_required = {
        "configuration_id",
        "dt",
        "n_transitions",
        "n_common4",
        "fraction_common4",
    }

    for frame, required, label in [
        (intervals, interval_required, "interval map"),
        (blocks, block_required, "block manifest"),
        (support, support_required, "support table"),
    ]:
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(
                f"{label}: missing required columns {missing}"
            )

    if pd.api.types.is_bool_dtype(intervals["reverse_sign"]):
        intervals["reverse_sign"] = intervals["reverse_sign"].fillna(False)
    elif pd.api.types.is_numeric_dtype(intervals["reverse_sign"]):
        intervals["reverse_sign"] = (
            pd.to_numeric(intervals["reverse_sign"], errors="coerce")
            .fillna(0).ne(0)
        )
    else:
        intervals["reverse_sign"] = (
            intervals["reverse_sign"].astype(str)
            .str.strip().str.lower()
            .isin(["true", "t", "1", "yes", "y"])
        )

    intervals["configuration_id"] = intervals[
        "configuration_id"
    ].astype(str)
    support["configuration_id"] = support[
        "configuration_id"
    ].astype(str)
    blocks["block_id"] = blocks["block_id"].astype(str)

    return intervals, blocks, support


def select_configurations(
    intervals: pd.DataFrame,
    expected_configs: int,
    max_configs: Optional[int],
) -> List[str]:
    ids = list(
        dict.fromkeys(
            intervals["configuration_id"].astype(str).tolist()
        )
    )

    if int(expected_configs) > 0:
        expected = [
            f"C{i:06d}"
            for i in range(1, int(expected_configs) + 1)
        ]
        if ids != expected:
            raise ValueError(
                "Step-8 transition bank configuration IDs do not match "
                f"C000001..C{int(expected_configs):06d}"
            )

    if max_configs is not None:
        ids = ids[: int(max_configs)]
    return ids


def load_blocks_in_memory(
    block_manifest: pd.DataFrame,
    interval_map: pd.DataFrame,
    selected_ids: Sequence[str],
) -> Dict[str, Dict[str, np.ndarray]]:
    used = set(
        interval_map.loc[
            interval_map["configuration_id"].isin(selected_ids),
            "block_id",
        ].astype(str)
    )

    manifest = block_manifest[
        block_manifest["block_id"].astype(str).isin(used)
    ].copy()

    if manifest["block_id"].nunique() != len(used):
        missing = sorted(
            used - set(manifest["block_id"].astype(str))
        )
        raise ValueError(
            f"Block manifest missing required blocks: {missing[:10]}"
        )

    cache: Dict[str, Dict[str, np.ndarray]] = {}

    for i, row in enumerate(
        manifest.sort_values("block_id").itertuples(index=False),
        start=1,
    ):
        bid = str(row.block_id)
        path = Path(str(row.block_path)).expanduser().resolve(strict=True)
        d = pd.read_parquet(path, columns=BLOCK_COLUMNS)

        arrays = {
            c: pd.to_numeric(
                d[c], errors="coerce"
            ).to_numpy(dtype=float)
            for c in BLOCK_COLUMNS
        }
        cache[bid] = arrays

        if i % 100 == 0 or i == len(manifest):
            print(
                f"[memory] loaded block {i}/{len(manifest)}"
            )

    return cache


# -----------------------------------------------------------------------------
# Configuration reconstruction
# -----------------------------------------------------------------------------

def reconstruct_configuration_frame(
    cid: str,
    interval_map: pd.DataFrame,
    block_cache: Dict[str, Dict[str, np.ndarray]],
    dt_values: Sequence[int],
) -> pd.DataFrame:
    rows = interval_map[
        interval_map["configuration_id"].astype(str).eq(cid)
        & pd.to_numeric(
            interval_map["dt"], errors="coerce"
        ).isin([int(v) for v in dt_values])
    ].copy()

    if rows.empty:
        raise ValueError(
            f"{cid}: no interval-map rows for requested lags"
        )

    parts: List[pd.DataFrame] = []

    for row in rows.itertuples(index=False):
        bid = str(row.block_id)
        if bid not in block_cache:
            raise KeyError(
                f"{cid}: block {bid} is not loaded"
            )

        block = block_cache[bid]
        n = len(block["dx_observed_rep1"])
        sign = -1.0 if bool(row.reverse_sign) else 1.0

        parts.append(
            pd.DataFrame(
                {
                    "t0": np.full(n, int(row.t0), dtype=np.int16),
                    "t1": np.full(n, int(row.t1), dtype=np.int16),
                    "dt": np.full(n, int(row.dt), dtype=np.int8),
                    "xmid_latent": block["xmid_latent"],
                    "xstar_latent": block["xstar_latent"],
                    "xmid_observed": block["xmid_observed"],
                    "dx_observed_rep1":
                        sign * block["dx_observed_rep1"],
                    "dx_observed_rep2":
                        sign * block["dx_observed_rep2"],
                }
            )
        )

    return pd.concat(parts, ignore_index=True)


# -----------------------------------------------------------------------------
# Global abundance grid
# -----------------------------------------------------------------------------

def derive_global_edges(
    selected_ids: Sequence[str],
    interval_map: pd.DataFrame,
    block_cache: Dict[str, Dict[str, np.ndarray]],
    args: argparse.Namespace,
) -> Tuple[np.ndarray, Dict[str, object]]:
    parts: List[np.ndarray] = []

    for i, cid in enumerate(selected_ids):
        rows = interval_map[
            interval_map["configuration_id"].astype(str).eq(str(cid))
            & pd.to_numeric(
                interval_map["dt"], errors="coerce"
            ).isin([int(v) for v in args.dt])
        ]

        values = []
        for row in rows.itertuples(index=False):
            arr = block_cache[str(row.block_id)][PRIMARY_CONDITIONING]
            arr = arr[np.isfinite(arr)]
            if arr.size:
                values.append(arr)

        if not values:
            continue

        x = np.concatenate(values)
        parts.append(
            sample_values(
                x,
                int(args.edge_sample_per_config),
                int(args.seed) + i,
            )
        )

    if not parts:
        raise ValueError(
            "No finite xmid_latent values available for global Step-8 bins"
        )

    pooled = np.concatenate(parts)
    x_min = float(
        np.quantile(pooled, float(args.range_q_low))
    )
    x_max = float(
        np.quantile(pooled, float(args.range_q_high))
    )

    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min >= x_max:
        raise ValueError(
            "Invalid Step-8 abundance range"
        )

    edges = np.linspace(
        x_min,
        x_max,
        int(args.n_bins) + 1,
    )

    return edges, {
        "edge_source":
            "equal-size sample per selected configuration from bank xmid_latent",
        "edge_sample_per_config": int(args.edge_sample_per_config),
        "primary_coordinate": PRIMARY_CONDITIONING,
        "sensitivity_coordinates_reuse_primary_edges": True,
        "candidate_coordinates": list(CONDITIONINGS),
        "q_low": float(args.range_q_low),
        "q_high": float(args.range_q_high),
        "n_bins": int(args.n_bins),
        "x_min": float(x_min),
        "x_max": float(x_max),
    }


def write_edges(path: Path, edges: np.ndarray) -> None:
    pd.DataFrame(
        {
            "bin_id": np.arange(1, len(edges), dtype=int),
            "x_lo": edges[:-1],
            "x_hi": edges[1:],
            "x_center": 0.5 * (edges[:-1] + edges[1:]),
        }
    ).to_csv(path, index=False)


def assign_bins(
    values: np.ndarray,
    edges: np.ndarray,
) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    out = np.zeros(len(x), dtype=np.int32)

    finite = np.isfinite(x)
    if not finite.any():
        return out

    xx = x[finite]
    idx = np.searchsorted(edges, xx, side="right") - 1
    idx[xx == edges[-1]] = len(edges) - 2
    valid = (idx >= 0) & (idx < len(edges) - 1)

    pos = np.where(finite)[0]
    out[pos[valid]] = idx[valid] + 1
    return out


# -----------------------------------------------------------------------------
# Covariance decomposition
# -----------------------------------------------------------------------------

def covariance_decomposition(
    dx_a: np.ndarray,
    dx_b: np.ndarray,
    ddof: int,
) -> Dict[str, float]:
    a = np.asarray(dx_a, dtype=float)
    b = np.asarray(dx_b, dtype=float)
    valid = np.isfinite(a) & np.isfinite(b)
    a = a[valid]
    b = b[valid]

    n = int(len(a))
    if n <= int(ddof):
        return {
            "n": n,
            "mean_dx_A": np.nan,
            "mean_dx_B": np.nan,
            "var_dx_A": np.nan,
            "var_dx_B": np.nan,
            "cross_cov": np.nan,
            "same_var_mean": np.nan,
            "replicate_specific_excess": np.nan,
            "shared_fraction": np.nan,
            "identity_rhs": np.nan,
            "identity_residual": np.nan,
        }

    mean_a = float(np.mean(a))
    mean_b = float(np.mean(b))
    var_a = float(np.var(a, ddof=int(ddof)))
    var_b = float(np.var(b, ddof=int(ddof)))

    denom = float(n - int(ddof))
    cov = float(
        np.sum((a - mean_a) * (b - mean_b))
        / denom
    )

    same = 0.5 * (var_a + var_b)
    excess = same - cov
    shared_fraction = (
        cov / same
        if np.isfinite(same) and same > 0
        else np.nan
    )

    avg = 0.5 * (a + b)
    diff = a - b
    rhs = float(
        np.var(avg, ddof=int(ddof))
        - 0.25 * np.var(diff, ddof=int(ddof))
    )

    return {
        "n": n,
        "mean_dx_A": mean_a,
        "mean_dx_B": mean_b,
        "var_dx_A": var_a,
        "var_dx_B": var_b,
        "cross_cov": cov,
        "same_var_mean": same,
        "replicate_specific_excess": excess,
        "shared_fraction": shared_fraction,
        "identity_rhs": rhs,
        "identity_residual": cov - rhs,
    }


def mask_under_supported_metrics(
    metrics: Dict[str, float],
    min_n: int,
) -> Tuple[Dict[str, float], bool]:
    valid = int(metrics.get("n", 0)) >= int(min_n)
    if valid:
        return dict(metrics), True

    out = dict(metrics)
    for key in [
        "mean_dx_A",
        "mean_dx_B",
        "var_dx_A",
        "var_dx_B",
        "cross_cov",
        "same_var_mean",
        "replicate_specific_excess",
        "shared_fraction",
        "identity_rhs",
        "identity_residual",
    ]:
        out[key] = np.nan
    return out, False


# -----------------------------------------------------------------------------
# Configuration metrics
# -----------------------------------------------------------------------------

def physical_interval_metrics(
    frame: pd.DataFrame,
    ddof: int,
    min_global_n: int,
) -> pd.DataFrame:
    rows = []

    for (dt, t0, t1), g in frame.groupby(
        ["dt", "t0", "t1"],
        sort=True,
    ):
        raw = covariance_decomposition(
            g["dx_observed_rep1"].to_numpy(float),
            g["dx_observed_rep2"].to_numpy(float),
            ddof=ddof,
        )
        metrics, valid = mask_under_supported_metrics(
            raw,
            min_n=min_global_n,
        )
        rows.append(
            {
                "dt": int(dt),
                "t0": int(t0),
                "t1": int(t1),
                **metrics,
                "meets_min_global_n": bool(valid),
            }
        )

    return pd.DataFrame(rows)


def configuration_dt_metrics(
    frame: pd.DataFrame,
    interval_metrics: pd.DataFrame,
    ddof: int,
    min_global_n: int,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    rows = []

    for dt in dt_values:
        g = frame[
            pd.to_numeric(frame["dt"], errors="coerce").eq(int(dt))
        ]

        raw = covariance_decomposition(
            g["dx_observed_rep1"].to_numpy(float),
            g["dx_observed_rep2"].to_numpy(float),
            ddof=ddof,
        )
        metrics, valid = mask_under_supported_metrics(
            raw,
            min_n=min_global_n,
        )

        ig = interval_metrics[
            pd.to_numeric(
                interval_metrics["dt"], errors="coerce"
            ).eq(int(dt))
            & interval_metrics[
                "meets_min_global_n"
            ].astype(bool)
        ]

        def mean_finite(column: str) -> float:
            if ig.empty:
                return np.nan
            x = finite_numeric(ig[column])
            return float(np.mean(x)) if x.size else np.nan

        rows.append(
            {
                "dt": int(dt),
                **metrics,
                "meets_min_global_n": bool(valid),
                "n_intervals_total": int(
                    frame.loc[
                        pd.to_numeric(
                            frame["dt"], errors="coerce"
                        ).eq(int(dt)),
                        ["t0", "t1"],
                    ]
                    .drop_duplicates()
                    .shape[0]
                ),
                "n_intervals_valid": int(len(ig)),
                "equal_interval_cross_cov":
                    mean_finite("cross_cov"),
                "equal_interval_same_var_mean":
                    mean_finite("same_var_mean"),
                "equal_interval_replicate_specific_excess":
                    mean_finite("replicate_specific_excess"),
            }
        )

    return pd.DataFrame(rows)


def binned_metrics_for_conditioning(
    frame: pd.DataFrame,
    conditioning: str,
    edges: np.ndarray,
    ddof: int,
    min_n: int,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    d = frame.copy()

    x = pd.to_numeric(
        d[conditioning], errors="coerce"
    ).to_numpy(float)
    a = pd.to_numeric(
        d["dx_observed_rep1"], errors="coerce"
    ).to_numpy(float)
    b = pd.to_numeric(
        d["dx_observed_rep2"], errors="coerce"
    ).to_numpy(float)

    finite = np.isfinite(x) & np.isfinite(a) & np.isfinite(b)
    d = d.loc[finite].copy()
    x = x[finite]

    d["bin_id"] = assign_bins(x, edges)
    d = d[d["bin_id"] > 0].copy()

    rows = []

    for dt in dt_values:
        dt_frame = d[
            pd.to_numeric(d["dt"], errors="coerce").eq(int(dt))
        ]

        for bin_id, g in dt_frame.groupby("bin_id", sort=True):
            bin_id = int(bin_id)
            metrics = covariance_decomposition(
                g["dx_observed_rep1"].to_numpy(float),
                g["dx_observed_rep2"].to_numpy(float),
                ddof=ddof,
            )

            rows.append(
                {
                    "conditioning": conditioning,
                    "dt": int(dt),
                    "bin_id": bin_id,
                    "x_lo": float(edges[bin_id - 1]),
                    "x_hi": float(edges[bin_id]),
                    "x_center": float(
                        0.5 * (edges[bin_id - 1] + edges[bin_id])
                    ),
                    "x_median": float(
                        np.median(
                            pd.to_numeric(
                                g[conditioning],
                                errors="coerce",
                            ).to_numpy(float)
                        )
                    ),
                    "n_intervals": int(
                        g[["t0", "t1"]]
                        .drop_duplicates()
                        .shape[0]
                    ),
                    **metrics,
                    "meets_min_n": bool(
                        int(metrics["n"]) >= int(min_n)
                    ),
                }
            )

    return pd.DataFrame(rows)


def abundance_curve_descriptors(
    binned: pd.DataFrame,
    min_valid_bins: int,
) -> pd.DataFrame:
    rows = []

    for (conditioning, dt), g in binned.groupby(
        ["conditioning", "dt"],
        sort=True,
    ):
        valid = g[
            g["meets_min_n"].astype(bool)
            & np.isfinite(
                pd.to_numeric(
                    g["cross_cov"], errors="coerce"
                ).to_numpy(float)
            )
        ].copy()

        row = {
            "conditioning": str(conditioning),
            "dt": int(dt),
            "n_bins_valid": int(len(valid)),
            "meets_min_valid_bins": bool(
                len(valid) >= int(min_valid_bins)
            ),
        }

        if len(valid) >= int(min_valid_bins):
            x = valid["x_center"].to_numpy(float)

            for metric in [
                "cross_cov",
                "same_var_mean",
                "replicate_specific_excess",
            ]:
                y = valid[metric].to_numpy(float)
                row[f"mean_{metric}_across_bins"] = float(
                    np.mean(y)
                )
                row[f"median_{metric}_across_bins"] = float(
                    np.median(y)
                )
                row[f"slope_{metric}_vs_x"] = safe_slope(x, y)

            y = valid["cross_cov"].to_numpy(float)
            row["mean_abs_cross_cov_across_bins"] = float(
                np.mean(np.abs(y))
            )
            row["fraction_bins_cross_cov_positive"] = float(
                np.mean(y > 0)
            )
        else:
            for key in [
                "mean_cross_cov_across_bins",
                "median_cross_cov_across_bins",
                "slope_cross_cov_vs_x",
                "mean_same_var_mean_across_bins",
                "median_same_var_mean_across_bins",
                "slope_same_var_mean_vs_x",
                "mean_replicate_specific_excess_across_bins",
                "median_replicate_specific_excess_across_bins",
                "slope_replicate_specific_excess_vs_x",
                "mean_abs_cross_cov_across_bins",
                "fraction_bins_cross_cov_positive",
            ]:
                row[key] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def conditioning_sensitivity(
    binned: pd.DataFrame,
) -> pd.DataFrame:
    primary = binned[
        binned["conditioning"].astype(str).eq(PRIMARY_CONDITIONING)
        & binned["meets_min_n"].astype(bool)
    ].copy()

    rows = []

    for sensitivity in [
        c
        for c in binned["conditioning"].astype(str).unique()
        if c != PRIMARY_CONDITIONING
    ]:
        sens = binned[
            binned["conditioning"].astype(str).eq(sensitivity)
            & binned["meets_min_n"].astype(bool)
        ]

        for dt in sorted(
            set(primary["dt"].astype(int))
            | set(sens["dt"].astype(int))
        ):
            a = primary[
                pd.to_numeric(
                    primary["dt"], errors="coerce"
                ).eq(int(dt))
            ][["bin_id", "cross_cov"]].rename(
                columns={"cross_cov": "primary"}
            )
            b = sens[
                pd.to_numeric(
                    sens["dt"], errors="coerce"
                ).eq(int(dt))
            ][["bin_id", "cross_cov"]].rename(
                columns={"cross_cov": "sensitivity"}
            )

            m = a.merge(b, on="bin_id", how="inner").dropna()
            if m.empty:
                rows.append(
                    {
                        "sensitivity_conditioning": sensitivity,
                        "dt": int(dt),
                        "n_matched_bins": 0,
                        "pearson_r": np.nan,
                        "spearman_r": np.nan,
                        "median_abs_difference": np.nan,
                        "mean_abs_difference": np.nan,
                        "sign_concordance": np.nan,
                    }
                )
                continue

            x = m["primary"].to_numpy(float)
            y = m["sensitivity"].to_numpy(float)
            delta = x - y

            rows.append(
                {
                    "sensitivity_conditioning": sensitivity,
                    "dt": int(dt),
                    "n_matched_bins": int(len(m)),
                    "pearson_r": safe_pearson(x, y),
                    "spearman_r": safe_spearman(x, y),
                    "median_abs_difference": float(
                        np.median(np.abs(delta))
                    ),
                    "mean_abs_difference": float(
                        np.mean(np.abs(delta))
                    ),
                    "sign_concordance": float(
                        np.mean(np.sign(x) == np.sign(y))
                    ),
                }
            )

    return pd.DataFrame(rows)


def configuration_temporal_slopes(
    dt_metrics: pd.DataFrame,
    min_dt_points: int,
) -> Dict[str, object]:
    out: Dict[str, object] = {}

    for metric in GLOBAL_SLOPE_METRICS:
        if metric not in dt_metrics.columns:
            continue

        x = pd.to_numeric(
            dt_metrics["dt"], errors="coerce"
        ).to_numpy(float)
        y = pd.to_numeric(
            dt_metrics[metric], errors="coerce"
        ).to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(y)

        out[f"{metric}_n_dt"] = int(valid.sum())
        out[f"{metric}_slope_per_pseudolag"] = (
            safe_slope(x[valid], y[valid])
            if int(valid.sum()) >= int(min_dt_points)
            else np.nan
        )

    return out


# -----------------------------------------------------------------------------
# Ensemble summaries
# -----------------------------------------------------------------------------

def summarize_numeric_by_groups(
    frame: pd.DataFrame,
    group_cols: Sequence[str],
    exclude_cols: Sequence[str] = (),
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()

    excluded = set(group_cols) | set(exclude_cols)
    numeric_cols = [
        c for c in frame.columns
        if c not in excluded
        and pd.api.types.is_numeric_dtype(frame[c])
    ]

    rows = []
    iterator = (
        frame.groupby(list(group_cols), sort=True)
        if group_cols
        else [((), frame)]
    )

    for keys, g in iterator:
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = {
            key: value
            for key, value in zip(group_cols, keys)
        }

        for column in numeric_cols:
            x = finite_numeric(g[column])
            if x.size == 0:
                continue

            if "configuration_id" in g.columns:
                valid = pd.to_numeric(
                    g[column], errors="coerce"
                ).notna()
                n_config = int(
                    g.loc[valid, "configuration_id"].nunique()
                )
            else:
                n_config = int(x.size)

            rows.append(
                {
                    **base,
                    "metric": column,
                    "n_configurations": n_config,
                    "mean": float(np.mean(x)),
                    "median": float(np.median(x)),
                    "sd": float(np.std(x, ddof=0)),
                    "randomization_q025": float(
                        np.quantile(x, 0.025)
                    ),
                    "randomization_q975": float(
                        np.quantile(x, 0.975)
                    ),
                    "min": float(np.min(x)),
                    "max": float(np.max(x)),
                }
            )

    return pd.DataFrame(rows)


def binned_randomization_envelope(
    binned: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for keys, g in binned.groupby(
        ["conditioning", "dt", "bin_id"],
        sort=True,
    ):
        conditioning, dt, bin_id = keys
        row = {
            "conditioning": conditioning,
            "dt": int(dt),
            "bin_id": int(bin_id),
            "n_configurations": int(
                g["configuration_id"].nunique()
            ),
            "n_configurations_valid": int(
                g.loc[
                    g["meets_min_n"].astype(bool),
                    "configuration_id",
                ].nunique()
            ),
        }

        for column in [
            "x_lo",
            "x_hi",
            "x_center",
            "x_median",
            "n",
            "n_intervals",
            "mean_dx_A",
            "mean_dx_B",
            "var_dx_A",
            "var_dx_B",
            "cross_cov",
            "same_var_mean",
            "replicate_specific_excess",
            "shared_fraction",
            "identity_residual",
        ]:
            if column not in g.columns:
                continue

            source = (
                g[g["meets_min_n"].astype(bool)]
                if column in {
                    "cross_cov",
                    "same_var_mean",
                    "replicate_specific_excess",
                    "shared_fraction",
                    "identity_residual",
                    "mean_dx_A",
                    "mean_dx_B",
                    "var_dx_A",
                    "var_dx_B",
                }
                else g
            )

            x = finite_numeric(source[column])
            if x.size == 0:
                continue

            row[column] = float(np.median(x))
            row[f"{column}_randomization_q025"] = float(
                np.quantile(x, 0.025)
            )
            row[f"{column}_randomization_q975"] = float(
                np.quantile(x, 0.975)
            )

        rows.append(row)

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Checkpointing
# -----------------------------------------------------------------------------

def load_checkpoint(
    path: Path,
    signature: str,
) -> Optional[Dict[str, object]]:
    if not path.is_file():
        return None
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if obj.get("status") != "complete":
        return None
    if obj.get("analysis_signature") != signature:
        return None
    return obj


def frame_records(frame: pd.DataFrame) -> List[Dict[str, object]]:
    if frame is None or frame.empty:
        return []
    return [
        {str(k): json_safe(v) for k, v in row.items()}
        for row in frame.to_dict(orient="records")
    ]


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Fast pseudo Step-8 fluctuation validation from the reusable "
            "pair-transition bank."
        ),
    )

    p.add_argument("--bank-root", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)

    p.add_argument("--expected-configs", type=int, default=2000)
    p.add_argument("--max-configs", type=int, default=None)

    p.add_argument(
        "--dt",
        type=int,
        nargs="+",
        default=[1, 2, 3, 4, 5],
    )

    p.add_argument("--n-bins", type=int, default=20)
    p.add_argument("--range-q-low", type=float, default=0.01)
    p.add_argument("--range-q-high", type=float, default=0.99)
    p.add_argument("--edge-sample-per-config", type=int, default=2000)

    p.add_argument("--min-n", type=int, default=10)
    p.add_argument("--min-valid-bins", type=int, default=6)
    p.add_argument("--min-global-n", type=int, default=20)
    p.add_argument("--min-dt-points", type=int, default=5)

    p.add_argument(
        "--ddof",
        type=int,
        choices=[0, 1],
        default=1,
    )
    p.add_argument("--seed", type=int, default=123)

    p.add_argument("--restart", action="store_true")
    p.add_argument("--continue-on-error", action="store_true")
    return p


def validate_args(args: argparse.Namespace) -> None:
    if args.expected_configs < 0:
        raise ValueError("--expected-configs must be >=0")
    if args.max_configs is not None and args.max_configs < 1:
        raise ValueError("--max-configs must be >=1")
    if not args.dt or len(set(args.dt)) != len(args.dt):
        raise ValueError("--dt must contain unique values")
    if any(int(v) < 1 or int(v) > 5 for v in args.dt):
        raise ValueError("--dt values must be between 1 and 5")
    if args.n_bins < 3:
        raise ValueError("--n-bins must be >=3")
    if not (0 <= args.range_q_low < args.range_q_high <= 1):
        raise ValueError("Invalid quantile range")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    args = build_parser().parse_args()
    validate_args(args)

    bank_root = args.bank_root.expanduser().resolve(strict=True)
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    checkpoint_dir = outdir / "configuration_checkpoints"
    if args.restart:
        shutil.rmtree(checkpoint_dir, ignore_errors=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    intervals, block_manifest, support_all = (
        load_bank_manifests(bank_root)
    )
    selected_ids = select_configurations(
        intervals,
        int(args.expected_configs),
        args.max_configs,
    )

    intervals_selected = intervals[
        intervals["configuration_id"].isin(selected_ids)
    ].copy()

    support_selected = support_all[
        support_all["configuration_id"].isin(selected_ids)
        & pd.to_numeric(
            support_all["dt"], errors="coerce"
        ).isin([int(v) for v in args.dt])
    ].copy()

    print("[INFO] configurations selected          :", len(selected_ids))
    print("[INFO] biological source subjects       : 1")
    print("[INFO] configuration_id as biological N : NO")
    print("[INFO] Step 2/4/5 reconstruction         : NONE")
    print("[INFO] primary support                   : common4")
    print("[INFO] primary conditioning              : xmid_latent")

    block_cache = load_blocks_in_memory(
        block_manifest,
        intervals_selected,
        selected_ids,
    )

    print("[INFO] unique blocks loaded in memory    :", len(block_cache))

    edges, edge_info = derive_global_edges(
        selected_ids,
        intervals_selected,
        block_cache,
        args,
    )
    edges_path = outdir / "10_global_bin_edges.csv"
    write_edges(edges_path, edges)

    manifest = pd.DataFrame(
        {"configuration_id": selected_ids}
    )
    manifest.to_csv(
        outdir / "00_configuration_manifest.csv",
        index=False,
    )

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "bank_root": path_identity(bank_root / "00_build_config.json"),
        "interval_map_sha256":
            sha256_file(bank_root / "00_configuration_interval_map.csv"),
        "block_manifest_sha256":
            sha256_file(bank_root / "00_transition_block_manifest.csv"),
        "support_sha256":
            sha256_file(bank_root / "00_configuration_support_by_dt.csv"),
        "configuration_ids": list(selected_ids),
        "parameters": {
            "dt": [int(v) for v in args.dt],
            "n_bins": int(args.n_bins),
            "range_q_low": float(args.range_q_low),
            "range_q_high": float(args.range_q_high),
            "edge_sample_per_config": int(
                args.edge_sample_per_config
            ),
            "min_n": int(args.min_n),
            "min_valid_bins": int(args.min_valid_bins),
            "min_global_n": int(args.min_global_n),
            "min_dt_points": int(args.min_dt_points),
            "ddof": int(args.ddof),
            "seed": int(args.seed),
        },
        "edges": [float(x) for x in edges],
    }
    signature = stable_hash(signature_payload)

    (outdir / "00_run_signature.json").write_text(
        json.dumps(
            {
                "analysis_signature": signature,
                "signature_payload": signature_payload,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    statuses = []

    for i, cid in enumerate(selected_ids, start=1):
        checkpoint = checkpoint_dir / f"{cid}.json"

        previous = (
            None
            if args.restart
            else load_checkpoint(checkpoint, signature)
        )
        if previous is not None:
            print("[resume-safe]", cid)
            statuses.append(
                {
                    "configuration_id": cid,
                    "status": "complete",
                    "execution": "resumed",
                }
            )
            continue

        if checkpoint.exists():
            checkpoint.unlink()

        print(f"[config {i}/{len(selected_ids)}] {cid}")

        try:
            frame = reconstruct_configuration_frame(
                cid,
                intervals_selected,
                block_cache,
                args.dt,
            )

            interval_metrics = physical_interval_metrics(
                frame,
                ddof=int(args.ddof),
                min_global_n=int(args.min_global_n),
            )
            dt_metrics = configuration_dt_metrics(
                frame,
                interval_metrics,
                ddof=int(args.ddof),
                min_global_n=int(args.min_global_n),
                dt_values=args.dt,
            )

            binned_parts = [
                binned_metrics_for_conditioning(
                    frame,
                    conditioning=c,
                    edges=edges,
                    ddof=int(args.ddof),
                    min_n=int(args.min_n),
                    dt_values=args.dt,
                )
                for c in CONDITIONINGS
            ]
            binned = pd.concat(
                binned_parts,
                ignore_index=True,
            )

            curve_desc = abundance_curve_descriptors(
                binned,
                min_valid_bins=int(args.min_valid_bins),
            )
            sensitivity = conditioning_sensitivity(binned)
            temporal = configuration_temporal_slopes(
                dt_metrics,
                min_dt_points=int(args.min_dt_points),
            )

            support = support_selected[
                support_selected[
                    "configuration_id"
                ].astype(str).eq(cid)
            ][
                [
                    "dt",
                    "n_transitions",
                    "n_common4",
                    "fraction_common4",
                ]
            ].copy()

            residual = finite_numeric(
                binned["identity_residual"]
            )
            max_residual = (
                float(np.max(np.abs(residual)))
                if residual.size else np.nan
            )

            payload = {
                "configuration_id": cid,
                "analysis_signature": signature,
                "dt_metrics": frame_records(dt_metrics),
                "interval_metrics": frame_records(interval_metrics),
                "binned": frame_records(binned),
                "curve_descriptors": frame_records(curve_desc),
                "conditioning_sensitivity":
                    frame_records(sensitivity),
                "temporal_slopes": json_safe(temporal),
                "support": frame_records(support),
                "audit": {
                    "rows_common4_reconstructed": int(len(frame)),
                    "max_abs_covariance_identity_residual":
                        max_residual,
                },
                "status": "complete",
            }

            checkpoint.write_text(
                json.dumps(
                    payload,
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )

            statuses.append(
                {
                    "configuration_id": cid,
                    "status": "complete",
                    "execution": "computed",
                }
            )

        except Exception as exc:
            statuses.append(
                {
                    "configuration_id": cid,
                    "status": "failed",
                    "execution": "computed",
                    "message": f"{type(exc).__name__}: {exc}",
                }
            )
            pd.DataFrame(statuses).to_csv(
                outdir / "00_run_status.csv",
                index=False,
            )
            if not args.continue_on_error:
                raise
            print("[ERROR]", cid, exc, file=sys.stderr)

        pd.DataFrame(statuses).to_csv(
            outdir / "00_run_status.csv",
            index=False,
        )

    # Reassemble.
    dt_rows = []
    interval_rows = []
    binned_rows = []
    curve_rows = []
    sensitivity_rows = []
    slope_rows = []
    support_rows = []
    audit_rows = []

    n_complete = 0

    for cid in selected_ids:
        obj = load_checkpoint(
            checkpoint_dir / f"{cid}.json",
            signature,
        )
        if obj is None:
            continue

        n_complete += 1

        for key, target in [
            ("dt_metrics", dt_rows),
            ("interval_metrics", interval_rows),
            ("binned", binned_rows),
            ("curve_descriptors", curve_rows),
            ("conditioning_sensitivity", sensitivity_rows),
            ("support", support_rows),
        ]:
            for rec in obj.get(key, []):
                target.append(
                    {"configuration_id": cid, **rec}
                )

        slope_rows.append(
            {
                "configuration_id": cid,
                **obj.get("temporal_slopes", {}),
            }
        )
        audit_rows.append(
            {
                "configuration_id": cid,
                **obj.get("audit", {}),
            }
        )

    if n_complete == 0:
        raise RuntimeError(
            "No current-signature configurations completed"
        )

    if (
        not args.continue_on_error
        and n_complete != len(selected_ids)
    ):
        raise RuntimeError(
            f"Completed={n_complete}; selected={len(selected_ids)}"
        )

    dt_metrics = pd.DataFrame(dt_rows)
    interval_metrics = pd.DataFrame(interval_rows)
    binned = pd.DataFrame(binned_rows)
    curve_desc = pd.DataFrame(curve_rows)
    sensitivity = pd.DataFrame(sensitivity_rows)
    slopes = pd.DataFrame(slope_rows)
    support = pd.DataFrame(support_rows)
    audit = pd.DataFrame(audit_rows)

    dt_metrics.to_csv(
        outdir / "01_configuration_dt_metrics.csv",
        index=False,
    )

    summarize_numeric_by_groups(
        dt_metrics,
        ["dt"],
        exclude_cols=["configuration_id"],
    ).to_csv(
        outdir / "02_ensemble_dt_summary.csv",
        index=False,
    )

    binned.to_csv(
        outdir / "03_configuration_binned_metrics_long.csv",
        index=False,
    )

    binned_randomization_envelope(binned).to_csv(
        outdir / "04_binned_randomization_envelope.csv",
        index=False,
    )

    slopes.to_csv(
        outdir / "05_configuration_temporal_slopes.csv",
        index=False,
    )

    summarize_numeric_by_groups(
        slopes,
        [],
        exclude_cols=["configuration_id"],
    ).to_csv(
        outdir / "06_temporal_slope_randomization_summary.csv",
        index=False,
    )

    interval_metrics.to_csv(
        outdir / "07_interval_metrics_by_configuration.csv",
        index=False,
    )

    sensitivity.to_csv(
        outdir / "08_conditioning_sensitivity_by_configuration.csv",
        index=False,
    )

    support.to_csv(
        outdir / "09_support_by_configuration_dt.csv",
        index=False,
    )

    curve_desc.to_csv(
        outdir / "11_configuration_curve_descriptors.csv",
        index=False,
    )

    audit.to_csv(
        outdir / "12_configuration_audit.csv",
        index=False,
    )

    run_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_signature": signature,
        "bank_root": str(bank_root),
        "outdir": str(outdir),
        "n_configurations_selected": int(len(selected_ids)),
        "n_current_signature_checkpoints": int(n_complete),
        "source_biological_subjects": 1,
        "configuration_is_biological_subject": False,
        "subject_bootstrap_used": False,
        "randomization_interval_not_biological_confidence_interval": True,
        "step2_aggregate_only_rerun": False,
        "step4_rerun": False,
        "step5_rerun": False,
        "primary_support": "common4",
        "primary_conditioning": PRIMARY_CONDITIONING,
        "conditioning_sensitivities": [
            "xstar_latent",
            "xmid_observed",
        ],
        "primary_estimands": {
            "cross_cov":
                "Cov(dx_observed_rep1,dx_observed_rep2 | conditioning,dt,common4)",
            "same_var_mean":
                "0.5*[Var(dx_observed_rep1)+Var(dx_observed_rep2)]",
            "replicate_specific_excess":
                "same_var_mean-cross_cov",
            "cross_cov_truncated_at_zero": False,
        },
        "ddof": int(args.ddof),
        "dt_values": [int(v) for v in args.dt],
        "binning": edge_info,
        "min_n_per_binned_cell": int(args.min_n),
        "min_valid_bins_for_curve_descriptor":
            int(args.min_valid_bins),
        "min_global_n": int(args.min_global_n),
        "min_dt_points_for_temporal_slope":
            int(args.min_dt_points),
        "pseudo_lag_interpretation":
            "arbitrary technical ordering; not biological duration",
        "seed": int(args.seed),
    }

    (outdir / "00_run_config.json").write_text(
        json.dumps(
            json_safe(run_config),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    readme = """# Step 8 — fast pair-transition-bank pseudo fluctuation validation

Primary estimand:
    Cov(dx_observed_rep1, dx_observed_rep2)

Support:
    common4 only

Primary conditioning:
    xmid_latent

Sensitivities:
    xstar_latent
    xmid_observed

The analyzer reuses unique pair-vs-pair transition blocks and does not rebuild
Step 2 aggregate-only, Step 4, or Step 5 for each pseudo configuration.

configuration_id is a randomization index, not biological N.
Randomization percentiles are not biological confidence intervals.
"""
    (outdir / "README_outputs.md").write_text(
        readme,
        encoding="utf-8",
    )

    # Lightweight output manifest.
    rows = []
    for path in sorted(outdir.rglob("*")):
        if (
            path.is_file()
            and "configuration_checkpoints" not in path.parts
        ):
            rows.append(
                {
                    "relative_path": str(
                        path.relative_to(outdir)
                    ),
                    "suffix": path.suffix.lower(),
                    "size_bytes": int(path.stat().st_size),
                }
            )
    pd.DataFrame(rows).to_csv(
        outdir / "00_pipeline_manifest.csv",
        index=False,
    )

    print("\n[DONE] fast Step 8 pseudo fluctuation validation")
    print("Selected configurations          :", len(selected_ids))
    print("Completed configurations         :", n_complete)
    print("Biological source subjects       : 1")
    print("configuration_id as biological N : NO")
    print("Step 2/4/5 rerun                 : NO")
    print("Primary support                  : common4")
    print("Primary conditioning             : xmid_latent")
    print("Cross covariance clipped         : NO")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
