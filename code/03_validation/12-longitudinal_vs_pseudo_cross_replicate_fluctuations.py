#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12-longitudinal_vs_pseudo_cross_replicate_fluctuations.py
==========================================================

Canonical ClonoDynamics Step 12.

PURPOSE
-------
Compare the genuine longitudinal cross-replicate fluctuation covariance from
Step 11 with the empirical pseudo-longitudinal technical-null distribution from
Step 8.

The SAME fluctuation estimand is used in both datasets:

    V_cross(x, dt)
        = Cov(Delta x_A, Delta x_B | xmid_latent, dt, common4)

with

    Delta x_A = dx_observed_rep1
    Delta x_B = dx_observed_rep2.

The primary support is `common4` in both datasets. No TT/TF/FT/FF restriction
is introduced.

WHY ABSOLUTE xmid MATCHING IS PRIMARY
-------------------------------------
Unlike the forward analysis, Step 12 does NOT percentile-normalize the
conditioning coordinate.

Cross-replicate covariance is strongly abundance-dependent in magnitude, and
both Step 8 and Step 11 already condition on the same model-based absolute
coordinate, xmid_latent. Percentile alignment would compare different absolute
abundance levels between the single-donor pseudo reference and the
10-subject longitudinal cohort.

Step 12 therefore restricts every temporal lag to the shared ABSOLUTE
xmid_latent support represented in both datasets and performs no extrapolation.

The shared support is determined separately for each dt.

PSEUDO STABLE SUPPORT
---------------------
Step 8 uses one fixed xmid_latent grid across pseudo configurations and
pseudo-lags. Some configuration x bin cells do not reach the Step-8 per-cell
minimum support.

For each dt, a pseudo native bin is considered stable when

    n_configurations_valid(bin, dt)
    --------------------------------
          full Step-8 ensemble N

        >= --min-pseudo-bin-valid-fraction

(default 0.80).

The denominator is always the complete Step-8 randomization ensemble. A
configuration with no valid bin at a lag therefore remains part of the support
denominator rather than being conditionally removed.

The final matched coordinates are the stable pseudo NATIVE bin centers that
also fall inside the valid longitudinal Step-11 xmid_latent support.

CONFIGURATION-LEVEL TECHNICAL NULL
----------------------------------
For each pseudo configuration and dt, Step 12 uses the Step-8
configuration-level binned xmid_latent curve.

Only cells with Step-8 `meets_min_n = True` contribute.

Missing pseudo native bins are NEVER interpolated or imputed. In production a
pseudo configuration contributes to the empirical null for a lag only when it
has a directly valid estimate at EVERY matched stable pseudo bin
(`--min-null-bin-coverage 1.0`).

The genuine longitudinal Step-11 curve is the only curve interpolated, because
Step 11 and Step 8 use different fixed xmid_latent grids. Longitudinal
interpolation is local between adjacent valid Step-11 native bins and never
bridges an unsupported internal bin.

PRIMARY MATCHED-SUPPORT STATISTIC
---------------------------------
For each dt, the real and pseudo curves are evaluated at the SAME matched
stable Step-8 native xmid_latent bin centers.

The primary scalar statistic is

    abundance_averaged_cross_cov
        = mean_b V_cross(x_b, dt)

with equal weight over the retained matched native abundance bins.

This avoids conflating the real/pseudo comparison with different transition
counts or different abundance compositions.

The primary empirical test is one-sided:

    H_A:
        real abundance-averaged cross covariance
            >
        pseudo technical-null distribution.

Finite-sample empirical P values use the +1 correction.

SECONDARY DESCRIPTORS
---------------------
The same matched-support interpolation also retains descriptive summaries of

    same_var_mean
    replicate_specific_excess
    shared_fraction

but the pre-specified Step-12 inferential test is `cross_cov`.

LONGITUDINAL UNCERTAINTY
------------------------
Step 11 provides subject x dt x abundance-bin sufficient statistics.

Step 12 reuses these statistics to reconstruct the matched-support
abundance-averaged real cross covariance under biological subject-cluster
bootstrap resampling.

The same biological subject draw is shared across all dt, preserving the
cross-lag dependence required for the optional joint summary.

JOINT ACROSS-LAG SUMMARY
------------------------
Step 12 also reports a descriptive/inferential global statistic:

    equal-lag mean of the five matched-support abundance-averaged cross_cov
    estimates.

For the pseudo null, only configurations eligible at every requested lag enter
this joint distribution.

This joint statistic summarizes the overall excess of shared fluctuation
covariance above the technical null. It does NOT test temporal accumulation;
temporal scaling belongs to Step 13.

INPUTS
------
--step11-dir
    Final Step-11 directory containing:

        00_run_config.json
        02_fluctuation_by_bin_dt.csv
        03_global_fluctuation_by_dt.csv
        08_subject_bin_sufficient.parquet

--pseudo-step8-dir
    Final Step-8 directory containing:

        00_run_config.json
        03_configuration_binned_metrics_long.csv
        04_binned_randomization_envelope.csv
        02_ensemble_dt_summary.csv

OUTPUT
------
<outdir>/
    00_run_config.json
    00_input_manifest.csv
    00_analysis_summary.csv
    00_pipeline_manifest.csv
    README_outputs.md

    01_shared_support_by_dt.csv
    02_pointwise_cross_covariance_comparison.csv
    03_pseudo_configuration_matched_support_metrics.csv
    04_empirical_tests_by_dt.csv
    05_longitudinal_matched_support_bootstrap.csv
    06_joint_comparison.csv
    07_support_audit.csv
    07b_stable_native_bins_by_dt.csv
    08_native_global_by_dt.csv

SMOKE TEST
----------
Use:

    --max-pseudo-configs 20

to exercise the full pipeline on a small subset.

In smoke-test mode:
    - stable shared support is still defined from the full Step-8 ensemble;
    - the script selects a small UNION of configurations that collectively
      spans the final matched support at every requested lag instead of taking
      the arbitrary first N;
    - the same configuration is NOT required to span every lag because the
      primary empirical null is lag-specific;
    - at least one eligible pseudo configuration is sufficient for execution;
    - the inferential minimum remains unchanged, so empirical P values are NA
      when fewer than --min-null-configurations are eligible.

PRODUCTION RUN
--------------
python3 12-longitudinal_vs_pseudo_cross_replicate_fluctuations.py \
    --step11-dir \
    ./results_4/11-cross_replicate_fluctuation_dynamics \
    --pseudo-step8-dir \
    ./results_4/8-pseudo_fluctuation_conditioning_validation \
    --outdir \
    ./results_4/12-longitudinal_vs_pseudo_cross_replicate_fluctuations \
    --dt-values 1 2 3 4 5 \
    --min-pseudo-bin-valid-fraction 0.80 \
    --min-shared-native-bins 8 \
    --min-null-bin-coverage 1.0 \
    --min-null-configurations 30
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import polars as pl
except Exception as exc:
    raise SystemExit(
        "ERROR: polars is required. Activate the ClonoDynamics environment."
    ) from exc


SCRIPT_VERSION = "v5-direct-native-pseudo-support-2026-09-17"
DDOF = 1
SUPPORT_TOL = 1e-8

PSEUDO_METRIC_COLUMNS = [
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
    "shared_fraction",
]


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


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


def stable_sha256(payload: Mapping[str, object]) -> str:
    import hashlib
    raw = json.dumps(
        json_safe(dict(payload)),
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


def interpolate_adjacent_valid_native(
    x_native: np.ndarray,
    y_native: np.ndarray,
    x_target: np.ndarray,
    tol: float = SUPPORT_TOL,
) -> np.ndarray:
    """
    Interpolate a longitudinal curve onto target coordinates using only
    adjacent native points.

    Internal missing native values are never bridged. This is required because
    the target Step-8 grid differs from the Step-11 grid, while preserving the
    no-gap-imputation principle used for the pseudo technical null.
    """
    x = np.asarray(x_native, dtype=float)
    y = np.asarray(y_native, dtype=float)
    target = np.asarray(x_target, dtype=float)

    finite_x = np.isfinite(x)
    x = x[finite_x]
    y = y[finite_x]

    if x.size < 2:
        return np.full(target.shape, np.nan)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    out = np.full(target.shape, np.nan, dtype=float)

    for i, t in enumerate(target):
        if not np.isfinite(t):
            continue

        # Exact/tolerance-level native-center match.
        exact = np.where(np.abs(x - t) <= float(tol))[0]
        if exact.size:
            value = y[int(exact[0])]
            if np.isfinite(value):
                out[i] = float(value)
            continue

        j = int(np.searchsorted(x, t, side="right") - 1)
        if j < 0 or j >= len(x) - 1:
            continue

        x0, x1 = float(x[j]), float(x[j + 1])
        y0, y1 = float(y[j]), float(y[j + 1])

        if not (np.isfinite(y0) and np.isfinite(y1)):
            continue
        if not (x0 - tol <= t <= x1 + tol):
            continue
        if x1 <= x0:
            continue

        w = (float(t) - x0) / (x1 - x0)
        out[i] = y0 + w * (y1 - y0)

    return out


def read_json(path: Path, required: bool = True) -> Dict[str, object]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(path)
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_csv_required(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"{label} is empty: {path}")
    return frame


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{label}: missing required columns {missing}; "
            f"available={list(frame.columns)}"
        )


def bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes", "y"])
    )


def finite_numeric(values) -> np.ndarray:
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    return x[np.isfinite(x)]


def interpolate_inside_support(
    x: np.ndarray,
    y: np.ndarray,
    grid: np.ndarray,
    support_tol: float = SUPPORT_TOL,
) -> np.ndarray:
    """
    Linear interpolation strictly inside native support, allowing only a tiny
    numerical tolerance at support boundaries.

    The Step-8 envelope and configuration-level curves originate from the same
    fixed bin grid, but CSV serialization / median reconstruction can differ at
    ~floating-point precision. Without a tolerance, a boundary that should be
    identical can be classified as infinitesimally outside support.

    Values lying within `support_tol` of the native boundary are clipped to the
    exact boundary before interpolation. This is numerical stabilization, not
    scientific extrapolation.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    grid = np.asarray(grid, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    if len(x) < 2:
        return np.full(grid.shape, np.nan)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    ux, inv = np.unique(x, return_inverse=True)
    if len(ux) != len(x):
        sums = np.zeros(len(ux), dtype=float)
        counts = np.zeros(len(ux), dtype=float)
        for i, j in enumerate(inv):
            sums[j] += y[i]
            counts[j] += 1
        x = ux
        y = sums / counts

    if len(x) < 2:
        return np.full(grid.shape, np.nan)

    xmin = float(np.min(x))
    xmax = float(np.max(x))
    tol = float(support_tol)

    out = np.full(grid.shape, np.nan)

    inside = (
        (grid >= xmin - tol)
        & (grid <= xmax + tol)
    )

    if inside.any():
        # Clip only tolerance-level boundary differences. This avoids numerical
        # pseudo-extrapolation while retaining the intended shared endpoints.
        g = np.clip(
            grid[inside],
            xmin,
            xmax,
        )
        out[inside] = np.interp(
            g,
            x,
            y,
        )

    return out


def summarize_values(values: np.ndarray) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "sd": np.nan,
            "q005": np.nan,
            "q025": np.nan,
            "median": np.nan,
            "q975": np.nan,
            "q995": np.nan,
        }

    q = np.quantile(
        x,
        [0.005, 0.025, 0.5, 0.975, 0.995],
    )

    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "sd": (
            float(np.std(x, ddof=1))
            if x.size > 1
            else 0.0
        ),
        "q005": float(q[0]),
        "q025": float(q[1]),
        "median": float(q[2]),
        "q975": float(q[3]),
        "q995": float(q[4]),
    }


def empirical_p_upper(
    null_values: np.ndarray,
    observed: float,
) -> float:
    x = np.asarray(null_values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0 or not np.isfinite(observed):
        return np.nan

    return float(
        (1 + np.sum(x >= float(observed)))
        / (x.size + 1)
    )


def write_manifest(outdir: Path) -> None:
    rows = []

    for path in sorted(outdir.rglob("*")):
        if path.is_file():
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


# =============================================================================
# Input resolution / validation
# =============================================================================

def resolve_inputs(
    step11_dir: Path,
    pseudo_step8_dir: Path,
) -> Dict[str, object]:
    step11_dir = step11_dir.expanduser().resolve(strict=True)
    step8_dir = pseudo_step8_dir.expanduser().resolve(strict=True)

    paths = {
        "step11_run_config": step11_dir / "00_run_config.json",
        "step11_run_signature": step11_dir / "00_run_signature.json",
        "step11_by_bin_dt": step11_dir / "02_fluctuation_by_bin_dt.csv",
        "step11_global_dt": step11_dir / "03_global_fluctuation_by_dt.csv",
        "step11_subject_sufficient":
            step11_dir / "08_subject_bin_sufficient.parquet",

        "step8_run_config": step8_dir / "00_run_config.json",
        "step8_run_signature": step8_dir / "00_run_signature.json",
        "step8_dt_summary": step8_dir / "02_ensemble_dt_summary.csv",
        "step8_config_binned":
            step8_dir / "03_configuration_binned_metrics_long.csv",
        "step8_binned_envelope":
            step8_dir / "04_binned_randomization_envelope.csv",
    }

    missing = [
        str(path)
        for path in paths.values()
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Required final Step-11/Step-8 files missing:\n  "
            + "\n  ".join(missing)
        )

    step11_cfg = read_json(paths["step11_run_config"])
    step11_sig = read_json(paths["step11_run_signature"])
    step8_cfg = read_json(paths["step8_run_config"])
    step8_sig = read_json(paths["step8_run_signature"])

    # Final Step 11 contract.
    if str(step11_cfg.get("primary_support")) != "common4":
        raise ValueError("Step 11 primary_support must be common4.")
    if str(step11_cfg.get("primary_conditioning")) != "xmid_latent":
        raise ValueError("Step 11 primary_conditioning must be xmid_latent.")
    if step11_cfg.get("operational_class_filter") is not None:
        raise ValueError(
            "Step 12 requires Step 11 without operational-class filtering."
        )
    if int(step11_cfg.get("ddof", DDOF)) != DDOF:
        raise ValueError("Step 11 ddof must equal 1.")

    s11_est = step11_cfg.get("primary_estimands", {})
    if s11_est.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError("Step 11 cross_cov must remain signed/untruncated.")

    s11_payload = step11_sig.get("signature_payload", {})
    if (
        s11_payload.get("primary_support") != "common4"
        or s11_payload.get("primary_conditioning") != "xmid_latent"
        or s11_payload.get("cross_cov_truncated_at_zero") is not False
    ):
        raise ValueError(
            "Step-11 run signature does not certify the final fluctuation contract."
        )

    # Final Step 8 contract.
    if str(step8_cfg.get("primary_support")) != "common4":
        raise ValueError("Step 8 primary_support must be common4.")
    if str(step8_cfg.get("primary_conditioning")) != "xmid_latent":
        raise ValueError("Step 8 primary_conditioning must be xmid_latent.")
    if int(step8_cfg.get("ddof", DDOF)) != DDOF:
        raise ValueError("Step 8 ddof must equal 1.")

    s8_est = step8_cfg.get("primary_estimands", {})
    if s8_est.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError("Step 8 cross_cov must remain signed/untruncated.")

    return {
        "step11_dir": step11_dir,
        "step8_dir": step8_dir,
        "paths": paths,
        "step11_config": step11_cfg,
        "step11_signature": step11_sig,
        "step8_config": step8_cfg,
        "step8_signature": step8_sig,
    }


# =============================================================================
# Step-8 pseudo loading
# =============================================================================

def discover_pseudo_configuration_ids(
    config_binned_path: Path,
) -> List[str]:
    lf = (
        pl.scan_csv(
            config_binned_path,
            infer_schema_length=10000,
        )
        .select(
            pl.col("configuration_id")
            .cast(pl.String, strict=False)
        )
        .unique()
        .sort("configuration_id")
    )

    try:
        frame = lf.collect(engine="streaming")
    except TypeError:
        frame = lf.collect(streaming=True)

    return [
        str(x)
        for x in frame["configuration_id"].to_list()
    ]



def load_full_pseudo_valid_support_cells(
    path: Path,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    """
    Load only the xmid_latent configuration x dt x bin cells that satisfy the
    Step-8 per-bin minimum support.

    This compact table is used to:
      1. count how many pseudo configurations have ANY valid binned support at
         each lag;
      2. define stable pseudo native support relative to that lag-specific
         available configuration set;
      3. choose structurally eligible configurations for smoke testing.

    The denominator is therefore NOT the full 2000 configurations when a lag
    itself has reduced usable support (especially dt=4-5).
    """
    required = [
        "configuration_id",
        "conditioning",
        "dt",
        "bin_id",
        "x_center",
        "meets_min_n",
    ]

    schema = set(
        pl.scan_csv(
            path,
            infer_schema_length=10000,
        )
        .collect_schema()
        .names()
    )
    missing = [c for c in required if c not in schema]
    if missing:
        raise ValueError(
            f"Step-8 configuration binned table missing {missing}"
        )

    lf = (
        pl.scan_csv(
            path,
            infer_schema_length=10000,
        )
        .filter(
            pl.col("conditioning")
            .cast(pl.String, strict=False)
            .eq("xmid_latent")
            & pl.col("dt")
            .cast(pl.Int64, strict=False)
            .is_in([int(x) for x in dt_values])
        )
        .select(
            [
                pl.col("configuration_id")
                .cast(pl.String, strict=False),
                pl.col("dt")
                .cast(pl.Int64, strict=False),
                pl.col("bin_id")
                .cast(pl.Int64, strict=False),
                pl.col("x_center")
                .cast(pl.Float64, strict=False),
                pl.col("meets_min_n")
                .cast(pl.Boolean, strict=False)
                .alias("meets_min_n"),
            ]
        )
        .filter(
            pl.col("meets_min_n").fill_null(False)
            & pl.col("x_center").is_finite()
        )
    )

    try:
        d = lf.collect(engine="streaming").to_pandas()
    except TypeError:
        d = lf.collect(streaming=True).to_pandas()

    if d.empty:
        raise ValueError(
            "No valid xmid_latent Step-8 configuration x bin cells found."
        )

    d["configuration_id"] = d["configuration_id"].astype(str)
    d["dt"] = pd.to_numeric(d["dt"], errors="coerce").astype(int)
    d["bin_id"] = pd.to_numeric(d["bin_id"], errors="coerce").astype(int)
    d["x_center"] = pd.to_numeric(d["x_center"], errors="coerce")
    return d


def pseudo_available_configurations_by_dt(
    valid_cells: pd.DataFrame,
    dt_values: Sequence[int],
) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for dt in dt_values:
        g = valid_cells[
            pd.to_numeric(
                valid_cells["dt"],
                errors="coerce",
            ).eq(int(dt))
        ]
        out[int(dt)] = int(
            g["configuration_id"].astype(str).nunique()
        )
        if out[int(dt)] < 1:
            raise ValueError(
                f"dt={dt}: no pseudo configuration has a valid xmid_latent bin."
            )
    return out


def select_smoke_test_configuration_ids(
    valid_cells: pd.DataFrame,
    support: pd.DataFrame,
    dt_values: Sequence[int],
    max_configs: int,
    min_native_bins: int,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Select a SMALL UNION of pseudo configurations that collectively exercises
    every requested lag on the final matched support.

    IMPORTANT
    ---------
    A smoke test does NOT require the SAME pseudo configuration to span all
    temporal lags. The primary Step-12 inference is lag-specific, and the
    eligible pseudo subset is allowed to differ by dt.

    The previous all-lags requirement was unnecessarily strict and can fail
    even when every lag has many valid configurations.

    Selection algorithm
    -------------------
    1. For each configuration, determine which requested lags it can cover
       completely on the matched support.
    2. Verify that every lag has at least one support-spanning configuration.
    3. Greedily choose configurations that cover the largest number of still
       uncovered lags.
    4. Fill remaining smoke-test slots with configurations covering the largest
       number of requested lags.

    Returns
    -------
    selected_ids
        Union of selected configuration IDs.

    audit
        Per-lag candidate and selected coverage counts.
    """
    max_configs = int(max_configs)
    if max_configs < 1:
        raise ValueError("max_configs must be >=1")

    support_map = {
        int(row.dt): (
            float(row.shared_support_min),
            float(row.shared_support_max),
        )
        for row in support.itertuples(index=False)
    }

    dt_set = [int(dt) for dt in dt_values]

    # Map configuration -> set of lags that it can fully span.
    coverage_by_config: Dict[str, set] = {}

    for cid, cg in valid_cells.groupby(
        "configuration_id",
        sort=True,
    ):
        covered = set()

        for dt in dt_set:
            g = cg[
                pd.to_numeric(
                    cg["dt"],
                    errors="coerce",
                ).eq(int(dt))
            ].copy()

            if len(g) < int(min_native_bins):
                continue

            x = pd.to_numeric(
                g["x_center"],
                errors="coerce",
            ).to_numpy(float)
            x = x[np.isfinite(x)]

            if x.size < int(min_native_bins):
                continue

            lo, hi = support_map[int(dt)]

            # Full interpolation coverage of the shared grid is possible only
            # when the valid native curve spans both matched-support boundaries.
            xmin = float(np.min(x))
            xmax = float(np.max(x))

            if (
                xmin <= float(lo) + SUPPORT_TOL
                and xmax >= float(hi) - SUPPORT_TOL
            ):
                covered.add(int(dt))

        if covered:
            coverage_by_config[str(cid)] = covered

    candidates_by_dt: Dict[int, List[str]] = {}
    for dt in dt_set:
        candidates = sorted(
            cid
            for cid, covered in coverage_by_config.items()
            if int(dt) in covered
        )
        candidates_by_dt[int(dt)] = candidates

        if not candidates:
            # Diagnostic summary: after the lag-relative 80% support rule, a
            # zero-candidate result should not normally occur. Report the native
            # support extremes to distinguish numerical from structural issues.
            dt_cells = valid_cells[
                pd.to_numeric(
                    valid_cells["dt"],
                    errors="coerce",
                ).eq(int(dt))
            ].copy()

            config_ranges = (
                dt_cells.groupby("configuration_id")["x_center"]
                .agg(["min", "max", "count"])
                .reset_index()
            )

            lo, hi = support_map[int(dt)]

            best_low = (
                float(config_ranges["min"].min())
                if not config_ranges.empty
                else np.nan
            )
            best_high = (
                float(config_ranges["max"].max())
                if not config_ranges.empty
                else np.nan
            )

            raise ValueError(
                f"Smoke-test selection found no pseudo configuration spanning "
                f"matched support for dt={dt}: shared=[{lo:.12g}, {hi:.12g}], "
                f"global valid-cell extent=[{best_low:.12g}, {best_high:.12g}], "
                f"support_tol={SUPPORT_TOL:g}. "
                "Inspect Step-8 configuration-level bin support."
            )

    # Greedy set-cover step: cover every lag at least once.
    selected: List[str] = []
    uncovered = set(dt_set)

    while uncovered:
        ranked = sorted(
            (
                (
                    len(coverage_by_config[cid] & uncovered),
                    len(coverage_by_config[cid]),
                    cid,
                )
                for cid in coverage_by_config
                if cid not in selected
            ),
            key=lambda x: (-x[0], -x[1], x[2]),
        )

        if not ranked or ranked[0][0] <= 0:
            raise RuntimeError(
                "Internal smoke-test selection error: remaining lags cannot be covered."
            )

        cid = ranked[0][2]
        selected.append(cid)
        uncovered -= coverage_by_config[cid]

        if len(selected) > max_configs:
            raise ValueError(
                f"--max-pseudo-configs={max_configs} is too small to obtain "
                "lag-wise smoke-test coverage. Increase it."
            )

    # Fill any remaining slots with configurations that cover many lags.
    ranked_all = sorted(
        coverage_by_config,
        key=lambda cid: (
            -len(coverage_by_config[cid]),
            cid,
        ),
    )

    for cid in ranked_all:
        if len(selected) >= max_configs:
            break
        if cid not in selected:
            selected.append(cid)

    audit_rows = []
    for dt in dt_set:
        selected_covering_dt = [
            cid
            for cid in selected
            if int(dt) in coverage_by_config.get(cid, set())
        ]
        audit_rows.append(
            {
                "dt": int(dt),
                "n_full_ensemble_support_spanning_candidates": int(
                    len(candidates_by_dt[int(dt)])
                ),
                "n_selected_configurations_covering_dt": int(
                    len(selected_covering_dt)
                ),
                "selected_configuration_ids_covering_dt": ";".join(
                    selected_covering_dt
                ),
            }
        )

    return selected, pd.DataFrame(audit_rows)


def load_pseudo_configuration_binned(
    path: Path,
    selected_ids: Sequence[str],
    dt_values: Sequence[int],
) -> pd.DataFrame:
    required = [
        "configuration_id",
        "conditioning",
        "dt",
        "bin_id",
        "x_center",
        "n",
        "cross_cov",
        "same_var_mean",
        "replicate_specific_excess",
        "shared_fraction",
        "meets_min_n",
    ]

    schema = set(
        pl.scan_csv(
            path,
            infer_schema_length=10000,
        )
        .collect_schema()
        .names()
    )
    missing = [
        c for c in required
        if c not in schema
    ]
    if missing:
        raise ValueError(
            f"Step-8 configuration binned table missing {missing}"
        )

    lf = (
        pl.scan_csv(
            path,
            infer_schema_length=10000,
        )
        .filter(
            pl.col("configuration_id")
            .cast(pl.String, strict=False)
            .is_in([str(x) for x in selected_ids])
            & pl.col("conditioning")
            .cast(pl.String, strict=False)
            .eq("xmid_latent")
            & pl.col("dt")
            .cast(pl.Int64, strict=False)
            .is_in([int(x) for x in dt_values])
        )
        .select(
            [
                pl.col("configuration_id")
                .cast(pl.String, strict=False),
                pl.col("dt")
                .cast(pl.Int64, strict=False),
                pl.col("bin_id")
                .cast(pl.Int64, strict=False),
                pl.col("x_center")
                .cast(pl.Float64, strict=False),
                pl.col("n")
                .cast(pl.Float64, strict=False),
                pl.col("cross_cov")
                .cast(pl.Float64, strict=False),
                pl.col("same_var_mean")
                .cast(pl.Float64, strict=False),
                pl.col("replicate_specific_excess")
                .cast(pl.Float64, strict=False),
                pl.col("shared_fraction")
                .cast(pl.Float64, strict=False),
                pl.col("meets_min_n")
                .cast(pl.Boolean, strict=False),
            ]
        )
    )

    try:
        d = lf.collect(engine="streaming").to_pandas()
    except TypeError:
        d = lf.collect(streaming=True).to_pandas()

    if d.empty:
        raise ValueError(
            "No xmid_latent Step-8 configuration-level binned rows remain."
        )

    d["configuration_id"] = d["configuration_id"].astype(str)
    return d


# =============================================================================
# Shared support
# =============================================================================

def build_shared_support(
    real: pd.DataFrame,
    pseudo_envelope: pd.DataFrame,
    dt_values: Sequence[int],
    full_pseudo_n: int,
    min_pseudo_bin_valid_fraction: float,
    min_shared_native_bins: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    require_columns(
        real,
        [
            "dt",
            "x_center",
            "cross_cov",
            "meets_min_n",
        ],
        "Step-11 pooled curve",
    )

    require_columns(
        pseudo_envelope,
        [
            "conditioning",
            "dt",
            "bin_id",
            "x_center",
            "n_configurations_valid",
            "cross_cov",
        ],
        "Step-8 binned envelope",
    )

    p = pseudo_envelope[
        pseudo_envelope["conditioning"]
        .astype(str)
        .eq("xmid_latent")
    ].copy()

    for c in [
        "dt",
        "bin_id",
        "x_center",
        "n_configurations_valid",
        "cross_cov",
    ]:
        p[c] = pd.to_numeric(
            p[c],
            errors="coerce",
        )

    support_rows = []
    bin_rows = []

    for dt in dt_values:
        r = real[
            pd.to_numeric(real["dt"], errors="coerce").eq(int(dt))
            & bool_series(real["meets_min_n"])
        ].copy()

        r["x_center"] = pd.to_numeric(
            r["x_center"], errors="coerce"
        )
        r["cross_cov"] = pd.to_numeric(
            r["cross_cov"], errors="coerce"
        )
        r = r[
            np.isfinite(r["x_center"])
            & np.isfinite(r["cross_cov"])
        ]

        q = p[
            p["dt"].eq(int(dt))
            & np.isfinite(p["x_center"])
            & np.isfinite(p["cross_cov"])
        ].copy()

        q["valid_fraction_full_ensemble"] = (
            q["n_configurations_valid"]
            / float(full_pseudo_n)
        )
        q["stable_pseudo_bin"] = (
            q["valid_fraction_full_ensemble"]
            >= float(min_pseudo_bin_valid_fraction)
        )

        if len(r) < 2:
            raise ValueError(
                f"dt={dt}: fewer than 2 valid longitudinal bins."
            )

        real_min = float(r["x_center"].min())
        real_max = float(r["x_center"].max())

        q["inside_longitudinal_support"] = (
            (q["x_center"] >= real_min - SUPPORT_TOL)
            & (q["x_center"] <= real_max + SUPPORT_TOL)
        )
        q["used_in_matched_analysis"] = (
            q["stable_pseudo_bin"].astype(bool)
            & q["inside_longitudinal_support"].astype(bool)
        )

        matched = q[
            q["used_in_matched_analysis"].astype(bool)
        ].sort_values("bin_id")

        if len(matched) < int(min_shared_native_bins):
            raise ValueError(
                f"dt={dt}: only {len(matched)} stable pseudo native bins "
                f"fall inside longitudinal support; "
                f"required={int(min_shared_native_bins)}."
            )

        for row in q.itertuples(index=False):
            bin_rows.append(
                {
                    "dt": int(dt),
                    "bin_id": int(row.bin_id),
                    "x_center": float(row.x_center),
                    "n_configurations_valid": int(
                        row.n_configurations_valid
                    ),
                    "valid_fraction_full_ensemble": float(
                        row.valid_fraction_full_ensemble
                    ),
                    "stable_pseudo_bin": bool(row.stable_pseudo_bin),
                    "inside_longitudinal_support": bool(
                        row.inside_longitudinal_support
                    ),
                    "used_in_matched_analysis": bool(
                        row.used_in_matched_analysis
                    ),
                }
            )

        support_rows.append(
            {
                "dt": int(dt),
                "real_support_min": real_min,
                "real_support_max": real_max,
                "matched_support_min": float(
                    matched["x_center"].min()
                ),
                "matched_support_max": float(
                    matched["x_center"].max()
                ),
                "n_real_valid_bins": int(len(r)),
                "n_pseudo_native_bins": int(len(q)),
                "n_pseudo_stable_bins": int(
                    q["stable_pseudo_bin"].sum()
                ),
                "n_matched_native_bins": int(len(matched)),
                "min_pseudo_bin_valid_fraction": float(
                    min_pseudo_bin_valid_fraction
                ),
                "full_pseudo_configurations": int(full_pseudo_n),
                "matched_bin_ids": ";".join(
                    str(int(x))
                    for x in matched["bin_id"].tolist()
                ),
            }
        )

    return (
        pd.DataFrame(support_rows),
        pd.DataFrame(bin_rows),
    )


# =============================================================================
# Real matched-support interpolation
# =============================================================================

def matched_native_bins_for_dt(
    stable_bins: pd.DataFrame,
    dt: int,
) -> pd.DataFrame:
    d = stable_bins[
        pd.to_numeric(stable_bins["dt"], errors="coerce").eq(int(dt))
        & bool_series(stable_bins["used_in_matched_analysis"])
    ].copy().sort_values("bin_id")

    if d.empty:
        raise ValueError(f"dt={dt}: no matched stable pseudo native bins.")
    return d


def real_curve_on_shared_grid(
    real: pd.DataFrame,
    stable_bins: pd.DataFrame,
    dt: int,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    matched = matched_native_bins_for_dt(
        stable_bins,
        dt=int(dt),
    )

    bin_ids = matched["bin_id"].astype(int).to_numpy()
    target = matched["x_center"].to_numpy(float)

    r = real[
        pd.to_numeric(real["dt"], errors="coerce").eq(int(dt))
        & bool_series(real["meets_min_n"])
    ].copy().sort_values("x_center")

    x_native = pd.to_numeric(
        r["x_center"], errors="coerce"
    ).to_numpy(float)

    out = pd.DataFrame(
        {
            "dt": int(dt),
            "bin_id": bin_ids,
            "x": target,
        }
    )

    for metric in PSEUDO_METRIC_COLUMNS:
        if metric not in r.columns:
            continue

        out[f"longitudinal_{metric}"] = (
            interpolate_adjacent_valid_native(
                x_native,
                pd.to_numeric(
                    r[metric], errors="coerce"
                ).to_numpy(float),
                target,
            )
        )

        lo_col = f"{metric}_ci025"
        hi_col = f"{metric}_ci975"

        if lo_col in r.columns and hi_col in r.columns:
            out[f"longitudinal_{metric}_ci025"] = (
                interpolate_adjacent_valid_native(
                    x_native,
                    pd.to_numeric(
                        r[lo_col], errors="coerce"
                    ).to_numpy(float),
                    target,
                )
            )
            out[f"longitudinal_{metric}_ci975"] = (
                interpolate_adjacent_valid_native(
                    x_native,
                    pd.to_numeric(
                        r[hi_col], errors="coerce"
                    ).to_numpy(float),
                    target,
                )
            )

    if not np.isfinite(
        out["longitudinal_cross_cov"].to_numpy(float)
    ).all():
        raise ValueError(
            f"dt={dt}: longitudinal Step-11 curve cannot be locally "
            "interpolated onto every matched Step-8 native center without "
            "bridging unsupported bins."
        )

    return bin_ids, target, out


# =============================================================================
# Pseudo matched-support null
# =============================================================================

def pseudo_null_for_dt(
    pseudo: pd.DataFrame,
    selected_ids: Sequence[str],
    dt: int,
    matched_bin_ids: np.ndarray,
    matched_x: np.ndarray,
    min_bin_coverage: float,
    execution_min_configurations: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = pseudo[
        pd.to_numeric(pseudo["dt"], errors="coerce").eq(int(dt))
        & bool_series(pseudo["meets_min_n"])
    ].copy()

    grouped = {
        cid: g.copy()
        for cid, g in d.groupby(
            "configuration_id",
            sort=False,
        )
    }

    bin_ids = [int(x) for x in matched_bin_ids]
    bin_to_pos = {
        int(b): i
        for i, b in enumerate(bin_ids)
    }

    metric_rows = []
    cross_matrix = np.full(
        (len(selected_ids), len(bin_ids)),
        np.nan,
        dtype=float,
    )

    for i, cid in enumerate(selected_ids):
        g = grouped.get(str(cid))

        if g is None:
            metric_rows.append(
                {
                    "configuration_id": str(cid),
                    "dt": int(dt),
                    "eligible_for_null": False,
                    "reason": "no_valid_cells",
                    "n_matched_bins_valid": 0,
                    "matched_bin_coverage_fraction": 0.0,
                }
            )
            continue

        g = g[
            pd.to_numeric(
                g["bin_id"], errors="coerce"
            ).isin(bin_ids)
        ].copy()

        if g["bin_id"].duplicated().any():
            raise ValueError(
                f"{cid} dt={dt}: duplicate Step-8 native bin IDs."
            )

        values_by_metric = {
            metric: np.full(
                len(bin_ids),
                np.nan,
                dtype=float,
            )
            for metric in PSEUDO_METRIC_COLUMNS
        }

        for row in g.itertuples(index=False):
            b = int(row.bin_id)
            if b not in bin_to_pos:
                continue
            j = bin_to_pos[b]
            for metric in PSEUDO_METRIC_COLUMNS:
                value = float(getattr(row, metric))
                if np.isfinite(value):
                    values_by_metric[metric][j] = value

        finite = np.isfinite(
            values_by_metric["cross_cov"]
        )
        n_valid = int(finite.sum())
        coverage = float(
            n_valid / len(bin_ids)
        )

        eligible = (
            coverage >= float(min_bin_coverage)
        )

        row_out = {
            "configuration_id": str(cid),
            "dt": int(dt),
            "eligible_for_null": bool(eligible),
            "reason": (
                "eligible"
                if eligible
                else "incomplete_direct_native_bin_coverage"
            ),
            "n_matched_bins_required": int(len(bin_ids)),
            "n_matched_bins_valid": n_valid,
            "matched_bin_coverage_fraction": coverage,
            "complete_matched_native_bins": bool(
                n_valid == len(bin_ids)
            ),
        }

        for metric in PSEUDO_METRIC_COLUMNS:
            values = values_by_metric[metric]
            if eligible and np.isfinite(values).all():
                row_out[
                    f"abundance_averaged_{metric}"
                ] = float(np.mean(values))
            elif eligible:
                # cross_cov coverage defines eligibility; secondary metrics can
                # still be missing and are reported as NA rather than imputed.
                finite_metric = values[np.isfinite(values)]
                row_out[
                    f"abundance_averaged_{metric}"
                ] = (
                    float(np.mean(finite_metric))
                    if finite_metric.size == len(values)
                    else np.nan
                )
            else:
                row_out[
                    f"abundance_averaged_{metric}"
                ] = np.nan

        metric_rows.append(row_out)

        if eligible:
            cross_matrix[i, :] = values_by_metric["cross_cov"]

    metrics = pd.DataFrame(metric_rows)

    eligible_n = int(
        metrics["eligible_for_null"]
        .astype(bool)
        .sum()
    )

    if eligible_n < int(execution_min_configurations):
        reason_counts = (
            metrics["reason"]
            .value_counts(dropna=False)
            .to_dict()
        )
        raise ValueError(
            f"dt={dt}: eligible pseudo configurations={eligible_n}, "
            f"required_for_execution={execution_min_configurations}; "
            f"reason_counts={reason_counts}."
        )

    pointwise = pd.DataFrame(
        {
            "dt": int(dt),
            "bin_id": bin_ids,
            "x": np.asarray(matched_x, dtype=float),
            "pseudo_cross_cov_median": np.nanmedian(
                cross_matrix,
                axis=0,
            ),
            "pseudo_cross_cov_q025": np.nanquantile(
                cross_matrix,
                0.025,
                axis=0,
            ),
            "pseudo_cross_cov_q975": np.nanquantile(
                cross_matrix,
                0.975,
                axis=0,
            ),
            "pseudo_n_configurations": np.sum(
                np.isfinite(cross_matrix),
                axis=0,
            ),
        }
    )

    return metrics, pointwise


# =============================================================================
# Step-11 subject bootstrap on matched support
# =============================================================================

def metrics_from_sums(
    n: float,
    sum_A: float,
    sum_B: float,
    sum_A2: float,
    sum_B2: float,
    sum_AB: float,
) -> Dict[str, float]:
    n = float(n)

    if not np.isfinite(n) or n <= DDOF:
        return {
            "cross_cov": np.nan,
            "same_var_mean": np.nan,
            "replicate_specific_excess": np.nan,
            "shared_fraction": np.nan,
        }

    denom = n - DDOF

    var_A = (
        float(sum_A2)
        - float(sum_A) ** 2 / n
    ) / denom

    var_B = (
        float(sum_B2)
        - float(sum_B) ** 2 / n
    ) / denom

    cross_cov = (
        float(sum_AB)
        - float(sum_A) * float(sum_B) / n
    ) / denom

    if var_A < 0 and abs(var_A) < 1e-12:
        var_A = 0.0
    if var_B < 0 and abs(var_B) < 1e-12:
        var_B = 0.0

    same = 0.5 * (
        var_A + var_B
    )
    excess = same - cross_cov

    shared = (
        cross_cov / same
        if np.isfinite(same)
        and same > 0
        else np.nan
    )

    return {
        "cross_cov": float(cross_cov),
        "same_var_mean": float(same),
        "replicate_specific_excess": float(excess),
        "shared_fraction": float(shared),
    }


def subject_sufficient_arrays(
    path: Path,
    subjects: Sequence[int],
    dt_values: Sequence[int],
    n_bins: int,
) -> Dict[str, np.ndarray]:
    d = pl.read_parquet(path)

    required = [
        "subject",
        "dt",
        "bin",
        "n",
        "sum_A",
        "sum_B",
        "sum_A2",
        "sum_B2",
        "sum_AB",
    ]
    missing = [
        c for c in required
        if c not in d.columns
    ]
    if missing:
        raise ValueError(
            f"Step-11 subject sufficient table missing {missing}"
        )

    s_index = {
        int(s): i
        for i, s in enumerate(subjects)
    }
    d_index = {
        int(dt): i
        for i, dt in enumerate(dt_values)
    }

    shape = (
        len(subjects),
        len(dt_values),
        int(n_bins),
    )

    cols = [
        "n",
        "sum_A",
        "sum_B",
        "sum_A2",
        "sum_B2",
        "sum_AB",
    ]

    arrays = {
        c: np.zeros(
            shape,
            dtype=float,
        )
        for c in cols
    }

    for row in d.iter_rows(named=True):
        subject = int(row["subject"])
        dt = int(row["dt"])
        b = int(row["bin"])

        if (
            subject not in s_index
            or dt not in d_index
            or b < 0
            or b >= n_bins
        ):
            continue

        si = s_index[subject]
        di = d_index[dt]

        for c in cols:
            arrays[c][si, di, b] = float(
                row[c]
            )

    return arrays


def longitudinal_bootstrap_matched_support(
    sufficient_path: Path,
    real_curve: pd.DataFrame,
    stable_bins: pd.DataFrame,
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
    min_n: int,
) -> pd.DataFrame:
    subjects = sorted(
        int(x)
        for x in pl.read_parquet(
            sufficient_path,
            columns=["subject"],
        )["subject"]
        .unique()
        .to_list()
    )

    bins = pd.to_numeric(
        real_curve["bin"],
        errors="coerce",
    ).dropna().astype(int)

    n_bins = int(
        bins.max() + 1
    )

    arrays = subject_sufficient_arrays(
        sufficient_path,
        subjects,
        dt_values,
        n_bins,
    )

    centers_by_dt = {}
    targets_by_dt = {}

    for dt in dt_values:
        d = real_curve[
            pd.to_numeric(
                real_curve["dt"],
                errors="coerce",
            ).eq(int(dt))
        ].copy().sort_values("bin")

        centers = np.full(
            n_bins,
            np.nan,
            dtype=float,
        )
        for row in d.itertuples(index=False):
            centers[int(row.bin)] = float(
                row.x_center
            )
        centers_by_dt[int(dt)] = centers

        targets_by_dt[int(dt)] = (
            matched_native_bins_for_dt(
                stable_bins,
                int(dt),
            )["x_center"].to_numpy(float)
        )

    rng = np.random.default_rng(
        int(seed)
    )
    rows = []

    for boot in range(int(n_bootstrap)):
        draw = rng.integers(
            0,
            len(subjects),
            size=len(subjects),
        )
        mult = np.bincount(
            draw,
            minlength=len(subjects),
        ).astype(float)

        lag_values = []

        for di, dt in enumerate(dt_values):
            sums = {
                c: np.tensordot(
                    mult,
                    arr[:, di, :],
                    axes=(0, 0),
                )
                for c, arr in arrays.items()
            }

            cross = np.full(
                n_bins,
                np.nan,
                dtype=float,
            )

            for b in range(n_bins):
                if (
                    not np.isfinite(sums["n"][b])
                    or sums["n"][b] < int(min_n)
                ):
                    continue

                m = metrics_from_sums(
                    sums["n"][b],
                    sums["sum_A"][b],
                    sums["sum_B"][b],
                    sums["sum_A2"][b],
                    sums["sum_B2"][b],
                    sums["sum_AB"][b],
                )
                cross[b] = m["cross_cov"]

            curve = interpolate_adjacent_valid_native(
                centers_by_dt[int(dt)],
                cross,
                targets_by_dt[int(dt)],
            )

            complete = bool(
                np.isfinite(curve).all()
            )
            value = (
                float(np.mean(curve))
                if complete
                else np.nan
            )

            rows.append(
                {
                    "bootstrap": int(boot),
                    "dt": int(dt),
                    "abundance_averaged_cross_cov": value,
                    "complete_matched_support": complete,
                }
            )
            lag_values.append(value)

        lag_values = np.asarray(
            lag_values,
            dtype=float,
        )

        rows.append(
            {
                "bootstrap": int(boot),
                "dt": 0,
                "abundance_averaged_cross_cov": (
                    float(np.mean(lag_values))
                    if np.isfinite(lag_values).all()
                    else np.nan
                ),
                "complete_matched_support": bool(
                    np.isfinite(lag_values).all()
                ),
            }
        )

    return pd.DataFrame(rows)


# =============================================================================
# Main comparison
# =============================================================================

def run_comparison(
    real: pd.DataFrame,
    pseudo: pd.DataFrame,
    stable_bins: pd.DataFrame,
    selected_ids: Sequence[str],
    dt_values: Sequence[int],
    args: argparse.Namespace,
    execution_min_configurations: int,
    inferential_min_configurations: int,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    all_pointwise = []
    all_pseudo_metrics = []
    test_rows = []
    audit_rows = []

    for dt in dt_values:
        matched_bin_ids, matched_x, real_grid = (
            real_curve_on_shared_grid(
                real,
                stable_bins,
                dt=int(dt),
            )
        )

        pseudo_metrics, pseudo_pointwise = (
            pseudo_null_for_dt(
                pseudo,
                selected_ids,
                dt=int(dt),
                matched_bin_ids=matched_bin_ids,
                matched_x=matched_x,
                min_bin_coverage=float(
                    args.min_null_bin_coverage
                ),
                execution_min_configurations=int(
                    execution_min_configurations
                ),
            )
        )

        merged = real_grid.merge(
            pseudo_pointwise,
            on=["dt", "bin_id", "x"],
            how="inner",
            validate="one_to_one",
        )

        merged[
            "longitudinal_minus_pseudo_cross_cov"
        ] = (
            merged["longitudinal_cross_cov"]
            - merged["pseudo_cross_cov_median"]
        )

        merged[
            "longitudinal_outside_pseudo_95"
        ] = (
            merged["longitudinal_cross_cov"]
            > merged["pseudo_cross_cov_q975"]
        ) | (
            merged["longitudinal_cross_cov"]
            < merged["pseudo_cross_cov_q025"]
        )

        real_values = merged[
            "longitudinal_cross_cov"
        ].to_numpy(float)

        if not np.isfinite(real_values).all():
            raise RuntimeError(
                f"dt={dt}: non-finite longitudinal matched-bin covariance."
            )

        real_value = float(
            np.mean(real_values)
        )

        eligible = pseudo_metrics[
            pseudo_metrics[
                "eligible_for_null"
            ].astype(bool)
        ].copy()

        null = pd.to_numeric(
            eligible[
                "abundance_averaged_cross_cov"
            ],
            errors="coerce",
        ).to_numpy(float)
        null = null[
            np.isfinite(null)
        ]

        summary = summarize_values(
            null
        )

        if int(summary["n"]) < int(
            inferential_min_configurations
        ):
            p_emp = np.nan
            status = (
                "insufficient_null_configurations"
            )
        else:
            p_emp = empirical_p_upper(
                null,
                real_value,
            )
            status = "complete"

        test_rows.append(
            {
                "dt": int(dt),
                "metric": (
                    "abundance_averaged_cross_cov"
                ),
                "longitudinal_value": real_value,
                "n_matched_native_bins": int(len(merged)),
                "pseudo_n": int(summary["n"]),
                "pseudo_mean": summary["mean"],
                "pseudo_sd": summary["sd"],
                "pseudo_q005": summary["q005"],
                "pseudo_q025": summary["q025"],
                "pseudo_median": summary["median"],
                "pseudo_q975": summary["q975"],
                "pseudo_q995": summary["q995"],
                "tail": "upper",
                "empirical_p": p_emp,
                "status": status,
            }
        )

        audit_rows.append(
            {
                "dt": int(dt),
                "matched_support_min": float(
                    np.min(matched_x)
                ),
                "matched_support_max": float(
                    np.max(matched_x)
                ),
                "n_matched_native_bins": int(
                    len(matched_x)
                ),
                "n_selected_pseudo_configurations": int(
                    len(selected_ids)
                ),
                "n_eligible_pseudo_configurations": int(
                    len(eligible)
                ),
                "fraction_eligible": float(
                    len(eligible)
                    / len(selected_ids)
                ),
                "pseudo_pointwise_n_min": int(
                    merged["pseudo_n_configurations"].min()
                ),
                "pseudo_pointwise_n_max": int(
                    merged["pseudo_n_configurations"].max()
                ),
                "n_longitudinal_points_above_pseudo_q975": int(
                    (
                        merged[
                            "longitudinal_cross_cov"
                        ]
                        > merged[
                            "pseudo_cross_cov_q975"
                        ]
                    ).sum()
                ),
                "fraction_longitudinal_points_above_pseudo_q975": float(
                    (
                        merged[
                            "longitudinal_cross_cov"
                        ]
                        > merged[
                            "pseudo_cross_cov_q975"
                        ]
                    ).mean()
                ),
            }
        )

        all_pointwise.append(
            merged
        )
        all_pseudo_metrics.append(
            pseudo_metrics
        )

    return (
        pd.concat(
            all_pointwise,
            ignore_index=True,
        ),
        pd.concat(
            all_pseudo_metrics,
            ignore_index=True,
        ),
        pd.DataFrame(test_rows),
        pd.DataFrame(audit_rows),
    )


# =============================================================================
# Joint across-lag comparison
# =============================================================================

def joint_comparison(
    tests: pd.DataFrame,
    pseudo_metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    dt_values: Sequence[int],
    inferential_min_configurations: int,
) -> pd.DataFrame:
    real_values = (
        tests.set_index("dt")
        .loc[
            [int(x) for x in dt_values],
            "longitudinal_value",
        ]
        .to_numpy(float)
    )

    real_joint = float(
        np.mean(real_values)
    )

    eligible = pseudo_metrics[
        pseudo_metrics[
            "eligible_for_null"
        ].astype(bool)
    ].copy()

    pivot = eligible.pivot_table(
        index="configuration_id",
        columns="dt",
        values="abundance_averaged_cross_cov",
        aggfunc="first",
    )

    needed = [
        int(x) for x in dt_values
    ]

    if all(
        dt in pivot.columns
        for dt in needed
    ):
        joint_null = (
            pivot[needed]
            .dropna()
            .mean(axis=1)
            .to_numpy(float)
        )
    else:
        joint_null = np.asarray(
            [],
            dtype=float,
        )

    null_summary = summarize_values(
        joint_null
    )

    if int(
        null_summary["n"]
    ) >= int(
        inferential_min_configurations
    ):
        p_emp = empirical_p_upper(
            joint_null,
            real_joint,
        )
        status = "complete"
    else:
        p_emp = np.nan
        status = (
            "insufficient_complete_configurations"
        )

    real_boot = bootstrap[
        bootstrap["dt"].eq(0)
    ]["abundance_averaged_cross_cov"]
    real_boot = finite_numeric(
        real_boot
    )

    return pd.DataFrame(
        [
            {
                "metric": (
                    "equal_lag_mean_abundance_averaged_cross_cov"
                ),
                "longitudinal_value": real_joint,
                "longitudinal_bootstrap_n": int(
                    real_boot.size
                ),
                "longitudinal_bootstrap_q025": (
                    float(
                        np.quantile(
                            real_boot,
                            0.025,
                        )
                    )
                    if real_boot.size
                    else np.nan
                ),
                "longitudinal_bootstrap_q975": (
                    float(
                        np.quantile(
                            real_boot,
                            0.975,
                        )
                    )
                    if real_boot.size
                    else np.nan
                ),
                "pseudo_n_complete_configurations": int(
                    null_summary["n"]
                ),
                "pseudo_q025": null_summary["q025"],
                "pseudo_median": null_summary["median"],
                "pseudo_q975": null_summary["q975"],
                "empirical_p": p_emp,
                "tail": "upper",
                "status": status,
            }
        ]
    )


# =============================================================================
# Native global descriptive output
# =============================================================================

def native_global_comparison(
    step11_global: pd.DataFrame,
    step8_dt_summary: pd.DataFrame,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    require_columns(
        step11_global,
        ["dt", "cross_cov", "same_var_mean"],
        "Step-11 global by dt",
    )

    require_columns(
        step8_dt_summary,
        [
            "dt",
            "metric",
            "median",
            "randomization_q025",
            "randomization_q975",
            "n_configurations",
        ],
        "Step-8 dt summary",
    )

    real = step11_global[
        pd.to_numeric(
            step11_global["dt"],
            errors="coerce",
        ).isin(
            [int(x) for x in dt_values]
        )
    ].copy()

    pseudo_cross = step8_dt_summary[
        step8_dt_summary["metric"]
        .astype(str)
        .eq("cross_cov")
        & pd.to_numeric(
            step8_dt_summary["dt"],
            errors="coerce",
        ).isin(
            [int(x) for x in dt_values]
        )
    ][
        [
            "dt",
            "n_configurations",
            "median",
            "randomization_q025",
            "randomization_q975",
        ]
    ].copy().rename(
        columns={
            "n_configurations":
                "pseudo_n_configurations_cross_cov",
            "median":
                "pseudo_cross_cov_median",
            "randomization_q025":
                "pseudo_cross_cov_q025",
            "randomization_q975":
                "pseudo_cross_cov_q975",
        }
    )

    return real.merge(
        pseudo_cross,
        on="dt",
        how="left",
    ).sort_values("dt")


# =============================================================================
# README
# =============================================================================

def write_readme(outdir: Path) -> None:
    text = """# Step 12 — longitudinal versus pseudo cross-replicate fluctuations

Primary estimand
----------------
Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt, common4)

The same signed cross-replicate covariance estimand is compared between the
genuine longitudinal cohort and the pseudo-longitudinal technical reference.

No TT/TF/FT/FF restriction is introduced.

Primary comparison
------------------
The real and pseudo curves are compared on matched ABSOLUTE xmid_latent support
separately for dt=1-5.

Stable pseudo support is defined relative to the COMPLETE Step-8 ensemble at
each lag. Configurations with no valid bin remain in the denominator, so the
stability criterion is not conditionally relaxed at long pseudo-lags.

No percentile normalization is performed because covariance magnitude is
abundance-dependent and percentile matching would compare different absolute
abundance levels between datasets.

Missing pseudo native bins are never interpolated or imputed. Production
pseudo null curves require direct validity at every matched stable pseudo bin.

Because Step 11 and Step 8 use different xmid_latent grids, only the genuine
longitudinal curve is interpolated onto matched Step-8 native centers.
Interpolation is local between adjacent valid Step-11 native bins and never
bridges an unsupported internal bin.

Pseudo inference
----------------
Each Step-8 configuration is a randomization unit conditional on one biological
source sample.

The primary lag-specific empirical statistic is the equal-bin average of
cross_cov over the matched stable Step-8 native xmid centers.

Empirical P values use the finite-sample +1 upper-tail test.

Pseudo randomization intervals are not biological confidence intervals.

Longitudinal uncertainty
------------------------
Biological uncertainty is quantified by subject-cluster bootstrap using the
Step-11 subject x dt x bin sufficient statistics.

The same subject draw is shared across temporal lags.

Joint summary
-------------
An equal-lag mean matched-support cross covariance is also reported using only
pseudo configurations eligible at every requested lag.

This summarizes overall real-versus-technical-null fluctuation magnitude and
must not be interpreted as a temporal-scaling test. The joint summary is
reported only when pseudo configurations are eligible at every requested lag;
the primary Step-12 tests remain the lag-specific comparisons. Temporal
accumulation is analyzed in Step 13.
"""
    (outdir / "README_outputs.md").write_text(
        text,
        encoding="utf-8",
    )


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 12: compare genuine longitudinal cross-replicate fluctuation "
            "covariance with the pseudo technical null on matched absolute "
            "xmid_latent support."
        ),
    )

    p.add_argument(
        "--step11-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--pseudo-step8-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    p.add_argument(
        "--dt-values",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5],
    )

    p.add_argument(
        "--min-pseudo-bin-valid-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum fraction of pseudo configurations having any valid "
            "binned support at that lag that must also have a valid estimate "
            "in a given bin for that bin to define stable pseudo native support."
        ),
    )

    p.add_argument(
        "--min-shared-native-bins",
        type=int,
        default=8,
        help="Minimum stable Step-8 native bins retained per lag.",
    )

    p.add_argument(
        "--min-null-bin-coverage",
        type=float,
        default=1.0,
        help=(
            "Required fraction of matched stable Step-8 native bins directly "
            "valid in one pseudo configuration. Production requires 1.0."
        ),
    )

    p.add_argument(
        "--min-null-configurations",
        type=int,
        default=30,
    )

    p.add_argument(
        "--max-pseudo-configs",
        type=int,
        default=None,
        help="Smoke testing only; omit for full production ensemble.",
    )

    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=None,
        help="Longitudinal subject bootstrap count; default inherits Step 11.",
    )

    p.add_argument(
        "--bootstrap-seed",
        type=int,
        default=None,
        help="Longitudinal subject bootstrap seed; default inherits Step 11.",
    )

    return p


def validate_args(args: argparse.Namespace) -> None:
    dt_values = sorted(
        set(int(x) for x in args.dt_values)
    )

    if any(dt < 1 for dt in dt_values):
        raise ValueError("--dt-values must all be >=1")
    if not (
        0 < float(
            args.min_pseudo_bin_valid_fraction
        ) <= 1
    ):
        raise ValueError(
            "--min-pseudo-bin-valid-fraction must be in (0,1]"
        )
    if int(args.min_shared_native_bins) < 3:
        raise ValueError("--min-shared-native-bins must be >=3")
    if not (
        0 < float(
            args.min_null_bin_coverage
        ) <= 1
    ):
        raise ValueError(
            "--min-null-bin-coverage must be in (0,1]"
        )
    if int(args.min_null_configurations) < 1:
        raise ValueError("--min-null-configurations must be >=1")
    if (
        args.max_pseudo_configs is not None
        and int(args.max_pseudo_configs) < 1
    ):
        raise ValueError("--max-pseudo-configs must be >=1")

    if (
        args.max_pseudo_configs is None
        and not np.isclose(
            float(args.min_null_bin_coverage),
            1.0,
            rtol=0.0,
            atol=1e-12,
        )
    ):
        raise ValueError(
            "Production Step-12 inference requires "
            "--min-null-bin-coverage 1.0 so pseudo configurations are never "
            "completed by interpolation or partial matched support."
        )


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)

    inputs = resolve_inputs(
        args.step11_dir,
        args.pseudo_step8_dir,
    )

    outdir = ensure_dir(
        args.outdir.expanduser().resolve()
    )

    dt_values = sorted(
        set(int(x) for x in args.dt_values)
    )

    step11_cfg = inputs["step11_config"]
    step8_cfg = inputs["step8_config"]

    # Bootstrap settings inherited from Step 11 unless overridden.
    step11_boot = step11_cfg.get(
        "bootstrap",
        {},
    )

    n_bootstrap = (
        int(args.n_bootstrap)
        if args.n_bootstrap is not None
        else int(
            step11_boot.get(
                "n_bootstrap",
                2000,
            )
        )
    )

    bootstrap_seed = (
        int(args.bootstrap_seed)
        if args.bootstrap_seed is not None
        else int(
            step11_boot.get(
                "seed",
                123,
            )
        )
    )

    # Full pseudo ensemble size.
    full_pseudo_n = int(
        step8_cfg.get(
            "n_configurations_selected",
            step8_cfg.get(
                "expected_configs_guard",
                2000,
            ),
        )
    )

    all_ids = discover_pseudo_configuration_ids(
        inputs["paths"]["step8_config_binned"]
    )

    if full_pseudo_n > 0 and len(all_ids) != full_pseudo_n:
        raise ValueError(
            f"Step-8 run config reports {full_pseudo_n} configurations "
            f"but configuration-level binned output contains {len(all_ids)}."
        )

    test_mode = args.max_pseudo_configs is not None

    inferential_min_configurations = int(
        args.min_null_configurations
    )
    execution_min_configurations = (
        1
        if test_mode
        else inferential_min_configurations
    )

    print("=== ClonoDynamics Step 12 ===")
    print("script version :", SCRIPT_VERSION)
    print("Step 11        :", inputs["step11_dir"])
    print("pseudo Step 8  :", inputs["step8_dir"])
    print("dt values      :", dt_values)
    print("full pseudo N  :", full_pseudo_n)
    print(
        "selected pseudo:",
        "pending support-aware smoke selection"
        if test_mode
        else len(all_ids),
    )
    print("primary support: common4")
    print("conditioning   : xmid_latent")
    print("comparison     : matched stable Step-8 native xmid bins")
    print("pseudo interpolation: NO")
    print("longitudinal interpolation: adjacent valid Step-11 bins only")
    print("percentile normalization: NO")

    if test_mode:
        print(
            "[TEST MODE] structural execution requires >=1 eligible "
            "configuration per lag."
        )
        print(
            "[TEST MODE] empirical P values still require",
            inferential_min_configurations,
            "eligible configurations.",
        )

    # -----------------------------------------------------------------
    # Load Step 11 / Step 8 summaries.
    # -----------------------------------------------------------------
    real = read_csv_required(
        inputs["paths"]["step11_by_bin_dt"],
        "Step-11 fluctuation by bin/dt",
    )
    real_global = read_csv_required(
        inputs["paths"]["step11_global_dt"],
        "Step-11 global fluctuation by dt",
    )
    pseudo_envelope = read_csv_required(
        inputs["paths"]["step8_binned_envelope"],
        "Step-8 binned randomization envelope",
    )
    pseudo_dt_summary = read_csv_required(
        inputs["paths"]["step8_dt_summary"],
        "Step-8 dt summary",
    )

    # -----------------------------------------------------------------
    # Stable full-ensemble support.
    #
    # The pseudo denominator is lag-specific: configurations with at least one
    # valid xmid_latent bin at that lag. This is essential for long pseudo-lags,
    # where not all 2000 randomizations have adequate common4 support.
    # -----------------------------------------------------------------
    full_valid_support_cells = load_full_pseudo_valid_support_cells(
        inputs["paths"]["step8_config_binned"],
        dt_values,
    )

    support, stable_bins = build_shared_support(
        real,
        pseudo_envelope,
        dt_values,
        full_pseudo_n=full_pseudo_n,
        min_pseudo_bin_valid_fraction=float(
            args.min_pseudo_bin_valid_fraction
        ),
        min_shared_native_bins=int(
            args.min_shared_native_bins
        ),
    )

    # Smoke-test selection is support-aware on the FINAL matched native bins.
    # Production always uses the complete Step-8 ensemble.
    if test_mode:
        coverage_by_config: Dict[str, set] = {}

        stable_map = {
            int(dt): set(
                matched_native_bins_for_dt(
                    stable_bins,
                    int(dt),
                )["bin_id"].astype(int).tolist()
            )
            for dt in dt_values
        }

        for cid, cg in full_valid_support_cells.groupby(
            "configuration_id",
            sort=True,
        ):
            covered = set()
            for dt in dt_values:
                present = set(
                    pd.to_numeric(
                        cg.loc[
                            pd.to_numeric(
                                cg["dt"], errors="coerce"
                            ).eq(int(dt)),
                            "bin_id",
                        ],
                        errors="coerce",
                    )
                    .dropna()
                    .astype(int)
                    .tolist()
                )
                if stable_map[int(dt)].issubset(present):
                    covered.add(int(dt))
            if covered:
                coverage_by_config[str(cid)] = covered

        selected_ids = []
        uncovered = set(int(dt) for dt in dt_values)

        while uncovered:
            ranked = sorted(
                (
                    (
                        len(covered & uncovered),
                        len(covered),
                        cid,
                    )
                    for cid, covered in coverage_by_config.items()
                    if cid not in selected_ids
                ),
                key=lambda x: (-x[0], -x[1], x[2]),
            )
            if not ranked or ranked[0][0] <= 0:
                raise ValueError(
                    "Smoke-test selection cannot directly cover all matched "
                    "stable native bins at every requested lag."
                )
            cid = ranked[0][2]
            selected_ids.append(cid)
            uncovered -= coverage_by_config[cid]

            if len(selected_ids) > int(args.max_pseudo_configs):
                raise ValueError(
                    f"--max-pseudo-configs={args.max_pseudo_configs} is too "
                    "small for direct matched-bin smoke coverage."
                )

        ranked_all = sorted(
            coverage_by_config,
            key=lambda cid: (
                -len(coverage_by_config[cid]),
                cid,
            ),
        )
        for cid in ranked_all:
            if len(selected_ids) >= int(args.max_pseudo_configs):
                break
            if cid not in selected_ids:
                selected_ids.append(cid)

        smoke_rows = []
        for dt in dt_values:
            covering = [
                cid
                for cid in selected_ids
                if int(dt) in coverage_by_config.get(cid, set())
            ]
            smoke_rows.append(
                {
                    "dt": int(dt),
                    "n_selected_configurations_covering_dt": int(
                        len(covering)
                    ),
                    "selected_configuration_ids_covering_dt":
                        ";".join(covering),
                }
            )
        smoke_selection_audit = pd.DataFrame(smoke_rows)

        print(
            "[TEST MODE] direct matched-bin pseudo configurations selected:",
            len(selected_ids),
        )
    else:
        selected_ids = all_ids
        smoke_selection_audit = pd.DataFrame()

    pseudo = load_pseudo_configuration_binned(
        inputs["paths"]["step8_config_binned"],
        selected_ids,
        dt_values,
    )

    # -----------------------------------------------------------------
    # Matched-support real-versus-pseudo comparison.
    # -----------------------------------------------------------------
    pointwise, pseudo_metrics, tests, audit = (
        run_comparison(
            real,
            pseudo,
            stable_bins,
            selected_ids,
            dt_values,
            args,
            execution_min_configurations,
            inferential_min_configurations,
        )
    )

    # -----------------------------------------------------------------
    # Longitudinal subject-cluster bootstrap of matched-support scalar.
    # -----------------------------------------------------------------
    bootstrap = longitudinal_bootstrap_matched_support(
        inputs["paths"]["step11_subject_sufficient"],
        real,
        stable_bins,
        dt_values,
        n_bootstrap=int(n_bootstrap),
        seed=int(bootstrap_seed),
        min_n=int(
            step11_cfg.get("binning", {}).get(
                "min_n_per_pooled_cell",
                50,
            )
        ),
    )

    bootstrap_summary_rows = []

    for dt in dt_values:
        values = finite_numeric(
            bootstrap[
                bootstrap["dt"].eq(int(dt))
            ][
                "abundance_averaged_cross_cov"
            ]
        )
        bootstrap_summary_rows.append(
            {
                "dt": int(dt),
                "n_bootstrap_valid": int(
                    values.size
                ),
                "bootstrap_median": (
                    float(np.median(values))
                    if values.size else np.nan
                ),
                "bootstrap_q025": (
                    float(
                        np.quantile(
                            values,
                            0.025,
                        )
                    )
                    if values.size else np.nan
                ),
                "bootstrap_q975": (
                    float(
                        np.quantile(
                            values,
                            0.975,
                        )
                    )
                    if values.size else np.nan
                ),
            }
        )

    bootstrap_summary = pd.DataFrame(
        bootstrap_summary_rows
    )

    tests = tests.merge(
        bootstrap_summary,
        on="dt",
        how="left",
        validate="one_to_one",
    )

    # -----------------------------------------------------------------
    # Joint across-lag summary.
    # -----------------------------------------------------------------
    joint = joint_comparison(
        tests,
        pseudo_metrics,
        bootstrap,
        dt_values,
        inferential_min_configurations,
    )

    # -----------------------------------------------------------------
    # Native global descriptive comparison.
    # -----------------------------------------------------------------
    native_global = native_global_comparison(
        real_global,
        pseudo_dt_summary,
        dt_values,
    )

    # -----------------------------------------------------------------
    # Analysis summary.
    # -----------------------------------------------------------------
    summary_rows = []

    for row in tests.itertuples(index=False):
        summary_rows.append(
            {
                "scope": f"dt{int(row.dt)}",
                "longitudinal_abundance_averaged_cross_cov":
                    float(row.longitudinal_value),
                "longitudinal_bootstrap_q025":
                    row.bootstrap_q025,
                "longitudinal_bootstrap_q975":
                    row.bootstrap_q975,
                "pseudo_n":
                    int(row.pseudo_n),
                "pseudo_median":
                    row.pseudo_median,
                "pseudo_q025":
                    row.pseudo_q025,
                "pseudo_q975":
                    row.pseudo_q975,
                "empirical_p":
                    row.empirical_p,
                "status":
                    row.status,
            }
        )

    j = joint.iloc[0]
    summary_rows.append(
        {
            "scope": "equal_lag_joint",
            "longitudinal_abundance_averaged_cross_cov":
                j["longitudinal_value"],
            "longitudinal_bootstrap_q025":
                j["longitudinal_bootstrap_q025"],
            "longitudinal_bootstrap_q975":
                j["longitudinal_bootstrap_q975"],
            "pseudo_n":
                j["pseudo_n_complete_configurations"],
            "pseudo_median":
                j["pseudo_median"],
            "pseudo_q025":
                j["pseudo_q025"],
            "pseudo_q975":
                j["pseudo_q975"],
            "empirical_p":
                j["empirical_p"],
            "status":
                j["status"],
        }
    )

    analysis_summary = pd.DataFrame(
        summary_rows
    )

    # -----------------------------------------------------------------
    # Write outputs.
    # -----------------------------------------------------------------
    support.to_csv(
        outdir / "01_shared_support_by_dt.csv",
        index=False,
    )
    stable_bins.to_csv(
        outdir / "07b_stable_native_bins_by_dt.csv",
        index=False,
    )

    if test_mode and not smoke_selection_audit.empty:
        smoke_selection_audit.to_csv(
            outdir / "00_smoke_selection_audit.csv",
            index=False,
        )

    pointwise.to_csv(
        outdir / "02_pointwise_cross_covariance_comparison.csv",
        index=False,
    )

    pseudo_metrics.to_csv(
        outdir / "03_pseudo_configuration_matched_support_metrics.csv",
        index=False,
    )

    tests.to_csv(
        outdir / "04_empirical_tests_by_dt.csv",
        index=False,
    )

    bootstrap.to_csv(
        outdir / "05_longitudinal_matched_support_bootstrap.csv",
        index=False,
    )

    joint.to_csv(
        outdir / "06_joint_comparison.csv",
        index=False,
    )

    audit.to_csv(
        outdir / "07_support_audit.csv",
        index=False,
    )

    native_global.to_csv(
        outdir / "08_native_global_by_dt.csv",
        index=False,
    )

    analysis_summary.to_csv(
        outdir / "00_analysis_summary.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Input manifest.
    # -----------------------------------------------------------------
    input_manifest = pd.DataFrame(
        [
            {
                "source": (
                    "Step11"
                    if key.startswith("step11_")
                    else "Step8"
                ),
                "input": key,
                "path": str(path),
                "exists": bool(path.exists()),
            }
            for key, path in inputs["paths"].items()
        ]
    )

    input_manifest.to_csv(
        outdir / "00_input_manifest.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Run config.
    # -----------------------------------------------------------------
    run_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "step11_dir": str(
            inputs["step11_dir"]
        ),
        "pseudo_step8_dir": str(
            inputs["step8_dir"]
        ),
        "outdir": str(outdir),
        "dt_values": dt_values,
        "primary_estimand": (
            "Cov(dx_observed_rep1, dx_observed_rep2 | "
            "xmid_latent, dt, common4)"
        ),
        "primary_support": "common4",
        "primary_conditioning": "xmid_latent",
        "operational_class_filter": None,
        "percentile_normalization": False,
        "matching": {
            "coordinate": "absolute xmid_latent",
            "support": (
                "stable Step-8 native bin centers inside valid Step-11 support"
            ),
            "pseudo_interpolation": False,
            "pseudo_imputation": False,
            "longitudinal_interpolation": (
                "adjacent valid Step-11 native bins only"
            ),
            "numerical_boundary_tolerance": float(
                SUPPORT_TOL
            ),
            "pseudo_stable_bin_min_valid_fraction":
                float(
                    args.min_pseudo_bin_valid_fraction
                ),
            "pseudo_stable_bin_fraction_denominator":
                "full Step-8 randomization ensemble",
            "min_shared_native_bins": int(
                args.min_shared_native_bins
            ),
        },
        "pseudo_null": {
            "unit": "Step-8 randomization configuration",
            "source_biological_subjects": 1,
            "configuration_is_biological_subject": False,
            "full_configuration_count": int(
                full_pseudo_n
            ),
            "selected_configuration_count": int(
                len(selected_ids)
            ),
            "min_matched_native_bin_coverage": float(
                args.min_null_bin_coverage
            ),
            "production_requires_complete_matched_native_bins": True,
            "min_configurations_for_inference": int(
                inferential_min_configurations
            ),
            "min_configurations_for_execution": int(
                execution_min_configurations
            ),
            "randomization_interval_not_biological_confidence_interval":
                True,
        },
        "longitudinal_bootstrap": {
            "unit": "biological subject",
            "n_bootstrap": int(
                n_bootstrap
            ),
            "seed": int(
                bootstrap_seed
            ),
            "same_subject_draw_shared_across_dt": True,
        },
        "primary_test": {
            "metric": (
                "equal-bin mean cross_cov over matched stable Step-8 native "
                "absolute xmid centers"
            ),
            "tail": "upper",
            "finite_sample_plus_one_correction": True,
        },
        "joint_summary": {
            "metric": (
                "equal-lag mean of matched-support "
                "abundance-averaged cross_cov"
            ),
            "requires_pseudo_eligibility_at_all_dt": True,
            "temporal_scaling_test": False,
        },
        "step11_script_version": step11_cfg.get(
            "script_version"
        ),
        "step11_analysis_signature": step11_cfg.get(
            "analysis_signature"
        ),
        "step8_script_version": step8_cfg.get(
            "script_version"
        ),
        "step8_analysis_signature": step8_cfg.get(
            "analysis_signature"
        ),
        "testing_subset_mode": bool(
            test_mode
        ),
        "smoke_test_selection": (
            "lag-wise support-aware union; same configuration not required at all dt"
            if test_mode
            else None
        ),
    }

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "step11_analysis_signature": step11_cfg.get(
            "analysis_signature"
        ),
        "step8_analysis_signature": step8_cfg.get(
            "analysis_signature"
        ),
        "step11_by_bin_dt": path_identity(
            inputs["paths"]["step11_by_bin_dt"]
        ),
        "step8_config_binned": path_identity(
            inputs["paths"]["step8_config_binned"]
        ),
        "step8_binned_envelope": path_identity(
            inputs["paths"]["step8_binned_envelope"]
        ),
        "selected_configuration_ids": list(selected_ids),
        "parameters": {
            "dt_values": [int(x) for x in dt_values],
            "min_pseudo_bin_valid_fraction": float(
                args.min_pseudo_bin_valid_fraction
            ),
            "min_shared_native_bins": int(
                args.min_shared_native_bins
            ),
            "min_null_bin_coverage": float(
                args.min_null_bin_coverage
            ),
            "min_null_configurations": int(
                inferential_min_configurations
            ),
            "n_bootstrap": int(n_bootstrap),
            "bootstrap_seed": int(bootstrap_seed),
        },
        "pseudo_missing_bin_interpolation": False,
        "pseudo_missing_bin_imputation": False,
    }
    analysis_signature = stable_sha256(signature_payload)
    run_config["analysis_signature"] = analysis_signature

    (
        outdir / "00_run_signature.json"
    ).write_text(
        json.dumps(
            {
                "analysis_signature": analysis_signature,
                "signature_payload": signature_payload,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    (
        outdir / "00_run_config.json"
    ).write_text(
        json.dumps(
            json_safe(run_config),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_readme(outdir)
    write_manifest(outdir)

    # -----------------------------------------------------------------
    # Console summary.
    # -----------------------------------------------------------------
    print("\n[DONE] Step 12 longitudinal vs pseudo cross-replicate fluctuations")

    for row in tests.itertuples(index=False):
        ptxt = (
            f"{row.empirical_p:.6g}"
            if np.isfinite(row.empirical_p)
            else "NA"
        )

        print(
            f"dt={int(row.dt)}: "
            f"real={row.longitudinal_value:.6f}, "
            f"pseudo median={row.pseudo_median:.6f}, "
            f"n_null={int(row.pseudo_n)}, "
            f"p={ptxt}"
        )

    jr = joint.iloc[0]
    jptxt = (
        f"{jr['empirical_p']:.6g}"
        if np.isfinite(jr["empirical_p"])
        else "NA"
    )

    print(
        "joint equal-lag:",
        f"real={jr['longitudinal_value']:.6f},",
        f"pseudo median={jr['pseudo_median']:.6f},",
        f"n_null={int(jr['pseudo_n_complete_configurations'])},",
        f"p={jptxt}",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
