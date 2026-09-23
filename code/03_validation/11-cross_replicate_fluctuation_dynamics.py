#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
11-cross_replicate_fluctuation_dynamics.py
===========================================

Canonical ClonoDynamics Step 11.

PURPOSE
-------
Quantify the replicate-consistent stochastic component of genuine longitudinal
TCRβ clonotype displacement using cross-replicate covariance.

The primary fluctuation estimand is

    V_cross(x, dt)
        = Cov(Delta x_A, Delta x_B | x_mid, dt, common4)

where

    Delta x_A = dx_observed_rep1
    Delta x_B = dx_observed_rep2

and the primary conditioning coordinate is

    x_mid = xmid_latent.

The support is `common4`: both technical replicates are observed at both
transition endpoints. TT/TF/FT/FF operational classes are NOT used as a primary
filter.

The within-replicate comparator is

    V_same
        = 0.5 * [Var(Delta x_A) + Var(Delta x_B)]

and

    replicate_specific_excess
        = V_same - V_cross.

`V_cross` remains signed and is never truncated at zero.

POSITION IN THE PIPELINE
------------------------
Step 7 validated replicate-decoupled forward estimands in the pseudo technical
null.
Step 8 validated cross-replicate covariance on common4 and established
xmid_latent as the primary conditioning coordinate.
Steps 9-10 applied and benchmarked the forward estimator in the genuine
longitudinal cohort.
Step 11 now applies the fluctuation estimand to the genuine longitudinal cohort.

INPUT
-----
Current Step-5 generic transition table:

    longitudinal_transitions.parquet

Required fields:

    subject
    t0
    t1
    dt
    common4
    xmid_latent
    dx_observed_rep1
    dx_observed_rep2

No latent displacement is used as the primary fluctuation outcome.

ABUNDANCE GRID
--------------
One fixed abundance grid is used across dt = 1-5.

Default:
    n_bins  = 20
    binning = equal_width
    q_low   = 0.01
    q_high  = 0.99

The default grid is derived from xmid_latent among all finite common4
transitions across the selected dt values. For equal-width binning, q_low and
q_high define the robust abundance-range endpoints.

An external bin-edge file may be supplied with --bin-edges-file.

PRIMARY CELL STATISTICS
-----------------------
For every dt x abundance bin the script reports:

    n
    n_subjects
    n_subject_intervals

    mean_dx_A
    mean_dx_B

    var_dx_A
    var_dx_B
    same_var_mean

    cross_cov

    replicate_specific_excess
        = same_var_mean - cross_cov

    shared_fraction
        = cross_cov / same_var_mean
        when same_var_mean > 0

The covariance / variance convention is ddof = 1.

COVARIANCE IDENTITY AUDIT
-------------------------
The script verifies numerically that

    Cov(A, B)
      = Var[(A+B)/2] - 0.25 Var(A-B)

for every valid pooled cell.

The residual is written as:

    identity_residual
        = cross_cov - identity_rhs.

SUBJECT-CLUSTER BOOTSTRAP
-------------------------
Uncertainty in pooled dt x abundance-bin estimates is quantified by resampling
BIOLOGICAL SUBJECTS with replacement.

Default:
    n_bootstrap   = 2000
    bootstrap_seed = 123

The pooled estimator is reconstructed exactly from compact subject-level
sufficient statistics. Thus the bootstrap preserves the transition-weighted
cohort estimand while accounting for clustering by biological subject.

Bootstrap 95% intervals are reported for:

    cross_cov
    same_var_mean
    replicate_specific_excess
    shared_fraction.

A bootstrap draw contributes to a dt x abundance-bin interval only when the
reconstructed cell satisfies the SAME `--min-n` threshold used by the point
estimate. Thus uncertainty is evaluated on the fixed point-estimate support
rather than allowing under-supported bootstrap cells to enter the interval.

SUBJECT-LEVEL OUTPUT
--------------------
For downstream temporal-scaling analysis, Step 11 also writes exact
subject x dt x abundance-bin covariance decompositions.

These subject-level values are NOT averaged to define the primary Step-11
cohort curve; they are retained for later subject-resolved temporal analyses.

INTERVAL-LEVEL OUTPUT
---------------------
The script additionally writes subject x physical-interval x abundance-bin
metrics:

    subject
    t0
    t1
    dt
    bin
    ...

This table provides the interval-position information needed by the later
calendar/common-start/common-end sensitivity analyses.

OUTPUT
------
<outdir>/
    00_run_config.json
    00_pipeline_manifest.csv
    README_outputs.md

    01_fluctuation_bin_edges.csv
    02_fluctuation_by_bin_dt.csv
    03_global_fluctuation_by_dt.csv
    04_subject_fluctuation_by_bin_dt.csv
    05_subject_interval_fluctuation_by_bin.csv
    06_support_by_dt.csv
    07_covariance_identity_audit.csv
    07b_cell_support_audit.csv

    08_subject_bin_sufficient.parquet
    09_subject_dt_sufficient.parquet

Optional:
    10_bootstrap_draws.parquet
        only when --keep-bootstrap-draws is supplied.

PRIMARY RUN
-----------
python3 11-cross_replicate_fluctuation_dynamics.py \
    --transitions \
    ./results_4/5-longitudinal_transition_assembly/longitudinal_transitions.parquet \
    --outdir \
    ./results_4/11-cross_replicate_fluctuation_dynamics \
    --dt-values 1 2 3 4 5 \
    --n-bins 20 \
    --binning equal_width \
    --q-low 0.01 \
    --q-high 0.99 \
    --min-n 50 \
    --n-bootstrap 2000 \
    --bootstrap-seed 123

INTERPRETATION
--------------
Step 11 asks whether genuine longitudinal clonotype displacement contains a
replicate-consistent fluctuation component after requiring all four observed
endpoint measurements.

It characterizes shared fluctuation covariance as a function of abundance and
time lag.

It does NOT yet compare the real cohort with the pseudo technical null; that
comparison belongs to Step 12.

It does NOT yet test temporal accumulation over a fixed abundance core; that
belongs to Step 13.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    import polars as pl
except Exception as exc:
    raise SystemExit(
        "ERROR: polars is required. Activate the ClonoDynamics environment."
    ) from exc


SCRIPT_VERSION = "v2-fixed-support-bootstrap-signature-2026-09-17"

DDOF = 1

REQUIRED_COLUMNS = [
    "subject",
    "t0",
    "t1",
    "dt",
    "common4",
    "xmid_latent",
    "dx_observed_rep1",
    "dx_observed_rep2",
]


# =============================================================================
# Generic helpers
# =============================================================================

def eprint(*args: Any, **kwargs: Any) -> None:
    import sys
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def collect_streaming(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)


def schema_names(path: Path) -> List[str]:
    return list(pl.scan_parquet(path).collect_schema().names())


def parse_bool_expr(column: str) -> pl.Expr:
    s = pl.col(column)
    return (
        pl.when(s.cast(pl.Boolean, strict=False).is_not_null())
        .then(s.cast(pl.Boolean, strict=False).fill_null(False))
        .otherwise(
            s.cast(pl.String, strict=False)
            .str.strip_chars()
            .str.to_lowercase()
            .is_in(["true", "t", "1", "yes", "y"])
            .fill_null(False)
        )
    )


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


def write_manifest(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "relative_path": str(path.relative_to(outdir)),
                    "suffix": path.suffix.lower(),
                    "size_bytes": int(path.stat().st_size),
                }
            )
    pl.DataFrame(rows).write_csv(
        outdir / "00_pipeline_manifest.csv"
    )


# =============================================================================
# Binning
# =============================================================================

def read_edges_file(path: Path) -> np.ndarray:
    path = Path(path)

    if path.suffix.lower() == ".csv":
        d = pl.read_csv(path)
    elif path.suffix.lower() in {".parquet", ".pq"}:
        d = pl.read_parquet(path)
    else:
        values = [
            float(line.strip())
            for line in path.read_text().splitlines()
            if line.strip()
        ]
        edges = np.unique(np.asarray(values, dtype=float))
        if edges.size < 3:
            raise ValueError("External bin-edge file contains too few edges.")
        return edges

    if "edge" in d.columns:
        edges = np.asarray(d["edge"].to_numpy(), dtype=float)
    elif {"x_left", "x_right"}.issubset(set(d.columns)):
        edges = np.r_[
            np.asarray(d["x_left"].to_numpy(), dtype=float)[:1],
            np.asarray(d["x_right"].to_numpy(), dtype=float),
        ]
    elif {"x_lo", "x_hi"}.issubset(set(d.columns)):
        edges = np.r_[
            np.asarray(d["x_lo"].to_numpy(), dtype=float)[:1],
            np.asarray(d["x_hi"].to_numpy(), dtype=float),
        ]
    else:
        raise ValueError(
            f"{path}: expected edge or x_left/x_right columns."
        )

    edges = np.unique(edges[np.isfinite(edges)])
    edges.sort()
    if edges.size < 3:
        raise ValueError("External bin-edge file contains too few finite edges.")
    return edges


def filtered_primary_lf(
    transitions: Path,
    dt_values: Sequence[int],
) -> pl.LazyFrame:
    return (
        pl.scan_parquet(transitions)
        .filter(
            pl.col("dt")
            .cast(pl.Int64, strict=False)
            .is_in([int(x) for x in dt_values])
        )
        .select(
            [
                pl.col("subject").cast(pl.Int64, strict=False),
                pl.col("t0").cast(pl.Int64, strict=False),
                pl.col("t1").cast(pl.Int64, strict=False),
                pl.col("dt").cast(pl.Int64, strict=False),
                parse_bool_expr("common4").alias("common4"),
                pl.col("xmid_latent")
                .cast(pl.Float64, strict=False)
                .alias("x"),
                pl.col("dx_observed_rep1")
                .cast(pl.Float64, strict=False)
                .alias("A"),
                pl.col("dx_observed_rep2")
                .cast(pl.Float64, strict=False)
                .alias("B"),
            ]
        )
        .filter(
            pl.col("common4")
            & pl.col("x").is_finite()
            & pl.col("A").is_finite()
            & pl.col("B").is_finite()
        )
    )


def estimate_edges(
    transitions: Path,
    dt_values: Sequence[int],
    n_bins: int,
    mode: str,
    q_low: float,
    q_high: float,
) -> np.ndarray:
    xlf = filtered_primary_lf(
        transitions,
        dt_values,
    ).select("x")

    if mode == "equal_width":
        q = collect_streaming(
            xlf.select(
                [
                    pl.col("x")
                    .quantile(float(q_low), interpolation="linear")
                    .alias("lo"),
                    pl.col("x")
                    .quantile(float(q_high), interpolation="linear")
                    .alias("hi"),
                ]
            )
        )
        lo = float(q["lo"][0])
        hi = float(q["hi"][0])
        edges = np.linspace(
            lo,
            hi,
            int(n_bins) + 1,
        )

    elif mode == "quantile":
        probs = np.linspace(
            float(q_low),
            float(q_high),
            int(n_bins) + 1,
        )
        exprs = [
            pl.col("x")
            .quantile(float(prob), interpolation="linear")
            .alias(f"q_{i:04d}")
            for i, prob in enumerate(probs)
        ]
        row = collect_streaming(
            xlf.select(exprs)
        ).row(0)
        edges = np.asarray(row, dtype=float)

    else:
        raise ValueError(f"Unknown binning mode: {mode}")

    edges = np.unique(edges[np.isfinite(edges)])
    edges.sort()
    if edges.size < 3:
        raise RuntimeError("Could not derive at least two abundance bins.")
    return edges


def save_edges(edges: np.ndarray, path: Path) -> None:
    pl.DataFrame(
        {
            "bin": np.arange(
                len(edges) - 1,
                dtype=np.int32,
            ),
            "x_left": edges[:-1],
            "x_right": edges[1:],
            "x_center": 0.5 * (
                edges[:-1] + edges[1:]
            ),
        }
    ).write_csv(path)


def bin_expr(
    x_col: str,
    edges: np.ndarray,
) -> pl.Expr:
    x = pl.col(x_col).cast(pl.Float64, strict=False)
    expr = pl.lit(None, dtype=pl.Int32)

    for i in range(len(edges) - 2, -1, -1):
        lo = float(edges[i])
        hi = float(edges[i + 1])

        if i == len(edges) - 2:
            cond = (x >= lo) & (x <= hi)
        else:
            cond = (x >= lo) & (x < hi)

        expr = (
            pl.when(cond)
            .then(pl.lit(i, dtype=pl.Int32))
            .otherwise(expr)
        )

    return expr.alias("bin")


# =============================================================================
# Sufficient statistics
# =============================================================================

SUFF_COLUMNS = [
    "n",
    "sum_A",
    "sum_B",
    "sum_A2",
    "sum_B2",
    "sum_AB",
]


def sufficient_agg_exprs() -> List[pl.Expr]:
    return [
        pl.len().alias("n"),
        pl.col("A").sum().alias("sum_A"),
        pl.col("B").sum().alias("sum_B"),
        (pl.col("A") * pl.col("A")).sum().alias("sum_A2"),
        (pl.col("B") * pl.col("B")).sum().alias("sum_B2"),
        (pl.col("A") * pl.col("B")).sum().alias("sum_AB"),
    ]


def build_sufficient_tables(
    transitions: Path,
    dt_values: Sequence[int],
    edges: np.ndarray,
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    base = (
        filtered_primary_lf(
            transitions,
            dt_values,
        )
        .with_columns(bin_expr("x", edges))
        .filter(pl.col("bin").is_not_null())
    )

    subject_bin = collect_streaming(
        base.group_by(
            ["subject", "dt", "bin"]
        )
        .agg(sufficient_agg_exprs())
        .sort(["subject", "dt", "bin"])
    )

    interval_bin = collect_streaming(
        base.group_by(
            ["subject", "t0", "t1", "dt", "bin"]
        )
        .agg(sufficient_agg_exprs())
        .sort(["subject", "t0", "t1", "dt", "bin"])
    )

    subject_dt = collect_streaming(
        base.group_by(
            ["subject", "dt"]
        )
        .agg(sufficient_agg_exprs())
        .sort(["subject", "dt"])
    )

    return subject_bin, interval_bin, subject_dt


# =============================================================================
# Statistical reconstruction
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

    if (
        not np.isfinite(n)
        or n <= DDOF
    ):
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

    mean_A = float(sum_A / n)
    mean_B = float(sum_B / n)

    denom = n - DDOF

    var_A = (
        float(sum_A2)
        - (float(sum_A) ** 2) / n
    ) / denom

    var_B = (
        float(sum_B2)
        - (float(sum_B) ** 2) / n
    ) / denom

    cross_cov = (
        float(sum_AB)
        - float(sum_A) * float(sum_B) / n
    ) / denom

    # Numerical tolerance only; no clipping of cross_cov.
    if var_A < 0 and abs(var_A) < 1e-12:
        var_A = 0.0
    if var_B < 0 and abs(var_B) < 1e-12:
        var_B = 0.0

    same_var_mean = 0.5 * (
        var_A + var_B
    )
    excess = same_var_mean - cross_cov

    shared_fraction = (
        cross_cov / same_var_mean
        if np.isfinite(same_var_mean)
        and same_var_mean > 0
        else np.nan
    )

    # Cov(A,B) = Var((A+B)/2) - 0.25 Var(A-B)
    var_mean_pair = 0.25 * (
        var_A + var_B + 2.0 * cross_cov
    )
    var_difference = (
        var_A + var_B - 2.0 * cross_cov
    )
    identity_rhs = (
        var_mean_pair
        - 0.25 * var_difference
    )
    identity_residual = (
        cross_cov - identity_rhs
    )

    return {
        "n": n,
        "mean_dx_A": mean_A,
        "mean_dx_B": mean_B,
        "var_dx_A": float(var_A),
        "var_dx_B": float(var_B),
        "cross_cov": float(cross_cov),
        "same_var_mean": float(same_var_mean),
        "replicate_specific_excess": float(excess),
        "shared_fraction": float(shared_fraction),
        "identity_rhs": float(identity_rhs),
        "identity_residual": float(identity_residual),
    }


def aggregate_sufficient_rows(
    frame: pl.DataFrame,
    group_cols: Sequence[str],
) -> pl.DataFrame:
    return (
        frame.lazy()
        .group_by(list(group_cols))
        .agg(
            [
                pl.col(c).sum().alias(c)
                for c in SUFF_COLUMNS
            ]
        )
        .sort(list(group_cols))
        .collect()
    )


def metrics_table_from_sufficient(
    frame: pl.DataFrame,
    group_cols: Sequence[str],
    extra_counts: Optional[pl.DataFrame] = None,
) -> pl.DataFrame:
    rows: List[Dict[str, object]] = []

    for row in frame.iter_rows(named=True):
        metrics = metrics_from_sums(
            row["n"],
            row["sum_A"],
            row["sum_B"],
            row["sum_A2"],
            row["sum_B2"],
            row["sum_AB"],
        )
        out = {
            c: row[c]
            for c in group_cols
        }
        out.update(metrics)
        rows.append(out)

    result = (
        pl.DataFrame(rows)
        if rows
        else pl.DataFrame()
    )

    if (
        extra_counts is not None
        and not result.is_empty()
    ):
        result = result.join(
            extra_counts,
            on=list(group_cols),
            how="left",
        )

    return result


# =============================================================================
# Subject bootstrap
# =============================================================================

def subject_bin_arrays(
    subject_bin: pl.DataFrame,
    subjects: Sequence[int],
    dt_values: Sequence[int],
    n_bins: int,
) -> Dict[str, np.ndarray]:
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

    arrays = {
        c: np.zeros(
            shape,
            dtype=float,
        )
        for c in SUFF_COLUMNS
    }

    for row in subject_bin.iter_rows(named=True):
        s = s_index[int(row["subject"])]
        d = d_index[int(row["dt"])]
        b = int(row["bin"])

        for c in SUFF_COLUMNS:
            arrays[c][s, d, b] = float(
                row[c]
            )

    return arrays


def bootstrap_subject_cells(
    subject_bin: pl.DataFrame,
    subjects: Sequence[int],
    dt_values: Sequence[int],
    n_bins: int,
    n_bootstrap: int,
    seed: int,
) -> pl.DataFrame:
    if int(n_bootstrap) <= 0:
        return pl.DataFrame()

    arrays = subject_bin_arrays(
        subject_bin,
        subjects,
        dt_values,
        n_bins,
    )

    rng = np.random.default_rng(
        int(seed)
    )

    records: List[Dict[str, object]] = []

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

        sums = {
            c: np.tensordot(
                mult,
                arr,
                axes=(0, 0),
            )
            for c, arr in arrays.items()
        }

        for di, dt in enumerate(dt_values):
            for b in range(int(n_bins)):
                m = metrics_from_sums(
                    sums["n"][di, b],
                    sums["sum_A"][di, b],
                    sums["sum_B"][di, b],
                    sums["sum_A2"][di, b],
                    sums["sum_B2"][di, b],
                    sums["sum_AB"][di, b],
                )

                records.append(
                    {
                        "bootstrap": int(boot),
                        "dt": int(dt),
                        "bin": int(b),
                        "n": m["n"],
                        "cross_cov": m["cross_cov"],
                        "same_var_mean": m["same_var_mean"],
                        "replicate_specific_excess": m[
                            "replicate_specific_excess"
                        ],
                        "shared_fraction": m["shared_fraction"],
                    }
                )

    return pl.DataFrame(records)


def bootstrap_ci_by_cell(
    boot: pl.DataFrame,
    min_n: int,
) -> pl.DataFrame:
    if boot.is_empty():
        return pl.DataFrame()

    metrics = [
        "cross_cov",
        "same_var_mean",
        "replicate_specific_excess",
        "shared_fraction",
    ]

    aggs: List[pl.Expr] = []

    for metric in metrics:
        aggs.extend(
            [
                pl.col(metric)
                .drop_nans()
                .quantile(
                    0.025,
                    interpolation="linear",
                )
                .alias(
                    f"{metric}_ci025"
                ),
                pl.col(metric)
                .drop_nans()
                .quantile(
                    0.975,
                    interpolation="linear",
                )
                .alias(
                    f"{metric}_ci975"
                ),
            ]
        )

    aggs.append(
        pl.col("cross_cov")
        .is_finite()
        .sum()
        .alias("n_bootstrap_cross_cov_finite")
    )

    return (
        boot.lazy()
        .with_columns(
            (
                pl.col("n") >= int(min_n)
            ).alias("bootstrap_meets_min_n")
        )
        .with_columns(
            [
                pl.when(pl.col("bootstrap_meets_min_n"))
                .then(pl.col(metric))
                .otherwise(None)
                .alias(metric)
                for metric in metrics
            ]
        )
        .group_by(["dt", "bin"])
        .agg(
            aggs
            + [
                pl.col("bootstrap_meets_min_n")
                .sum()
                .alias("n_bootstrap_meets_min_n"),
                pl.len().alias("n_bootstrap_total"),
            ]
        )
        .with_columns(
            (
                pl.col("n_bootstrap_meets_min_n")
                / pl.col("n_bootstrap_total")
            ).alias("bootstrap_valid_fraction")
        )
        .sort(["dt", "bin"])
        .collect()
    )


def bootstrap_global_dt(
    subject_dt: pl.DataFrame,
    subjects: Sequence[int],
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
    min_n: int,
) -> pl.DataFrame:
    """
    Subject-cluster bootstrap for global (abundance-pooled) dt metrics.
    """
    if int(n_bootstrap) <= 0:
        return pl.DataFrame()

    s_index = {int(s): i for i, s in enumerate(subjects)}
    d_index = {int(dt): i for i, dt in enumerate(dt_values)}

    arrays = {
        c: np.zeros(
            (len(subjects), len(dt_values)),
            dtype=float,
        )
        for c in SUFF_COLUMNS
    }

    for row in subject_dt.iter_rows(named=True):
        si = s_index[int(row["subject"])]
        di = d_index[int(row["dt"])]
        for c in SUFF_COLUMNS:
            arrays[c][si, di] = float(row[c])

    rng = np.random.default_rng(int(seed))
    records: List[Dict[str, object]] = []

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

        sums = {
            c: np.tensordot(
                mult,
                arr,
                axes=(0, 0),
            )
            for c, arr in arrays.items()
        }

        for di, dt in enumerate(dt_values):
            m = metrics_from_sums(
                sums["n"][di],
                sums["sum_A"][di],
                sums["sum_B"][di],
                sums["sum_A2"][di],
                sums["sum_B2"][di],
                sums["sum_AB"][di],
            )
            records.append(
                {
                    "bootstrap": int(boot),
                    "dt": int(dt),
                    "n": m["n"],
                    "cross_cov": m["cross_cov"],
                    "same_var_mean": m["same_var_mean"],
                    "replicate_specific_excess":
                        m["replicate_specific_excess"],
                    "shared_fraction": m["shared_fraction"],
                    "bootstrap_meets_min_n": bool(
                        np.isfinite(m["n"])
                        and m["n"] >= int(min_n)
                    ),
                }
            )

    return pl.DataFrame(records)


def bootstrap_global_dt_ci(
    boot: pl.DataFrame,
) -> pl.DataFrame:
    if boot.is_empty():
        return pl.DataFrame()

    metrics = [
        "cross_cov",
        "same_var_mean",
        "replicate_specific_excess",
        "shared_fraction",
    ]
    aggs: List[pl.Expr] = []

    for metric in metrics:
        aggs.extend(
            [
                pl.when(pl.col("bootstrap_meets_min_n"))
                .then(pl.col(metric))
                .otherwise(None)
                .drop_nulls()
                .drop_nans()
                .quantile(0.025, interpolation="linear")
                .alias(f"{metric}_ci025"),
                pl.when(pl.col("bootstrap_meets_min_n"))
                .then(pl.col(metric))
                .otherwise(None)
                .drop_nulls()
                .drop_nans()
                .quantile(0.975, interpolation="linear")
                .alias(f"{metric}_ci975"),
            ]
        )

    return (
        boot.lazy()
        .group_by("dt")
        .agg(
            aggs
            + [
                pl.col("bootstrap_meets_min_n")
                .sum()
                .alias("n_bootstrap_meets_min_n"),
                pl.len().alias("n_bootstrap_total"),
            ]
        )
        .with_columns(
            (
                pl.col("n_bootstrap_meets_min_n")
                / pl.col("n_bootstrap_total")
            ).alias("bootstrap_valid_fraction")
        )
        .sort("dt")
        .collect()
    )


# =============================================================================
# Pooled / subject / interval outputs
# =============================================================================

def pooled_by_bin_dt(
    subject_bin: pl.DataFrame,
    interval_bin: pl.DataFrame,
    edges: np.ndarray,
    min_n: int,
    bootstrap_ci: pl.DataFrame,
) -> pl.DataFrame:
    pooled_suff = aggregate_sufficient_rows(
        subject_bin,
        ["dt", "bin"],
    )

    subject_counts = (
        subject_bin.lazy()
        .group_by(["dt", "bin"])
        .agg(
            pl.col("subject")
            .n_unique()
            .alias("n_subjects")
        )
        .collect()
    )

    interval_counts = (
        interval_bin.lazy()
        .group_by(["dt", "bin"])
        .agg(
            pl.struct(
                ["subject", "t0", "t1"]
            )
            .n_unique()
            .alias("n_subject_intervals")
        )
        .collect()
    )

    extra = subject_counts.join(
        interval_counts,
        on=["dt", "bin"],
        how="outer",
        coalesce=True,
    )

    metrics = metrics_table_from_sufficient(
        pooled_suff,
        ["dt", "bin"],
        extra_counts=extra,
    ).with_columns(
        (
            pl.col("n") >= int(min_n)
        ).alias("meets_min_n")
    )

    edge_df = pl.DataFrame(
        {
            "bin": np.arange(
                len(edges) - 1,
                dtype=np.int32,
            ),
            "x_left": edges[:-1],
            "x_right": edges[1:],
            "x_center": 0.5 * (
                edges[:-1] + edges[1:]
            ),
        }
    )

    metrics = metrics.join(
        edge_df,
        on="bin",
        how="left",
    )

    if not bootstrap_ci.is_empty():
        metrics = metrics.join(
            bootstrap_ci,
            on=["dt", "bin"],
            how="left",
        )

    return metrics.sort(["dt", "bin"])


def subject_metrics_by_bin_dt(
    subject_bin: pl.DataFrame,
    edges: np.ndarray,
) -> pl.DataFrame:
    metrics = metrics_table_from_sufficient(
        subject_bin,
        ["subject", "dt", "bin"],
    )

    edge_df = pl.DataFrame(
        {
            "bin": np.arange(
                len(edges) - 1,
                dtype=np.int32,
            ),
            "x_left": edges[:-1],
            "x_right": edges[1:],
            "x_center": 0.5 * (
                edges[:-1] + edges[1:]
            ),
        }
    )

    return (
        metrics.join(
            edge_df,
            on="bin",
            how="left",
        )
        .sort(["subject", "dt", "bin"])
    )


def interval_metrics_by_bin(
    interval_bin: pl.DataFrame,
    edges: np.ndarray,
) -> pl.DataFrame:
    metrics = metrics_table_from_sufficient(
        interval_bin,
        ["subject", "t0", "t1", "dt", "bin"],
    )

    edge_df = pl.DataFrame(
        {
            "bin": np.arange(
                len(edges) - 1,
                dtype=np.int32,
            ),
            "x_left": edges[:-1],
            "x_right": edges[1:],
            "x_center": 0.5 * (
                edges[:-1] + edges[1:]
            ),
        }
    )

    return (
        metrics.join(
            edge_df,
            on="bin",
            how="left",
        )
        .sort(
            ["subject", "t0", "t1", "bin"]
        )
    )


def global_by_dt(
    subject_dt: pl.DataFrame,
    bootstrap_ci: Optional[pl.DataFrame] = None,
) -> pl.DataFrame:
    pooled = aggregate_sufficient_rows(
        subject_dt,
        ["dt"],
    )

    counts = (
        subject_dt.lazy()
        .group_by("dt")
        .agg(
            pl.col("subject")
            .n_unique()
            .alias("n_subjects")
        )
        .collect()
    )

    result = (
        metrics_table_from_sufficient(
            pooled,
            ["dt"],
            extra_counts=counts,
        )
        .sort("dt")
    )

    if (
        bootstrap_ci is not None
        and not bootstrap_ci.is_empty()
    ):
        result = result.join(
            bootstrap_ci,
            on="dt",
            how="left",
        )

    return result.sort("dt")


def support_by_dt(
    transitions: Path,
    dt_values: Sequence[int],
) -> pl.DataFrame:
    lf = (
        pl.scan_parquet(transitions)
        .filter(
            pl.col("dt")
            .cast(pl.Int64, strict=False)
            .is_in([int(x) for x in dt_values])
        )
        .select(
            [
                pl.col("subject")
                .cast(pl.Int64, strict=False),
                pl.col("t0")
                .cast(pl.Int64, strict=False),
                pl.col("t1")
                .cast(pl.Int64, strict=False),
                pl.col("dt")
                .cast(pl.Int64, strict=False),
                parse_bool_expr("common4")
                .alias("common4"),
                pl.col("xmid_latent")
                .cast(pl.Float64, strict=False)
                .alias("x"),
                pl.col("dx_observed_rep1")
                .cast(pl.Float64, strict=False)
                .alias("A"),
                pl.col("dx_observed_rep2")
                .cast(pl.Float64, strict=False)
                .alias("B"),
            ]
        )
        .with_columns(
            (
                pl.col("common4")
                & pl.col("x").is_finite()
                & pl.col("A").is_finite()
                & pl.col("B").is_finite()
            ).alias("primary_complete_case")
        )
        .group_by("dt")
        .agg(
            [
                pl.len().alias("n_all_transitions"),
                pl.col("common4")
                .sum()
                .alias("n_common4"),
                pl.col("primary_complete_case")
                .sum()
                .alias("n_primary_complete_case"),
                pl.col("subject")
                .n_unique()
                .alias("n_subjects"),
                pl.struct(
                    ["subject", "t0", "t1"]
                )
                .n_unique()
                .alias("n_subject_intervals"),
            ]
        )
        .with_columns(
            [
                (
                    pl.col("n_common4")
                    / pl.col("n_all_transitions")
                ).alias("fraction_common4"),
                (
                    pl.col("n_primary_complete_case")
                    / pl.col("n_all_transitions")
                ).alias(
                    "fraction_primary_complete_case"
                ),
            ]
        )
        .sort("dt")
    )

    return collect_streaming(lf)


def cell_support_audit(
    pooled: pl.DataFrame,
    dt_values: Sequence[int],
    min_n: int,
) -> pl.DataFrame:
    """
    Compact audit of abundance-bin support by lag.
    """
    rows: List[Dict[str, object]] = []

    for dt in dt_values:
        g = pooled.filter(pl.col("dt") == int(dt))
        if g.is_empty():
            rows.append(
                {
                    "dt": int(dt),
                    "n_bins_present": 0,
                    "n_bins_meets_min_n": 0,
                    "fraction_bins_meets_min_n": 0.0,
                    "min_n": np.nan,
                    "median_n": np.nan,
                    "max_n": np.nan,
                    "min_subjects": np.nan,
                    "median_subjects": np.nan,
                    "max_subjects": np.nan,
                    "min_bootstrap_valid_fraction": np.nan,
                }
            )
            continue

        n = g["n"].to_numpy().astype(float)
        subjects = (
            g["n_subjects"].to_numpy().astype(float)
            if "n_subjects" in g.columns
            else np.full(g.height, np.nan)
        )
        meets = (
            g["meets_min_n"].to_numpy().astype(bool)
            if "meets_min_n" in g.columns
            else n >= int(min_n)
        )
        boot_frac = (
            g["bootstrap_valid_fraction"].to_numpy().astype(float)
            if "bootstrap_valid_fraction" in g.columns
            else np.full(g.height, np.nan)
        )

        finite_boot = boot_frac[np.isfinite(boot_frac)]

        rows.append(
            {
                "dt": int(dt),
                "n_bins_present": int(g.height),
                "n_bins_meets_min_n": int(np.sum(meets)),
                "fraction_bins_meets_min_n": float(
                    np.mean(meets)
                ),
                "min_n": float(np.nanmin(n)),
                "median_n": float(np.nanmedian(n)),
                "max_n": float(np.nanmax(n)),
                "min_subjects": float(np.nanmin(subjects)),
                "median_subjects": float(np.nanmedian(subjects)),
                "max_subjects": float(np.nanmax(subjects)),
                "min_bootstrap_valid_fraction": (
                    float(np.min(finite_boot))
                    if finite_boot.size else np.nan
                ),
            }
        )

    return pl.DataFrame(rows)


# =============================================================================
# Audit
# =============================================================================

def identity_audit(
    pooled: pl.DataFrame,
    subject_metrics: pl.DataFrame,
    interval_metrics: pl.DataFrame,
) -> pl.DataFrame:
    rows = []

    for level, frame in [
        ("pooled_dt_bin", pooled),
        ("subject_dt_bin", subject_metrics),
        ("subject_interval_bin", interval_metrics),
    ]:
        if frame.is_empty():
            continue

        residual = (
            frame["identity_residual"]
            .to_numpy()
            .astype(float)
        )
        residual = residual[
            np.isfinite(residual)
        ]

        if residual.size == 0:
            continue

        rows.append(
            {
                "level": level,
                "n_cells": int(residual.size),
                "median_abs_identity_residual": float(
                    np.median(np.abs(residual))
                ),
                "q975_abs_identity_residual": float(
                    np.quantile(
                        np.abs(residual),
                        0.975,
                    )
                ),
                "max_abs_identity_residual": float(
                    np.max(np.abs(residual))
                ),
            }
        )

    return pl.DataFrame(rows)


# =============================================================================
# README
# =============================================================================

def write_readme(outdir: Path) -> None:
    text = """# Step 11 — real cross-replicate fluctuation dynamics

Primary estimand
----------------
V_cross = Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt, common4)

Support
-------
common4 only: both observed technical replicates are positive at both endpoints.

No TT/TF/FT/FF operational-class restriction is applied.

Comparator
----------
same_var_mean = 0.5 * [Var(dx_observed_rep1) + Var(dx_observed_rep2)]

replicate_specific_excess = same_var_mean - cross_cov

cross_cov is signed and is never truncated at zero.

Abundance conditioning
----------------------
One fixed xmid_latent abundance grid is reused across dt=1-5.

Uncertainty
-----------
Pooled dt x bin uncertainty is estimated with a biological subject-cluster
bootstrap. The cohort estimator remains transition-weighted inside each
bootstrap sample.

A bootstrap draw contributes to a cell interval only if its reconstructed
transition count satisfies the same min-n threshold used by the point estimate.
Global-by-lag metrics are bootstrapped with the same biological-subject
resampling scheme.

Downstream support
------------------
04_subject_fluctuation_by_bin_dt.csv is intended for subject-resolved temporal
scaling in Step 13.

05_subject_interval_fluctuation_by_bin.csv is intended for interval-position,
common-start and common-end controls in Step 14.

Step 11 does not compare the genuine cohort with the pseudo technical null.
That comparison belongs to Step 12.
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
            "Step 11: real cross-replicate fluctuation dynamics on common4."
        ),
    )

    p.add_argument(
        "--transitions",
        required=True,
        type=Path,
        help="Current Step-5 longitudinal_transitions.parquet.",
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--dataset-label",
        default="healthy_longitudinal",
    )

    p.add_argument(
        "--dt-values",
        type=int,
        nargs="+",
        default=[1, 2, 3, 4, 5],
    )

    p.add_argument(
        "--n-bins",
        type=int,
        default=20,
    )
    p.add_argument(
        "--binning",
        choices=["equal_width", "quantile"],
        default="equal_width",
    )
    p.add_argument(
        "--q-low",
        type=float,
        default=0.01,
    )
    p.add_argument(
        "--q-high",
        type=float,
        default=0.99,
    )
    p.add_argument(
        "--bin-edges-file",
        type=Path,
        default=None,
    )

    p.add_argument(
        "--min-n",
        type=int,
        default=50,
        help="Minimum pooled transitions per dt x abundance bin.",
    )

    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=2000,
    )
    p.add_argument(
        "--bootstrap-seed",
        type=int,
        default=123,
    )
    p.add_argument(
        "--keep-bootstrap-draws",
        action="store_true",
    )

    return p


def validate_args(args: argparse.Namespace) -> None:
    dt_values = sorted(set(int(x) for x in args.dt_values))

    if any(x < 1 for x in dt_values):
        raise ValueError("--dt-values must all be >=1")
    if int(args.n_bins) < 3:
        raise ValueError("--n-bins must be >=3")
    if not (
        0 <= float(args.q_low)
        < float(args.q_high)
        <= 1
    ):
        raise ValueError("Require 0 <= q-low < q-high <= 1")
    if int(args.min_n) < 2:
        raise ValueError("--min-n must be >=2")
    if int(args.n_bootstrap) < 0:
        raise ValueError("--n-bootstrap must be >=0")


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)

    transitions = (
        args.transitions
        .expanduser()
        .resolve(strict=True)
    )
    outdir = ensure_dir(
        args.outdir
        .expanduser()
        .resolve()
    )

    dt_values = sorted(
        set(int(x) for x in args.dt_values)
    )

    available = set(
        schema_names(transitions)
    )
    missing = [
        c for c in REQUIRED_COLUMNS
        if c not in available
    ]
    if missing:
        raise ValueError(
            "Current Step-5 transition table is missing required Step-11 "
            f"columns: {missing}"
        )


    present_dt = set(
        collect_streaming(
            pl.scan_parquet(transitions)
            .select(
                pl.col("dt")
                .cast(pl.Int64, strict=False)
                .unique()
            )
        )["dt"]
        .drop_nulls()
        .to_list()
    )
    absent_dt = [
        int(dt)
        for dt in dt_values
        if int(dt) not in present_dt
    ]
    if absent_dt:
        raise ValueError(
            f"Requested Step-11 lag values absent from transitions: {absent_dt}"
        )

    print("=== ClonoDynamics Step 11 ===")
    print("script version :", SCRIPT_VERSION)
    print("dataset        :", args.dataset_label)
    print("transitions    :", transitions)
    print("outdir         :", outdir)
    print("dt values      :", dt_values)
    print("primary support: common4")
    print("conditioning   : xmid_latent")
    print("outcomes       : dx_observed_rep1 / dx_observed_rep2")
    print("cross_cov signed / untruncated: YES")

    # ---------------------------------------------------------------
    # Fixed abundance grid.
    # ---------------------------------------------------------------
    if args.bin_edges_file is not None:
        edges = read_edges_file(
            args.bin_edges_file
            .expanduser()
            .resolve(strict=True)
        )
        edge_source = str(
            args.bin_edges_file
            .expanduser()
            .resolve(strict=True)
        )
    else:
        edges = estimate_edges(
            transitions,
            dt_values=dt_values,
            n_bins=int(args.n_bins),
            mode=str(args.binning),
            q_low=float(args.q_low),
            q_high=float(args.q_high),
        )
        edge_source = (
            "common4 xmid_latent across selected dt values"
        )

    n_bins = len(edges) - 1

    edges_path = (
        outdir / "01_fluctuation_bin_edges.csv"
    )
    save_edges(
        edges,
        edges_path,
    )

    print("abundance bins :", n_bins)
    print(
        "abundance range:",
        f"{edges[0]:.6f}",
        "to",
        f"{edges[-1]:.6f}",
    )

    # ---------------------------------------------------------------
    # Compact sufficient statistics.
    # ---------------------------------------------------------------
    print("[1/6] Building subject/bin and interval/bin sufficient statistics...")
    subject_bin, interval_bin, subject_dt = (
        build_sufficient_tables(
            transitions,
            dt_values=dt_values,
            edges=edges,
        )
    )

    subject_bin.write_parquet(
        outdir / "08_subject_bin_sufficient.parquet",
        compression="zstd",
        statistics=True,
    )
    subject_dt.write_parquet(
        outdir / "09_subject_dt_sufficient.parquet",
        compression="zstd",
        statistics=True,
    )

    subjects = sorted(
        int(x)
        for x in subject_bin[
            "subject"
        ].unique().to_list()
    )

    if len(subjects) < 2:
        raise RuntimeError(
            "Need at least two biological subjects."
        )

    print("subjects       :", len(subjects))

    # ---------------------------------------------------------------
    # Biological subject-cluster bootstrap.
    # ---------------------------------------------------------------
    print(
        "[2/6] Subject-cluster bootstrap:",
        int(args.n_bootstrap),
        "seed=",
        int(args.bootstrap_seed),
    )

    boot = bootstrap_subject_cells(
        subject_bin,
        subjects=subjects,
        dt_values=dt_values,
        n_bins=n_bins,
        n_bootstrap=int(args.n_bootstrap),
        seed=int(args.bootstrap_seed),
    )

    boot_ci = bootstrap_ci_by_cell(
        boot,
        min_n=int(args.min_n),
    )

    global_boot = bootstrap_global_dt(
        subject_dt,
        subjects=subjects,
        dt_values=dt_values,
        n_bootstrap=int(args.n_bootstrap),
        seed=int(args.bootstrap_seed),
        min_n=int(args.min_n),
    )
    global_boot_ci = bootstrap_global_dt_ci(
        global_boot
    )

    if (
        args.keep_bootstrap_draws
        and not boot.is_empty()
    ):
        boot.write_parquet(
            outdir / "10_bootstrap_draws.parquet",
            compression="zstd",
            statistics=True,
        )

    # ---------------------------------------------------------------
    # Pooled primary curves.
    # ---------------------------------------------------------------
    print("[3/6] Reconstructing pooled dt x abundance-bin metrics...")
    pooled = pooled_by_bin_dt(
        subject_bin,
        interval_bin,
        edges,
        min_n=int(args.min_n),
        bootstrap_ci=boot_ci,
    )
    pooled_path = (
        outdir / "02_fluctuation_by_bin_dt.csv"
    )
    pooled.write_csv(pooled_path)

    # ---------------------------------------------------------------
    # Global dt summaries.
    # ---------------------------------------------------------------
    print("[4/6] Computing global dt summaries...")
    global_dt = global_by_dt(
        subject_dt,
        bootstrap_ci=global_boot_ci,
    )
    global_path = (
        outdir / "03_global_fluctuation_by_dt.csv"
    )
    global_dt.write_csv(global_path)

    # ---------------------------------------------------------------
    # Subject and interval resolved metrics.
    # ---------------------------------------------------------------
    print("[5/6] Computing subject- and interval-resolved metrics...")
    subject_metrics = subject_metrics_by_bin_dt(
        subject_bin,
        edges,
    )
    subject_path = (
        outdir
        / "04_subject_fluctuation_by_bin_dt.csv"
    )
    subject_metrics.write_csv(
        subject_path
    )

    interval_metrics = interval_metrics_by_bin(
        interval_bin,
        edges,
    )
    interval_path = (
        outdir
        / "05_subject_interval_fluctuation_by_bin.csv"
    )
    interval_metrics.write_csv(
        interval_path
    )

    support = support_by_dt(
        transitions,
        dt_values,
    )
    support_path = (
        outdir / "06_support_by_dt.csv"
    )
    support.write_csv(support_path)

    # ---------------------------------------------------------------
    # Identity audit.
    # ---------------------------------------------------------------
    print("[6/6] Covariance identity audit...")
    audit = identity_audit(
        pooled,
        subject_metrics,
        interval_metrics,
    )
    audit_path = (
        outdir / "07_covariance_identity_audit.csv"
    )
    audit.write_csv(audit_path)

    cell_audit = cell_support_audit(
        pooled,
        dt_values=dt_values,
        min_n=int(args.min_n),
    )
    cell_audit_path = (
        outdir / "07b_cell_support_audit.csv"
    )
    cell_audit.write_csv(cell_audit_path)

    # ---------------------------------------------------------------
    # Run config.
    # ---------------------------------------------------------------
    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "input": path_identity(transitions),
        "dataset_label": str(args.dataset_label),
        "dt_values": [int(x) for x in dt_values],
        "primary_support": "common4",
        "primary_conditioning": "xmid_latent",
        "cross_cov_truncated_at_zero": False,
        "ddof": int(DDOF),
        "binning": {
            "edge_source": edge_source,
            "edges": [float(x) for x in edges],
            "requested_mode": str(args.binning),
            "q_low": float(args.q_low),
            "q_high": float(args.q_high),
            "min_n": int(args.min_n),
        },
        "bootstrap": {
            "unit": "biological subject",
            "n_bootstrap": int(args.n_bootstrap),
            "seed": int(args.bootstrap_seed),
            "bootstrap_cell_min_n_matches_point_estimate": True,
        },
    }
    analysis_signature = stable_sha256(signature_payload)

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

    run_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_signature": analysis_signature,
        "dataset_label": str(
            args.dataset_label
        ),
        "transitions": str(
            transitions
        ),
        "outdir": str(
            outdir
        ),
        "dt_values": dt_values,
        "primary_support": "common4",
        "operational_class_filter": None,
        "primary_conditioning": "xmid_latent",
        "primary_estimands": {
            "cross_cov": (
                "Cov(dx_observed_rep1, dx_observed_rep2 | "
                "xmid_latent, dt, common4)"
            ),
            "cross_cov_truncated_at_zero": False,
            "same_var_mean": (
                "0.5*[Var(dx_observed_rep1)+Var(dx_observed_rep2)]"
            ),
            "replicate_specific_excess": (
                "same_var_mean-cross_cov"
            ),
        },
        "ddof": int(DDOF),
        "binning": {
            "source": edge_source,
            "n_bins": int(n_bins),
            "requested_mode": str(
                args.binning
            ),
            "q_low": float(
                args.q_low
            ),
            "q_high": float(
                args.q_high
            ),
            "x_min": float(
                edges[0]
            ),
            "x_max": float(
                edges[-1]
            ),
            "min_n_per_pooled_cell": int(
                args.min_n
            ),
            "edges_file": str(
                edges_path
            ),
        },
        "bootstrap": {
            "unit": "biological subject",
            "n_subjects": int(
                len(subjects)
            ),
            "n_bootstrap": int(
                args.n_bootstrap
            ),
            "seed": int(
                args.bootstrap_seed
            ),
            "cohort_estimand": (
                "transition-weighted pooled covariance reconstructed "
                "within each subject-bootstrap draw"
            ),
            "cell_support_policy": (
                "bootstrap draw contributes to a dt x bin CI only when "
                "reconstructed n >= point-estimate min_n"
            ),
            "global_dt_bootstrap": True,
        },
        "downstream_contract": {
            "step12": (
                "02_fluctuation_by_bin_dt.csv + 00_run_signature.json "
                "for support-matched real-vs-pseudo fluctuation comparison"
            ),
            "step13": (
                "04_subject_fluctuation_by_bin_dt.csv for "
                "subject-resolved temporal scaling"
            ),
            "step14": (
                "05_subject_interval_fluctuation_by_bin.csv for "
                "calendar/common-start/common-end controls"
            ),
        },
    }

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

    # ---------------------------------------------------------------
    # Console summary.
    # ---------------------------------------------------------------
    print("\n[DONE] Step 11 real cross-replicate fluctuation dynamics")
    print("Subjects               :", len(subjects))
    print("Primary support        : common4")
    print("Conditioning           : xmid_latent")
    print("Cross covariance signed: YES")
    print("Outputs:")
    for path in [
        edges_path,
        pooled_path,
        global_path,
        subject_path,
        interval_path,
        support_path,
        audit_path,
        cell_audit_path,
        outdir / "08_subject_bin_sufficient.parquet",
        outdir / "09_subject_dt_sufficient.parquet",
    ]:
        print(" -", path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
