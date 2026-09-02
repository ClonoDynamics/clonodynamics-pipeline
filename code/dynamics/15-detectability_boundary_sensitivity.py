#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15-detectability_boundary_sensitivity.py
========================================

ClonoDynamics Step 15: sensitivity of genuine longitudinal dynamics to the
fixed upstream detectability-state boundary.

Scientific role
---------------
The upstream detectability/observability labels are treated as fixed data.
Step 15 DOES NOT change the Step-2 empirical observability alpha and DOES NOT
reclassify endpoints at alternative thresholds.

Instead, it asks whether the main longitudinal conclusions depend materially on
which already-labelled detectability-state transitions are included.

The sensitivity is estimand-specific.

1. Forward-drift sensitivity
----------------------------
Primary forward drift remains replicate-decoupled and conditioned on the
initial latent state x0.

The comparison is:

    TT
        detectable at both t0 and t1;

    T0PLUS
        detectable at t0 with t1 unrestricted
        = TT + TF.

This preserves the meaning of x0 as a detectable initial-state coordinate while
removing the future-detectability requirement.

The exact Step-9 reciprocal cross-fit geometry is reused:

    AB
        condition on replicate A at t0;
        calculate displacement from replicate B.

    BA
        condition on replicate B at t0;
        calculate displacement from replicate A.

AB and BA are correlated views of the same transitions and are combined with
equal fold weight, never concatenated as independent observations.

T0PLUS is additionally decomposed into its exact disjoint TT and TF
contributions:

    E[dx | T0 detectable]
        = P(TT | T0 detectable) E[dx | TT]
        + P(TF | T0 detectable) E[dx | TF].

TF is derived as the exact T0PLUS-minus-TT remainder within each fold/bin before
the reciprocal-fold combination.

2. Fluctuation-magnitude sensitivity
------------------------------------
The transition-centred Step-11 estimand remains:

    Var(dx_latent | xstar_latent, dt)      [primary]
    MSD(dx_latent | xstar_latent, dt)      [complementary]

The comparison is:

    TT
        detectable at both endpoints;

    NON_FF
        TT + TF + FT, excluding FF.

TT and NON_FF are evaluated on one common absolute xstar support with identical
bin definitions.

3. Temporal-scaling sensitivity
-------------------------------
The current Step-13 analysis is run separately on the Step-11 TT and NON_FF
outputs. Model comparison is summarized only on the largest contiguous
cross-mode common temporal core.

No paired inferential test of model parameters across TT and NON_FF is claimed.
This is a robustness/sensitivity comparison.

4. Boundary-crossing descriptors
--------------------------------
Using the Step-5 transition table, Step 15 reports:

    crossing burden among all transitions
        C_all = (n_TF + n_FT) / n_all

    crossing burden among NON_FF transitions
        C_nonFF = (n_TF + n_FT) / (n_TT + n_TF + n_FT)

    crossing imbalance
        J = (n_FT - n_TF) / (n_TF + n_FT).

J is computed with explicit signed-integer arithmetic and validated to remain in
[-1, 1]. These are descriptive quantities, not transition rates.

Pipeline boundary
-----------------
Step 15 consumes:

    Step 4
        latent trajectories, including replicate-specific endpoint summaries
        needed for Step-9 AB/BA decoupling;

    Step 5
        latent transitions with fixed TT/TF/FT/FF labels, xstar_latent and
        dx_latent;

and reuses current implementations of:

    Step 9
        replicate-decoupled forward drift;

    Step 11
        longitudinal fluctuation dynamics;

    Step 13
        temporal fluctuation scaling.

Step 15 does not alter Step-2 alpha and does not perform an alpha sweep.
A numerical alpha sweep, if ever required, is a separate second-level
sensitivity analysis.

Current public outputs
----------------------
    00_run_config.json
    01_input_manifest.csv
    02_class_composition_by_dt.csv
    03_fluctuation_mode_support_by_dt.csv
    04_fluctuation_common_xstar_support.csv
    05_crossing_metrics_by_dt_bin.csv/parquet
    06_forward_mode_support.csv
    07_forward_bin_edges.csv
    08_forward_TT_vs_T0plus_curves.csv/parquet
    09_forward_T0plus_minus_TT_by_bin.csv/parquet
    10_forward_sensitivity_summary.csv
    10a_forward_pair_build_audit.csv
    10b_forward_TT_TF_decomposition_by_bin.csv/parquet
    10c_forward_decomposition_summary.csv
    11_step11_commands.csv
    12_step13_commands.csv
    13_combined_step11_fluctuation_curves.csv/parquet       [when available]
    14_fluctuation_mode_difference_vs_TT.csv/parquet        [when available]
    15_fluctuation_mode_robustness_summary.csv               [when available]
    16_combined_step13_model_comparison_by_bin.csv/parquet  [when available]
    17_cross_mode_common_temporal_core.csv                   [when available]
    18_temporal_scaling_mode_comparison_core.csv/parquet     [when available]
    19_temporal_scaling_sensitivity_summary.csv              [when available]
    20_step15_summary.csv
    DEBUG_REPORT.txt

Recommended run
---------------
python3 ./code/analysis/15-detectability_boundary_sensitivity.py \
    --trajectories \
    ./dataset_longitudinal_clonodynamics_results/longitudinal/4-latent_trajectory_construction/latent_trajectories_long.parquet \
    --transitions \
    ./dataset_longitudinal_clonodynamics_results/longitudinal/5-latent_transition_construction/latent_transitions.parquet \
    --outdir \
    ./dataset_longitudinal_clonodynamics_results/longitudinal/15-detectability_boundary_sensitivity \
    --step9-script ./code/analysis/9-replicate_decoupled_forward_drift.py \
    --step11-script ./code/analysis/11-longitudinal_fluctuation_dynamics.py \
    --step13-script ./code/analysis/13-temporal_fluctuation_scaling.py \
    --dt-values 1,2,3,4,5 \
    --forward-dt 1 \
    --forward-n-bins 20 \
    --fluctuation-n-bins 30 \
    --min-n 30 \
    --min-subjects 2 \
    --n-bootstrap 2000 \
    --seed 123 \
    --run-step11 \
    --run-step13 \
    --streaming

Python >= 3.9
Polars >= 1.30
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import math
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import polars as pl
except Exception:
    pl = None  # type: ignore


SCRIPT_VERSION = "v3-step15-detectability-boundary-sensitivity-2026-08-25"
FORWARD_MODES = ("TT", "T0PLUS")
FLUCTUATION_MODES = ("TT", "NON_FF")
OBS_CLASSES = ("TT", "TF", "FT", "FF")
CURVE_METRICS = ("mean_dx", "median_dx", "ppos", "msd", "var_dx", "mad_dx")


# =============================================================================
# CLI
# =============================================================================


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 15: test sensitivity of longitudinal dynamics to TT selection "
            "without changing the upstream detectability threshold."
        ),
    )

    parser.add_argument(
        "--trajectories",
        required=True,
        type=Path,
        help="Step-4 genuine-longitudinal latent trajectory table used by Step 9.",
    )
    parser.add_argument(
        "--transitions",
        required=True,
        type=Path,
        help="Step-5 genuine-longitudinal latent transition table used by Step 11.",
    )
    parser.add_argument("--outdir", "--out-dir", dest="outdir", required=True, type=Path)
    parser.add_argument("--dataset-label", default="healthy")
    parser.add_argument(
        "--detectability-alpha",
        type=float,
        default=0.05,
        help=(
            "Upstream threshold represented by existing observable/obs_class labels. "
            "Step 15 records but does not alter this threshold."
        ),
    )

    parser.add_argument(
        "--dt-values",
        default="1,2,3,4,5",
        help="Comma-separated temporal lags for fluctuation and scaling sensitivity.",
    )
    parser.add_argument(
        "--forward-dt",
        type=int,
        default=1,
        help="Temporal lag for the forward-drift boundary sensitivity.",
    )

    parser.add_argument("--forward-n-bins", type=int, default=20)
    parser.add_argument("--forward-q-low", type=float, default=0.005)
    parser.add_argument("--forward-q-high", type=float, default=0.995)
    parser.add_argument(
        "--fluctuation-n-bins",
        type=int,
        default=30,
        help="Common absolute xstar bins used by both TT and NON_FF Step-11 runs.",
    )
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--min-subjects", type=int, default=2)
    parser.add_argument("--n-bootstrap", "--n-boot", dest="n_bootstrap", type=int, default=2000)
    parser.add_argument("--ci", type=float, default=95.0)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--core-min-subject-fraction",
        type=float,
        default=1.0,
        help="Step-13 structural subject-availability fraction.",
    )
    parser.add_argument("--min-dt-points", type=int, default=5)

    parser.add_argument(
        "--step9-script",
        type=Path,
        default=Path("./code/analysis/9-replicate_decoupled_forward_drift.py"),
        help="Current Step-9 analysis script; imported for the exact replicate-decoupled algebra.",
    )
    parser.add_argument(
        "--step11-script",
        type=Path,
        default=Path("./code/analysis/11-longitudinal_fluctuation_dynamics.py"),
    )
    parser.add_argument(
        "--step13-script",
        type=Path,
        default=Path("./code/analysis/13-temporal_fluctuation_scaling.py"),
    )
    parser.add_argument(
        "--python-executable",
        default=sys.executable,
        help="Python executable used for Step-11/13 subprocesses.",
    )

    parser.add_argument(
        "--run-step11",
        action="store_true",
        help="Execute current Step 11 for TT and NON_FF on the common xstar support.",
    )
    parser.add_argument(
        "--run-step13",
        action="store_true",
        help="Execute current Step 13 for TT and NON_FF. Implies --run-step11.",
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help=(
            "Reuse existing Step-11/13 mode directories when their expected outputs exist; "
            "otherwise execute the requested downstream step."
        ),
    )

    parser.add_argument("--streaming", action="store_true")
    parser.add_argument("--compression", default="zstd")
    parser.add_argument(
        "--keep-forward-transitions",
        action="store_true",
        help="Keep temporary Step-15 TT/T0PLUS wide forward transition tables.",
    )
    parser.add_argument(
        "--restart-forward",
        action="store_true",
        help="Delete and rebuild resumable Step-15 forward transition parts.",
    )

    return parser


# =============================================================================
# General utilities
# =============================================================================


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def parse_dt_values(text: str) -> List[int]:
    values = []
    for token in str(text).split(","):
        token = token.strip()
        if token:
            values.append(int(token))
    values = sorted(set(values))
    if not values:
        raise ValueError("--dt-values must contain at least one positive integer.")
    if any(v <= 0 for v in values):
        raise ValueError("Every --dt-values entry must be positive.")
    return values


def validate_args(args: argparse.Namespace) -> None:
    if not args.trajectories.exists():
        raise FileNotFoundError(args.trajectories)
    if not args.transitions.exists():
        raise FileNotFoundError(args.transitions)
    if not args.step9_script.exists():
        raise FileNotFoundError(args.step9_script)
    if args.run_step11 or args.run_step13:
        if not args.step11_script.exists():
            raise FileNotFoundError(args.step11_script)
    if args.run_step13 and not args.step13_script.exists():
        raise FileNotFoundError(args.step13_script)
    if args.run_step13:
        args.run_step11 = True
    if args.forward_dt <= 0:
        raise ValueError("--forward-dt must be positive.")
    if args.forward_n_bins < 3:
        raise ValueError("--forward-n-bins must be >= 3.")
    if args.fluctuation_n_bins < 3:
        raise ValueError("--fluctuation-n-bins must be >= 3.")
    if not (0.0 <= args.forward_q_low < args.forward_q_high <= 1.0):
        raise ValueError("Require 0 <= forward-q-low < forward-q-high <= 1.")
    if args.min_n < 1:
        raise ValueError("--min-n must be >= 1.")
    if args.min_subjects < 1:
        raise ValueError("--min-subjects must be >= 1.")
    if args.n_bootstrap <= 0:
        raise ValueError("--n-bootstrap must be > 0.")
    if not (0.0 < args.ci < 100.0):
        raise ValueError("--ci must be between 0 and 100.")
    if not (0.0 < args.detectability_alpha < 1.0):
        raise ValueError("--detectability-alpha must be in (0,1).")
    if not (0.0 < args.core_min_subject_fraction <= 1.0):
        raise ValueError("--core-min-subject-fraction must be in (0,1].")
    if args.min_dt_points < 2:
        raise ValueError("--min-dt-points must be >= 2.")
    if args.min_dt_points > len(args.dt_values_parsed):
        raise ValueError("--min-dt-points cannot exceed the number of --dt-values.")


def collect_frame(lf: pl.LazyFrame, streaming: bool) -> pl.DataFrame:
    if streaming:
        try:
            return lf.collect(engine="streaming")
        except TypeError:
            return lf.collect(streaming=True)
    return lf.collect()


def scan_table(path: Path) -> pl.LazyFrame:
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pl.scan_parquet(str(path))
    if suffix == ".csv":
        return pl.scan_csv(str(path), infer_schema_length=10000)
    if suffix in {".tsv", ".txt"}:
        return pl.scan_csv(str(path), separator="\t", infer_schema_length=10000)
    if suffix in {".feather", ".arrow", ".ipc"}:
        return pl.scan_ipc(str(path))
    raise ValueError("Unsupported input format: {}".format(path))


def read_table(path: Path) -> pl.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pl.read_parquet(str(path))
    if suffix == ".csv":
        return pl.read_csv(str(path), infer_schema_length=10000)
    if suffix in {".tsv", ".txt"}:
        return pl.read_csv(str(path), separator="\t", infer_schema_length=10000)
    if suffix in {".feather", ".arrow", ".ipc"}:
        return pl.read_ipc(str(path))
    raise ValueError("Unsupported input format: {}".format(path))


def write_csv_parquet(df: pl.DataFrame, stem: Path, compression: str) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(stem.with_suffix(".csv"))
    df.write_parquet(stem.with_suffix(".parquet"), compression=compression)


def find_output_table(folder: Path, stem: str) -> Optional[Path]:
    for suffix in (".parquet", ".pq", ".csv", ".tsv", ".txt"):
        candidate = folder / (stem + suffix)
        if candidate.exists():
            return candidate
    return None


def quote_command(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run_command(command: Sequence[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("COMMAND\n")
        handle.write(quote_command(command) + "\n\nOUTPUT\n")
        handle.flush()
        process = subprocess.run(
            list(command),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if process.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code {}. See {}".format(
                process.returncode, log_path
            )
        )


def schema_names(path: Path) -> List[str]:
    return list(scan_table(path).collect_schema().names())


def largest_contiguous_run(values: Iterable[int]) -> List[int]:
    ids = sorted(set(int(v) for v in values))
    if not ids:
        return []
    runs: List[List[int]] = []
    current = [ids[0]]
    for value in ids[1:]:
        if value == current[-1] + 1:
            current.append(value)
        else:
            runs.append(current)
            current = [value]
    runs.append(current)
    runs.sort(key=lambda run: (len(run), -run[0]), reverse=True)
    return runs[0]


def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    aa = a[mask]
    bb = b[mask]
    if aa.size < 2 or np.std(aa) <= 0 or np.std(bb) <= 0:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    xx = x[mask]
    yy = y[mask]
    if xx.size < 2 or np.std(xx) <= 0:
        return float("nan")
    return float(np.polyfit(xx, yy, 1)[0])


# =============================================================================
# Step-5 detectability composition and fluctuation support
# =============================================================================


def prepare_transition_base(
    path: Path,
    dt_values: Sequence[int],
    streaming: bool,
) -> pl.DataFrame:
    lf = scan_table(path)
    names = list(lf.collect_schema().names())
    required = ["subject", "aaSeqCDR3", "t0", "t1", "obs_class", "xstar_latent"]
    missing = [c for c in required if c not in names]
    if missing:
        raise ValueError("Step-5 transition table missing required columns: {}".format(missing))

    dt_expr = (
        pl.col("dt").cast(pl.Int32, strict=False)
        if "dt" in names
        else (
            pl.col("t1").cast(pl.Int32, strict=False)
            - pl.col("t0").cast(pl.Int32, strict=False)
        )
    )
    if "dx_latent" in names:
        dx_expr = pl.col("dx_latent").cast(pl.Float64, strict=False)
    else:
        missing_dx = [c for c in ("x0_latent", "x1_latent") if c not in names]
        if missing_dx:
            raise ValueError(
                "Step-5 transition table requires dx_latent or both x0_latent/x1_latent; "
                "missing {}".format(missing_dx)
            )
        dx_expr = (
            pl.col("x1_latent").cast(pl.Float64, strict=False)
            - pl.col("x0_latent").cast(pl.Float64, strict=False)
        )

    base = (
        lf.select(
            [
                pl.col("subject").cast(pl.String, strict=False).alias("subject"),
                pl.col("aaSeqCDR3").cast(pl.String, strict=False).alias("aaSeqCDR3"),
                pl.col("t0").cast(pl.Int32, strict=False).alias("t0"),
                pl.col("t1").cast(pl.Int32, strict=False).alias("t1"),
                dt_expr.alias("dt"),
                pl.col("obs_class")
                .cast(pl.String, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
                .alias("obs_class"),
                pl.col("xstar_latent").cast(pl.Float64, strict=False).alias("xstar_latent"),
                dx_expr.alias("dx_latent"),
            ]
        )
        .filter(
            pl.col("subject").is_not_null()
            & pl.col("aaSeqCDR3").is_not_null()
            & pl.col("dt").is_in([int(v) for v in dt_values])
            & pl.col("obs_class").is_in(list(OBS_CLASSES))
            & pl.col("xstar_latent").is_finite()
            & pl.col("dx_latent").is_finite()
        )
    )
    data = collect_frame(base, streaming)
    if data.height == 0:
        raise ValueError("No Step-5 transitions remain after Step-15 filtering.")

    represented = set(int(v) for v in data["dt"].unique().to_list())
    missing_dt = sorted(set(int(v) for v in dt_values) - represented)
    if missing_dt:
        raise ValueError("Requested dt values absent from Step-5 transitions: {}".format(missing_dt))
    return data


def class_composition_by_dt(data: pl.DataFrame) -> pl.DataFrame:
    counts = (
        data.group_by(["dt", "obs_class"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.struct(["subject", "aaSeqCDR3"]).n_unique().alias("n_subject_clonotype_pairs"),
            ]
        )
        .sort(["dt", "obs_class"])
    )
    totals = counts.group_by("dt").agg(pl.col("n_transitions").sum().alias("n_total"))
    return (
        counts.join(totals, on="dt", how="left")
        .with_columns((pl.col("n_transitions") / pl.col("n_total")).alias("fraction"))
        .sort(["dt", "obs_class"])
    )


def fluctuation_mode_support(data: pl.DataFrame, dt_values: Sequence[int]) -> pl.DataFrame:
    rows: List[pl.DataFrame] = []
    for mode in FLUCTUATION_MODES:
        if mode == "TT":
            sub = data.filter(pl.col("obs_class") == "TT")
        else:
            sub = data.filter(pl.col("obs_class") != "FF")
        by_dt = (
            sub.group_by("dt")
            .agg(
                [
                    pl.col("xstar_latent").min().alias("xstar_min"),
                    pl.col("xstar_latent").max().alias("xstar_max"),
                    pl.len().alias("n_transitions"),
                    pl.col("subject").n_unique().alias("n_subjects"),
                ]
            )
            .sort("dt")
            .with_columns(pl.lit(mode).alias("transition_mode"))
        )
        represented = set(int(v) for v in by_dt["dt"].to_list())
        missing = sorted(set(int(v) for v in dt_values) - represented)
        if missing:
            raise ValueError("Mode {} lacks dt values {}".format(mode, missing))
        rows.append(by_dt)
    return pl.concat(rows, how="diagonal_relaxed").select(
        ["transition_mode", "dt", "xstar_min", "xstar_max", "n_transitions", "n_subjects"]
    ).sort(["transition_mode", "dt"])


def resolve_fluctuation_common_support(support: pl.DataFrame) -> Tuple[float, float]:
    x_min = float(support["xstar_min"].max())
    x_max = float(support["xstar_max"].min())
    if not math.isfinite(x_min) or not math.isfinite(x_max) or x_max <= x_min:
        raise ValueError(
            "No common absolute xstar support across TT/NON_FF x dt: [{}, {}]".format(
                x_min, x_max
            )
        )
    return x_min, x_max


def common_xstar_edges(x_min: float, x_max: float, n_bins: int) -> np.ndarray:
    return np.linspace(float(x_min), float(x_max), int(n_bins) + 1, dtype=float)


def bin_expr(column: str, edges: np.ndarray) -> pl.Expr:
    # Equal-width edges are used for Step-15 common-support descriptors.
    x_min = float(edges[0])
    x_max = float(edges[-1])
    n_bins = len(edges) - 1
    width = (x_max - x_min) / float(n_bins)
    raw = ((pl.col(column) - x_min) / width).floor().cast(pl.Int32)
    return (
        pl.when((pl.col(column) >= x_min) & (pl.col(column) <= x_max))
        .then(
            pl.when(pl.col(column) == x_max)
            .then(pl.lit(n_bins - 1, dtype=pl.Int32))
            .otherwise(raw)
        )
        .otherwise(pl.lit(None, dtype=pl.Int32))
        .alias("bin_id")
    )


def build_crossing_metrics(data: pl.DataFrame, edges: np.ndarray) -> pl.DataFrame:
    edge_df = pl.DataFrame(
        {
            "bin_id": np.arange(len(edges) - 1, dtype=np.int32),
            "bin_left": edges[:-1],
            "bin_right": edges[1:],
            "bin_mid": 0.5 * (edges[:-1] + edges[1:]),
        }
    )
    binned = (
        data.with_columns(bin_expr("xstar_latent", edges))
        .filter(pl.col("bin_id").is_not_null())
        .group_by(["dt", "bin_id"])
        .agg(
            [
                pl.len().alias("n_all"),
                pl.col("subject").n_unique().alias("n_subjects"),
                (pl.col("obs_class") == "TT").sum().alias("n_TT"),
                (pl.col("obs_class") == "TF").sum().alias("n_TF"),
                (pl.col("obs_class") == "FT").sum().alias("n_FT"),
                (pl.col("obs_class") == "FF").sum().alias("n_FF"),
            ]
        )
        .with_columns(
            [
                (pl.col("n_TT") + pl.col("n_TF") + pl.col("n_FT")).alias("n_NON_FF"),
                (pl.col("n_TF") + pl.col("n_FT")).alias("n_crossing"),
            ]
        )
        .with_columns(
            [
                (pl.col("n_crossing") / pl.col("n_all")).alias("crossing_burden_all"),
                pl.when(pl.col("n_NON_FF") > 0)
                .then(pl.col("n_crossing") / pl.col("n_NON_FF"))
                .otherwise(None)
                .alias("crossing_burden_non_ff"),
                pl.when(pl.col("n_crossing") > 0)
                .then(
                    (
                        pl.col("n_FT").cast(pl.Int64)
                        - pl.col("n_TF").cast(pl.Int64)
                    )
                    / pl.col("n_crossing").cast(pl.Float64)
                )
                .otherwise(None)
                .alias("crossing_imbalance"),
                (pl.col("n_TT") / pl.col("n_all")).alias("tt_fraction_all"),
            ]
        )
    )
    dt_grid = data.select(pl.col("dt").unique().sort()).join(edge_df, how="cross")
    out = dt_grid.join(binned, on=["dt", "bin_id"], how="left")
    count_cols = [
        "n_all", "n_subjects", "n_TT", "n_TF", "n_FT", "n_FF",
        "n_NON_FF", "n_crossing",
    ]
    existing_counts = [c for c in count_cols if c in out.columns]
    if existing_counts:
        out = out.with_columns([pl.col(c).fill_null(0) for c in existing_counts])

    # J=(FT-TF)/(FT+TF) is mathematically bounded in [-1, 1].
    # Explicit signed casting above prevents unsigned-integer underflow when TF > FT.
    bad_j = out.filter(
        pl.col("crossing_imbalance").is_not_null()
        & (
            (pl.col("crossing_imbalance") < -1.0 - 1e-12)
            | (pl.col("crossing_imbalance") > 1.0 + 1e-12)
        )
    )
    if bad_j.height:
        raise RuntimeError(
            "Internal error: crossing_imbalance outside [-1,1] after signed casting."
        )
    return out.sort(["dt", "bin_id"])


# =============================================================================
# Step-9 import and forward TT vs T0PLUS sensitivity
# =============================================================================


def import_step9_module(path: Path):
    spec = importlib.util.spec_from_file_location("clonodynamics_step9_for_step15", str(path))
    if spec is None or spec.loader is None:
        raise ImportError("Could not import Step-9 script: {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    required_functions = [
        "schema_names",
        "discover_subject_times",
        "load_subject",
        "build_pair",
        "merge_parts",
        "point_stats_one",
        "subject_sufficient_stats_one",
        "combine_cross_point",
    ]
    missing = [name for name in required_functions if not hasattr(module, name)]
    if missing:
        raise AttributeError("Step-9 script lacks required functions: {}".format(missing))
    return module


def t0plus_endpoint_frame(subject_df: pl.DataFrame, time_value: int, suffix: str) -> pl.DataFrame:
    """Step-9-compatible endpoint frame: require detectability only at t0."""
    filt = pl.col("time") == int(time_value)
    if str(suffix) == "t0":
        filt = filt & pl.col("observable")
    endpoint = subject_df.filter(filt).drop(["subject", "time", "observable"])
    rename_map = {
        c: "{}_{}".format(c, suffix)
        for c in endpoint.columns
        if c != "aaSeqCDR3"
    }
    return endpoint.rename(rename_map)



def build_forward_mode_parts(
    step9,
    trajectories: Path,
    mode: str,
    parts_dir: Path,
    dt: int,
    compression: str,
    restart: bool,
) -> Tuple[List[Path], List[Dict[str, Any]]]:
    mode = str(mode).upper()
    if mode not in FORWARD_MODES:
        raise ValueError("Unsupported forward mode: {}".format(mode))

    available = step9.schema_names(trajectories)
    missing = [c for c in step9.REQUIRED_COLUMNS if c not in available]
    if missing:
        raise ValueError(
            "Step-4 trajectory table missing Step-9 replicate-decoupling columns: {}".format(
                missing
            )
        )

    if restart and parts_dir.exists():
        shutil.rmtree(parts_dir)
    ensure_dir(parts_dir)

    subject_times = step9.discover_subject_times(trajectories)
    original_endpoint = step9.endpoint_frame
    if mode == "T0PLUS":
        step9.endpoint_frame = t0plus_endpoint_frame

    parts: List[Path] = []
    audit: List[Dict[str, Any]] = []
    try:
        for subject in sorted(subject_times):
            times = sorted(subject_times[subject])
            time_set = set(times)
            valid_pairs = [(t, t + int(dt)) for t in times if (t + int(dt)) in time_set]
            if not valid_pairs:
                continue
            sdf = step9.load_subject(trajectories, subject, available)
            for t0, t1 in valid_pairs:
                part = parts_dir / "subject_{}_t{}_t{}.parquet".format(subject, t0, t1)
                if part.exists() and not restart:
                    parts.append(part)
                    continue
                pair = step9.build_pair(sdf, subject, t0, t1)
                n = int(pair.height)
                audit.append(
                    {
                        "transition_mode": mode,
                        "subject": int(subject),
                        "t0": int(t0),
                        "t1": int(t1),
                        "dt": int(dt),
                        "n_forward_usable": n,
                    }
                )
                if n == 0:
                    continue
                pair.write_parquet(part, compression=compression, statistics=True)
                parts.append(part)
                del pair
                gc.collect()
            del sdf
            gc.collect()
    finally:
        step9.endpoint_frame = original_endpoint

    if not parts:
        raise RuntimeError("No forward transition parts were generated for mode {}".format(mode))
    return parts, audit


def forward_mode_support(
    mode_paths: Dict[str, Path],
    q_low: float,
    q_high: float,
    streaming: bool,
) -> Tuple[pl.DataFrame, float, float]:
    rows: List[Dict[str, Any]] = []
    for mode, path in mode_paths.items():
        values = collect_frame(
            pl.scan_parquet(str(path))
            .select(pl.col("x0_joint").cast(pl.Float64, strict=False).alias("x"))
            .filter(pl.col("x").is_finite()),
            streaming,
        )["x"].to_numpy().astype(float)
        if values.size == 0:
            raise ValueError("No finite x0_joint values for forward mode {}".format(mode))
        lo = float(np.quantile(values, q_low))
        hi = float(np.quantile(values, q_high))
        rows.append(
            {
                "transition_mode": mode,
                "n_transitions": int(values.size),
                "x0_min": float(np.min(values)),
                "x0_q_low": lo,
                "x0_median": float(np.median(values)),
                "x0_q_high": hi,
                "x0_max": float(np.max(values)),
                "q_low": float(q_low),
                "q_high": float(q_high),
            }
        )
    table = pl.DataFrame(rows).sort("transition_mode")
    common_min = float(table["x0_q_low"].max())
    common_max = float(table["x0_q_high"].min())
    if not math.isfinite(common_min) or not math.isfinite(common_max) or common_max <= common_min:
        raise ValueError("No valid common forward x0 support.")
    return table, common_min, common_max


def forward_reference_edges(
    tt_path: Path,
    common_min: float,
    common_max: float,
    n_bins: int,
    streaming: bool,
) -> np.ndarray:
    vals = collect_frame(
        pl.scan_parquet(str(tt_path))
        .select(pl.col("x0_joint").cast(pl.Float64, strict=False).alias("x"))
        .filter(
            pl.col("x").is_finite()
            & (pl.col("x") >= float(common_min))
            & (pl.col("x") <= float(common_max))
        ),
        streaming,
    )["x"].to_numpy().astype(float)
    if vals.size < n_bins:
        raise ValueError("Too few TT forward transitions to define {} bins.".format(n_bins))
    probs = np.linspace(0.0, 1.0, int(n_bins) + 1)
    edges = np.quantile(vals, probs)
    edges[0] = float(common_min)
    edges[-1] = float(common_max)
    edges = np.unique(edges.astype(float))
    if edges.size < 4:
        raise ValueError("Forward TT reference quantile edges collapsed to < 3 bins.")
    return edges


def paired_forward_bootstrap(
    suff: Dict[Tuple[str, str], pl.DataFrame],
    n_bins: int,
    n_bootstrap: int,
    seed: int,
    ci: float,
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, np.ndarray]]:
    subjects = sorted(
        {
            str(v)
            for df in suff.values()
            for v in df["subject"].unique().to_list()
        }
    )
    if len(subjects) < 2:
        raise RuntimeError("Need at least two subjects for paired forward bootstrap.")
    subj_to_i = {s: i for i, s in enumerate(subjects)}

    arrays: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for key, df in suff.items():
        n_arr = np.zeros((len(subjects), n_bins), dtype=float)
        sum_arr = np.zeros_like(n_arr)
        pos_arr = np.zeros_like(n_arr)
        for row in df.iter_rows(named=True):
            i = subj_to_i[str(row["subject"])]
            b = int(row["bin"])
            if 0 <= b < n_bins:
                n_arr[i, b] = float(row["n"])
                sum_arr[i, b] = float(row["sum_dx"])
                pos_arr[i, b] = float(row["n_pos"])
        arrays[key] = (n_arr, sum_arr, pos_arr)

    combined_draws: Dict[str, np.ndarray] = {
        mode: np.full((n_bootstrap, n_bins), np.nan, dtype=float) for mode in FORWARD_MODES
    }
    ppos_draws: Dict[str, np.ndarray] = {
        mode: np.full((n_bootstrap, n_bins), np.nan, dtype=float) for mode in FORWARD_MODES
    }

    rng = np.random.default_rng(int(seed))
    for boot in range(int(n_bootstrap)):
        draw = rng.integers(0, len(subjects), size=len(subjects))
        multiplicity = np.bincount(draw, minlength=len(subjects)).astype(float)[:, None]
        for mode in FORWARD_MODES:
            fold_means = []
            fold_ppos = []
            for fold in ("AB", "BA"):
                n_arr, sum_arr, pos_arr = arrays[(mode, fold)]
                n = np.sum(multiplicity * n_arr, axis=0)
                s = np.sum(multiplicity * sum_arr, axis=0)
                p = np.sum(multiplicity * pos_arr, axis=0)
                mean = np.divide(s, n, out=np.full(n_bins, np.nan), where=n > 0)
                ppos = np.divide(p, n, out=np.full(n_bins, np.nan), where=n > 0)
                fold_means.append(mean)
                fold_ppos.append(ppos)
            combined_draws[mode][boot, :] = 0.5 * (fold_means[0] + fold_means[1])
            ppos_draws[mode][boot, :] = 0.5 * (fold_ppos[0] + fold_ppos[1])

    alpha = (100.0 - float(ci)) / 200.0
    summaries: Dict[str, Dict[str, np.ndarray]] = {}
    for mode in FORWARD_MODES:
        summaries[mode] = {
            "mean_ci_low": np.nanquantile(combined_draws[mode], alpha, axis=0),
            "mean_ci_high": np.nanquantile(combined_draws[mode], 1.0 - alpha, axis=0),
            "ppos_ci_low": np.nanquantile(ppos_draws[mode], alpha, axis=0),
            "ppos_ci_high": np.nanquantile(ppos_draws[mode], 1.0 - alpha, axis=0),
        }

    delta = combined_draws["T0PLUS"] - combined_draws["TT"]
    delta_summary = {
        "delta_ci_low": np.nanquantile(delta, alpha, axis=0),
        "delta_ci_high": np.nanquantile(delta, 1.0 - alpha, axis=0),
        "delta_median": np.nanmedian(delta, axis=0),
    }
    return summaries, delta_summary



def build_forward_tt_tf_decomposition(
    point_folds: Dict[str, Dict[str, pl.DataFrame]],
    edges: np.ndarray,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Decompose T0PLUS=TT+TF algebraically within each replicate-decoupled fold.

    TT is a strict subset of T0PLUS after the same fold-specific finite-value and
    x0-bin filters. Therefore the TF component is identified exactly as the
    disjoint remainder:

        n_TF      = n_T0PLUS - n_TT
        sum_dx_TF = n_T0PLUS*mean_T0PLUS - n_TT*mean_TT

    This avoids constructing a second TF transition table and guarantees that the
    decomposition targets exactly the same analyzed T0PLUS rows. AB and BA are
    then combined with equal fold weight; their counts are never pooled as an
    independent sample size.
    """
    for mode in ("TT", "T0PLUS"):
        if mode not in point_folds:
            raise ValueError("Missing forward decomposition mode: {}".format(mode))
        for fold in ("AB", "BA"):
            if fold not in point_folds[mode]:
                raise ValueError(
                    "Missing forward decomposition fold {} for mode {}".format(fold, mode)
                )

    n_bins = len(edges) - 1

    def row_lookup(df: pl.DataFrame) -> Dict[int, Dict[str, Any]]:
        return {int(row["bin"]): row for row in df.iter_rows(named=True)}

    fold_rows: List[Dict[str, Any]] = []
    for fold in ("AB", "BA"):
        tt_lookup = row_lookup(point_folds["TT"][fold])
        tp_lookup = row_lookup(point_folds["T0PLUS"][fold])

        for b in range(n_bins):
            r_tt = tt_lookup.get(b, {})
            r_tp = tp_lookup.get(b, {})

            n_tt = int(r_tt.get("n", 0) or 0)
            n_tp = int(r_tp.get("n", 0) or 0)
            if n_tt > n_tp:
                raise RuntimeError(
                    "Forward decomposition failed: n_TT > n_T0PLUS in fold {} bin {}.".format(
                        fold, b
                    )
                )
            n_tf = int(n_tp - n_tt)

            mean_tt = float(r_tt.get("mean_dx", np.nan)) if n_tt > 0 else np.nan
            mean_tp = float(r_tp.get("mean_dx", np.nan)) if n_tp > 0 else np.nan

            sum_tt = float(n_tt * mean_tt) if n_tt > 0 and np.isfinite(mean_tt) else 0.0
            sum_tp = float(n_tp * mean_tp) if n_tp > 0 and np.isfinite(mean_tp) else np.nan
            if np.isfinite(sum_tp):
                sum_tf = float(sum_tp - sum_tt)
            else:
                sum_tf = np.nan
            mean_tf = float(sum_tf / n_tf) if n_tf > 0 and np.isfinite(sum_tf) else np.nan

            p_tt = float(n_tt / n_tp) if n_tp > 0 else np.nan
            p_tf = float(n_tf / n_tp) if n_tp > 0 else np.nan

            contrib_tt = p_tt * mean_tt if np.isfinite(p_tt) and np.isfinite(mean_tt) else np.nan
            if n_tf == 0 and n_tp > 0:
                contrib_tf = 0.0
            else:
                contrib_tf = p_tf * mean_tf if np.isfinite(p_tf) and np.isfinite(mean_tf) else np.nan

            reconstructed = (
                contrib_tt + contrib_tf
                if np.isfinite(contrib_tt) and np.isfinite(contrib_tf)
                else np.nan
            )
            residual = (
                mean_tp - reconstructed
                if np.isfinite(mean_tp) and np.isfinite(reconstructed)
                else np.nan
            )

            fold_rows.append(
                {
                    "bin": int(b),
                    "x_left": float(edges[b]),
                    "x_right": float(edges[b + 1]),
                    "x_center": float(0.5 * (edges[b] + edges[b + 1])),
                    "fold": fold,
                    "n_TT": n_tt,
                    "n_TF_derived": n_tf,
                    "n_T0PLUS": n_tp,
                    "p_TT_given_T0_detectable": p_tt,
                    "p_TF_given_T0_detectable": p_tf,
                    "mean_dx_TT": mean_tt,
                    "mean_dx_TF_derived": mean_tf,
                    "mean_dx_T0PLUS_direct": mean_tp,
                    "TT_contribution_to_T0PLUS_mean": contrib_tt,
                    "TF_contribution_to_T0PLUS_mean": contrib_tf,
                    "mean_dx_T0PLUS_reconstructed": reconstructed,
                    "reconstruction_residual": residual,
                }
            )

    fold_df = pl.DataFrame(fold_rows).sort(["fold", "bin"])

    def finite_mean(values: Sequence[Any]) -> float:
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        return float(np.mean(arr)) if arr.size else float("nan")

    combined_rows: List[Dict[str, Any]] = []
    for b in range(n_bins):
        rows = list(
            fold_df.filter(pl.col("bin") == int(b)).iter_rows(named=True)
        )
        combined_rows.append(
            {
                "bin": int(b),
                "x_left": float(edges[b]),
                "x_right": float(edges[b + 1]),
                "x_center": float(0.5 * (edges[b] + edges[b + 1])),
                "fold": "combined",
                "n_TT": None,
                "n_TF_derived": None,
                "n_T0PLUS": None,
                "p_TT_given_T0_detectable": finite_mean(
                    [r["p_TT_given_T0_detectable"] for r in rows]
                ),
                "p_TF_given_T0_detectable": finite_mean(
                    [r["p_TF_given_T0_detectable"] for r in rows]
                ),
                "mean_dx_TT": finite_mean([r["mean_dx_TT"] for r in rows]),
                "mean_dx_TF_derived": finite_mean(
                    [r["mean_dx_TF_derived"] for r in rows]
                ),
                "mean_dx_T0PLUS_direct": finite_mean(
                    [r["mean_dx_T0PLUS_direct"] for r in rows]
                ),
                "TT_contribution_to_T0PLUS_mean": finite_mean(
                    [r["TT_contribution_to_T0PLUS_mean"] for r in rows]
                ),
                "TF_contribution_to_T0PLUS_mean": finite_mean(
                    [r["TF_contribution_to_T0PLUS_mean"] for r in rows]
                ),
                "mean_dx_T0PLUS_reconstructed": finite_mean(
                    [r["mean_dx_T0PLUS_reconstructed"] for r in rows]
                ),
                "reconstruction_residual": finite_mean(
                    [r["reconstruction_residual"] for r in rows]
                ),
            }
        )

    combined_df = pl.DataFrame(combined_rows)
    decomposition = pl.concat([fold_df, combined_df], how="diagonal_relaxed").sort(
        ["bin", "fold"]
    )

    c = combined_df.sort("bin")
    residual = c["reconstruction_residual"].to_numpy().astype(float)
    p_tf = c["p_TF_given_T0_detectable"].to_numpy().astype(float)
    mean_tf = c["mean_dx_TF_derived"].to_numpy().astype(float)
    tf_contrib = c["TF_contribution_to_T0PLUS_mean"].to_numpy().astype(float)
    tt_contrib = c["TT_contribution_to_T0PLUS_mean"].to_numpy().astype(float)
    direct = c["mean_dx_T0PLUS_direct"].to_numpy().astype(float)
    recon = c["mean_dx_T0PLUS_reconstructed"].to_numpy().astype(float)
    x = c["x_center"].to_numpy().astype(float)

    finite_res = np.isfinite(residual)
    summary = pl.DataFrame(
        [
            {
                "decomposition": "T0PLUS = TT + TF (TF derived as exact remainder)",
                "n_bins": int(n_bins),
                "mean_p_TF_given_T0_detectable": float(np.nanmean(p_tf)),
                "median_p_TF_given_T0_detectable": float(np.nanmedian(p_tf)),
                "mean_dx_TF_across_bins": float(np.nanmean(mean_tf)),
                "tf_signed_slope": linear_slope(x, mean_tf),
                "mean_abs_TT_contribution": float(np.nanmean(np.abs(tt_contrib))),
                "mean_abs_TF_contribution": float(np.nanmean(np.abs(tf_contrib))),
                "mean_abs_reconstruction_error": (
                    float(np.nanmean(np.abs(residual))) if np.any(finite_res) else float("nan")
                ),
                "max_abs_reconstruction_error": (
                    float(np.nanmax(np.abs(residual))) if np.any(finite_res) else float("nan")
                ),
                "pearson_direct_vs_reconstructed": safe_pearson(direct, recon),
                "fold_combination": "equal AB/BA weight; folds not pooled as independent N",
            }
        ]
    )
    return decomposition, summary


def build_forward_analysis(
    step9,
    trajectories: Path,
    outdir: Path,
    dt: int,
    n_bins: int,
    q_low: float,
    q_high: float,
    n_bootstrap: int,
    seed: int,
    ci: float,
    compression: str,
    streaming: bool,
    restart: bool,
    keep_transitions: bool,
) -> Tuple[Any, ...]:
    temp_root = ensure_dir(outdir / "_forward_boundary_work")
    mode_paths: Dict[str, Path] = {}
    audit_rows: List[Dict[str, Any]] = []

    for mode in FORWARD_MODES:
        parts_dir = temp_root / "parts_{}".format(mode.lower())
        parts, audit = build_forward_mode_parts(
            step9=step9,
            trajectories=trajectories,
            mode=mode,
            parts_dir=parts_dir,
            dt=dt,
            compression=compression,
            restart=restart,
        )
        audit_rows.extend(audit)
        merged = temp_root / "forward_{}_dt{}.parquet".format(mode.lower(), dt)
        if merged.exists() and restart:
            merged.unlink()
        if not merged.exists():
            step9.merge_parts(parts, merged, compression, 100000)
        mode_paths[mode] = merged

    support, common_min, common_max = forward_mode_support(
        {mode: mode_paths[mode] for mode in FORWARD_MODES}, q_low, q_high, streaming
    )
    support = support.with_columns(
        [
            pl.lit(common_min).alias("common_x0_min"),
            pl.lit(common_max).alias("common_x0_max"),
            pl.lit("intersection_of_mode_specific_x0_quantile_ranges").alias("support_rule"),
        ]
    )

    edges = forward_reference_edges(
        mode_paths["TT"], common_min, common_max, n_bins, streaming
    )
    edge_table = pl.DataFrame(
        {
            "edge_index": np.arange(edges.size, dtype=np.int32),
            "edge": edges,
            "reference_mode": ["TT"] * edges.size,
            "definition": ["TT quantile edges within common TT/T0PLUS x0 support"] * edges.size,
        }
    )

    point_frames: List[pl.DataFrame] = []
    suff: Dict[Tuple[str, str], pl.DataFrame] = {}
    mode_points: Dict[str, pl.DataFrame] = {}
    point_folds: Dict[str, Dict[str, pl.DataFrame]] = {}

    for mode in FORWARD_MODES:
        point_fold: Dict[str, pl.DataFrame] = {}
        for fold, rep in (("AB", "cross_AB"), ("BA", "cross_BA")):
            point = step9.point_stats_one(mode_paths[mode], edges, rep)
            point_fold[fold] = point
            suff[(mode, fold)] = step9.subject_sufficient_stats_one(
                mode_paths[mode], edges, rep
            )
        point_folds[mode] = point_fold

        combined = step9.combine_cross_point(point_fold["AB"], point_fold["BA"], edges)
        combined = combined.with_columns(
            [
                pl.lit(mode).alias("transition_mode"),
                pl.lit(
                    "T at t0 and t1"
                    if mode == "TT"
                    else "T at t0; t1 unrestricted"
                ).alias("mode_definition"),
            ]
        )
        mode_points[mode] = combined
        point_frames.append(combined)

    n_actual_bins = len(edges) - 1
    boot_summary, delta_boot = paired_forward_bootstrap(
        suff=suff,
        n_bins=n_actual_bins,
        n_bootstrap=n_bootstrap,
        seed=seed,
        ci=ci,
    )

    curve_frames: List[pl.DataFrame] = []
    for mode in FORWARD_MODES:
        base = mode_points[mode]
        ci_df = pl.DataFrame(
            {
                "bin": np.arange(n_actual_bins, dtype=np.int32),
                "mean_dx_ci_low": boot_summary[mode]["mean_ci_low"],
                "mean_dx_ci_high": boot_summary[mode]["mean_ci_high"],
                "p_dx_gt0_ci_low": boot_summary[mode]["ppos_ci_low"],
                "p_dx_gt0_ci_high": boot_summary[mode]["ppos_ci_high"],
            }
        )
        curve_frames.append(base.join(ci_df, on="bin", how="left"))
    curves = pl.concat(curve_frames, how="diagonal_relaxed").sort(["transition_mode", "bin"])

    tt = curves.filter(pl.col("transition_mode") == "TT").select(
        [
            "bin",
            "x_left",
            "x_right",
            "x_center",
            pl.col("mean_dx").alias("mean_dx_TT"),
            pl.col("p_dx_gt0").alias("p_dx_gt0_TT"),
        ]
    )
    inc = curves.filter(pl.col("transition_mode") == "T0PLUS").select(
        [
            "bin",
            pl.col("mean_dx").alias("mean_dx_T0PLUS"),
            pl.col("p_dx_gt0").alias("p_dx_gt0_T0PLUS"),
        ]
    )
    diff = (
        tt.join(inc, on="bin", how="inner")
        .with_columns(
            [
                (pl.col("mean_dx_T0PLUS") - pl.col("mean_dx_TT")).alias("delta_mean_dx_T0PLUS_minus_TT"),
                (pl.col("p_dx_gt0_T0PLUS") - pl.col("p_dx_gt0_TT")).alias("delta_ppos_T0PLUS_minus_TT"),
            ]
        )
        .join(
            pl.DataFrame(
                {
                    "bin": np.arange(n_actual_bins, dtype=np.int32),
                    "delta_mean_dx_boot_median": delta_boot["delta_median"],
                    "delta_mean_dx_ci_low": delta_boot["delta_ci_low"],
                    "delta_mean_dx_ci_high": delta_boot["delta_ci_high"],
                }
            ),
            on="bin",
            how="left",
        )
        .sort("bin")
    )

    x = diff["x_center"].to_numpy().astype(float)
    y_tt = diff["mean_dx_TT"].to_numpy().astype(float)
    y_inc = diff["mean_dx_T0PLUS"].to_numpy().astype(float)
    d = y_inc - y_tt
    summary = pl.DataFrame(
        [
            {
                "comparison": "T0PLUS_vs_TT",
                "n_matched_bins": int(np.sum(np.isfinite(y_tt) & np.isfinite(y_inc))),
                "pearson_r_curve": safe_pearson(y_tt, y_inc),
                "mean_abs_difference": float(np.nanmean(np.abs(d))),
                "rms_difference": float(np.sqrt(np.nanmean(d ** 2))),
                "max_abs_difference": float(np.nanmax(np.abs(d))),
                "tt_signed_slope": linear_slope(x, y_tt),
                "t0plus_signed_slope": linear_slope(x, y_inc),
                "mean_abs_drift_TT": float(np.nanmean(np.abs(y_tt))),
                "mean_abs_drift_T0PLUS": float(np.nanmean(np.abs(y_inc))),
                "common_x0_min": float(common_min),
                "common_x0_max": float(common_max),
                "n_bins": int(n_actual_bins),
                "bootstrap_replicates": int(n_bootstrap),
            }
        ]
    )

    decomposition, decomposition_summary = build_forward_tt_tf_decomposition(
        point_folds=point_folds,
        edges=edges,
    )

    if keep_transitions:
        keep_dir = ensure_dir(outdir / "forward_transition_tables")
        for mode, source in mode_paths.items():
            target = keep_dir / source.name
            if target.exists():
                target.unlink()
            shutil.copy2(source, target)

    if not keep_transitions:
        shutil.rmtree(temp_root, ignore_errors=True)

    return (
        support, edge_table, curves, diff, summary, audit_rows,
        decomposition, decomposition_summary,
    )


# =============================================================================
# Step-11 and Step-13 command construction / combination
# =============================================================================


def build_step11_command(
    args: argparse.Namespace,
    mode: str,
    output_dir: Path,
    x_min: float,
    x_max: float,
    dt_values: Sequence[int],
) -> List[str]:
    command = [
        str(args.python_executable),
        str(args.step11_script),
        "--transitions",
        str(args.transitions),
        "--dataset-label",
        "{}_step15_{}".format(args.dataset_label, mode.lower()),
        "--outdir",
        str(output_dir),
        "--representation",
        "latent",
        "--conditioning",
        "xstar",
        "--mode",
        mode,
        "--dt",
    ] + [str(v) for v in dt_values] + [
        "--shared-range",
        "intersection",
        "--x-min",
        "{:.17g}".format(float(x_min)),
        "--x-max",
        "{:.17g}".format(float(x_max)),
        "--n-bins",
        str(int(args.fluctuation_n_bins)),
        "--min-n",
        str(int(args.min_n)),
        "--min-subjects",
        str(int(args.min_subjects)),
        "--n-bootstrap",
        str(int(args.n_bootstrap)),
        "--ci",
        str(float(args.ci)),
        "--seed",
        str(int(args.seed)),
        "--compression",
        str(args.compression),
    ]
    if args.streaming:
        command.append("--streaming")
    if mode != "TT":
        command.append("--allow-nonprimary")
    return command


def build_step13_command(
    args: argparse.Namespace,
    mode: str,
    step11_dir: Path,
    output_dir: Path,
    dt_values: Sequence[int],
) -> List[str]:
    command = [
        str(args.python_executable),
        str(args.step13_script),
        "--input-dir",
        str(step11_dir),
        "--outdir",
        str(output_dir),
        "--dataset-label",
        "{}_step15_{}".format(args.dataset_label, mode.lower()),
        "--required-mode",
        mode,
        "--dt-values",
        ",".join(str(v) for v in dt_values),
        "--metrics",
        "var_dx,msd",
        "--estimators",
        "transition_weighted,equal_subject_weighted",
        "--bin-scope",
        "core",
        "--core-min-subject-fraction",
        str(float(args.core_min_subject_fraction)),
        "--min-dt-points",
        str(int(args.min_dt_points)),
        "--n-bootstrap",
        str(int(args.n_bootstrap)),
        "--seed",
        str(int(args.seed)),
        "--ci",
        str(float(args.ci)),
    ]
    if mode != "TT":
        command.append("--allow-nonprimary")
    return command


def combine_step11_curves(mode_dirs: Dict[str, Path]) -> Optional[pl.DataFrame]:
    frames: List[pl.DataFrame] = []
    for mode, folder in mode_dirs.items():
        table = find_output_table(folder, "10_binned_dynamics_long")
        if table is None:
            continue
        frames.append(read_table(table).with_columns(pl.lit(mode).alias("transition_mode")))
    if not frames:
        return None
    return pl.concat(frames, how="diagonal_relaxed").sort(["transition_mode", "dt", "bin_id"])


def build_fluctuation_mode_difference(
    combined: pl.DataFrame,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    rows: List[pl.DataFrame] = []
    summary_rows: List[Dict[str, Any]] = []

    cols = set(combined.columns)
    metrics = [m for m in CURVE_METRICS if m in cols]
    for metric in metrics:
        tt = combined.filter(pl.col("transition_mode") == "TT").select(
            [
                "dt",
                "bin_id",
                "bin_left",
                "bin_mid",
                "bin_right",
                pl.col(metric).alias("TT"),
            ]
        )
        inc = combined.filter(pl.col("transition_mode") == "NON_FF").select(
            ["dt", "bin_id", pl.col(metric).alias("NON_FF")]
        )
        joined = (
            tt.join(inc, on=["dt", "bin_id"], how="inner")
            .with_columns(
                [
                    pl.lit(metric).alias("metric"),
                    (pl.col("NON_FF") - pl.col("TT")).alias("difference_NON_FF_minus_TT"),
                    (pl.col("NON_FF") - pl.col("TT")).abs().alias("absolute_difference"),
                ]
            )
            .sort(["dt", "bin_id"])
        )
        rows.append(joined)

        for dt in sorted(int(v) for v in joined["dt"].unique().to_list()):
            g = joined.filter(pl.col("dt") == dt)
            a = g["TT"].to_numpy().astype(float)
            b = g["NON_FF"].to_numpy().astype(float)
            d = b - a
            summary_rows.append(
                {
                    "metric": metric,
                    "dt": int(dt),
                    "n_matched_bins": int(np.sum(np.isfinite(a) & np.isfinite(b))),
                    "pearson_r": safe_pearson(a, b),
                    "mean_difference_NON_FF_minus_TT": float(np.nanmean(d)),
                    "mae": float(np.nanmean(np.abs(d))),
                    "rmse": float(np.sqrt(np.nanmean(d ** 2))),
                    "mean_TT": float(np.nanmean(a)),
                    "mean_NON_FF": float(np.nanmean(b)),
                    "mean_ratio_NON_FF_over_TT": (
                        float(np.nanmean(b) / np.nanmean(a))
                        if np.isfinite(np.nanmean(a)) and abs(np.nanmean(a)) > 1e-12
                        else float("nan")
                    ),
                }
            )

    diff = pl.concat(rows, how="diagonal_relaxed") if rows else pl.DataFrame()
    summary = pl.DataFrame(summary_rows) if summary_rows else pl.DataFrame()
    if summary.height:
        summary = summary.sort(["metric", "dt"])
    return diff, summary


def combine_step13_model_tables(mode_dirs: Dict[str, Path]) -> Optional[pl.DataFrame]:
    frames: List[pl.DataFrame] = []
    for mode, folder in mode_dirs.items():
        table = find_output_table(folder, "05_model_comparison_by_bin")
        if table is None:
            continue
        frames.append(read_table(table).with_columns(pl.lit(mode).alias("transition_mode")))
    if not frames:
        return None
    return pl.concat(frames, how="diagonal_relaxed").sort(
        ["metric", "estimator", "transition_mode", "bin_id"]
    )


def cross_mode_temporal_core(model_table: pl.DataFrame) -> pl.DataFrame:
    rows: List[Dict[str, Any]] = []
    metrics = sorted(str(v) for v in model_table["metric"].unique().to_list())
    estimators = sorted(str(v) for v in model_table["estimator"].unique().to_list())

    for metric in metrics:
        for estimator in estimators:
            sub = model_table.filter(
                (pl.col("metric") == metric) & (pl.col("estimator") == estimator)
            )
            mode_bins: Dict[str, set] = {}
            for mode in FLUCTUATION_MODES:
                mode_bins[mode] = set(
                    int(v)
                    for v in sub.filter(pl.col("transition_mode") == mode)["bin_id"].to_list()
                )
            shared = set.intersection(*[mode_bins[m] for m in FLUCTUATION_MODES])
            selected = largest_contiguous_run(shared)
            if not selected:
                continue
            ref = sub.filter(
                (pl.col("transition_mode") == "TT") & pl.col("bin_id").is_in(selected)
            ).sort("bin_id")
            for row in ref.iter_rows(named=True):
                rows.append(
                    {
                        "metric": metric,
                        "estimator": estimator,
                        "bin_id": int(row["bin_id"]),
                        "bin_left": float(row["bin_left"]),
                        "bin_mid": float(row["bin_mid"]),
                        "bin_right": float(row["bin_right"]),
                        "n_shared_bins": int(len(selected)),
                        "core_rule": "largest_contiguous_intersection_of_TT_and_NON_FF_step13_tested_bins",
                    }
                )
    return pl.DataFrame(rows).sort(["metric", "estimator", "bin_id"]) if rows else pl.DataFrame()


def temporal_core_table(model_table: pl.DataFrame, core: pl.DataFrame) -> pl.DataFrame:
    frames: List[pl.DataFrame] = []
    for metric in sorted(str(v) for v in core["metric"].unique().to_list()):
        for estimator in sorted(str(v) for v in core.filter(pl.col("metric") == metric)["estimator"].unique().to_list()):
            bins = core.filter(
                (pl.col("metric") == metric) & (pl.col("estimator") == estimator)
            )["bin_id"].to_list()
            frames.append(
                model_table.filter(
                    (pl.col("metric") == metric)
                    & (pl.col("estimator") == estimator)
                    & pl.col("bin_id").is_in([int(v) for v in bins])
                )
            )
    if not frames:
        return pl.DataFrame()
    return pl.concat(frames, how="diagonal_relaxed").sort(
        ["metric", "estimator", "transition_mode", "bin_id"]
    )


def summarize_temporal_sensitivity(core_table: pl.DataFrame) -> pl.DataFrame:
    rows: List[Dict[str, Any]] = []
    if core_table.height == 0:
        return pl.DataFrame()
    for metric in sorted(str(v) for v in core_table["metric"].unique().to_list()):
        for estimator in sorted(str(v) for v in core_table.filter(pl.col("metric") == metric)["estimator"].unique().to_list()):
            for mode in FLUCTUATION_MODES:
                g = core_table.filter(
                    (pl.col("metric") == metric)
                    & (pl.col("estimator") == estimator)
                    & (pl.col("transition_mode") == mode)
                )
                if g.height == 0:
                    continue
                pref = g["observed_preferred_model"].cast(pl.String).to_list()
                slopes = g["signed_linear_slope_median"].to_numpy().astype(float)
                const_pref = g["constant_preference_fraction"].to_numpy().astype(float)
                ci_below = g["signed_linear_slope_ci_below_zero"].cast(pl.Boolean, strict=False).fill_null(False)
                ci_above = g["signed_linear_slope_ci_above_zero"].cast(pl.Boolean, strict=False).fill_null(False)
                rows.append(
                    {
                        "metric": metric,
                        "estimator": estimator,
                        "transition_mode": mode,
                        "n_common_core_bins": int(g.height),
                        "observed_constant_best_n": int(sum(1 for v in pref if str(v) == "constant")),
                        "median_constant_preference_fraction": float(np.nanmedian(const_pref)),
                        "median_signed_linear_slope": float(np.nanmedian(slopes)),
                        "n_slope_ci_below_zero": int(ci_below.sum()),
                        "n_slope_ci_above_zero": int(ci_above.sum()),
                    }
                )
    return pl.DataFrame(rows).sort(["metric", "estimator", "transition_mode"])


# =============================================================================
# Main
# =============================================================================


def main() -> int:
    args = build_argparser().parse_args()
    if pl is None:
        raise RuntimeError("polars is required to execute Step 15.")
    args.dt_values_parsed = parse_dt_values(args.dt_values)
    validate_args(args)

    outdir = ensure_dir(args.outdir.resolve())
    logs_dir = ensure_dir(outdir / "logs")
    step11_root = ensure_dir(outdir / "step11_by_mode")
    step13_root = ensure_dir(outdir / "step13_by_mode")

    print("=== ClonoDynamics Step 15: detectability-boundary sensitivity ===")
    print("script version:", SCRIPT_VERSION)
    print("threshold labels reclassified: NO")
    print("assumed upstream alpha:", args.detectability_alpha)
    print("forward comparison: TT vs T0PLUS=TT+TF")
    print("fluctuation comparison: TT vs NON_FF=TT+TF+FT")

    manifest = pl.DataFrame(
        [
            {"role": "step3_trajectories", "path": str(args.trajectories.resolve())},
            {"role": "step5_transitions", "path": str(args.transitions.resolve())},
            {"role": "step9_script", "path": str(args.step9_script.resolve())},
            {"role": "step11_script", "path": str(args.step11_script.resolve())},
            {"role": "step13_script", "path": str(args.step13_script.resolve())},
        ]
    )
    manifest.write_csv(outdir / "01_input_manifest.csv")

    # -------------------------------------------------------------------------
    # Step-5 class composition, common xstar support and crossing descriptors.
    # -------------------------------------------------------------------------
    transition_data = prepare_transition_base(
        args.transitions,
        args.dt_values_parsed,
        args.streaming,
    )
    class_by_dt = class_composition_by_dt(transition_data)
    class_by_dt.write_csv(outdir / "02_class_composition_by_dt.csv")

    fluct_support = fluctuation_mode_support(transition_data, args.dt_values_parsed)
    fluct_support.write_csv(outdir / "03_fluctuation_mode_support_by_dt.csv")
    xstar_min, xstar_max = resolve_fluctuation_common_support(fluct_support)
    fluct_support_summary = pl.DataFrame(
        [
            {
                "support_modes": "TT,NON_FF",
                "dt_values": ",".join(str(v) for v in args.dt_values_parsed),
                "xstar_common_min": xstar_min,
                "xstar_common_max": xstar_max,
                "n_bins": int(args.fluctuation_n_bins),
                "support_rule": "intersection_of_exact_mode_x_dt_ranges",
            }
        ]
    )
    fluct_support_summary.write_csv(outdir / "04_fluctuation_common_xstar_support.csv")

    xstar_edges = common_xstar_edges(xstar_min, xstar_max, args.fluctuation_n_bins)
    crossing = build_crossing_metrics(transition_data, xstar_edges)
    write_csv_parquet(crossing, outdir / "05_crossing_metrics_by_dt_bin", args.compression)

    # -------------------------------------------------------------------------
    # Replicate-decoupled forward TT vs T0PLUS sensitivity.
    # -------------------------------------------------------------------------
    step9 = import_step9_module(args.step9_script.resolve())
    (
        fwd_support,
        fwd_edges,
        fwd_curves,
        fwd_diff,
        fwd_summary,
        fwd_audit,
        fwd_decomposition,
        fwd_decomposition_summary,
    ) = build_forward_analysis(
        step9=step9,
        trajectories=args.trajectories.resolve(),
        outdir=outdir,
        dt=int(args.forward_dt),
        n_bins=int(args.forward_n_bins),
        q_low=float(args.forward_q_low),
        q_high=float(args.forward_q_high),
        n_bootstrap=int(args.n_bootstrap),
        seed=int(args.seed),
        ci=float(args.ci),
        compression=str(args.compression),
        streaming=bool(args.streaming),
        restart=bool(args.restart_forward),
        keep_transitions=bool(args.keep_forward_transitions),
    )
    fwd_support.write_csv(outdir / "06_forward_mode_support.csv")
    fwd_edges.write_csv(outdir / "07_forward_bin_edges.csv")
    write_csv_parquet(fwd_curves, outdir / "08_forward_TT_vs_T0plus_curves", args.compression)
    write_csv_parquet(fwd_diff, outdir / "09_forward_T0plus_minus_TT_by_bin", args.compression)
    fwd_summary.write_csv(outdir / "10_forward_sensitivity_summary.csv")
    if fwd_audit:
        pl.DataFrame(fwd_audit).write_csv(outdir / "10a_forward_pair_build_audit.csv")
    write_csv_parquet(
        fwd_decomposition,
        outdir / "10b_forward_TT_TF_decomposition_by_bin",
        args.compression,
    )
    fwd_decomposition_summary.write_csv(
        outdir / "10c_forward_decomposition_summary.csv"
    )

    # -------------------------------------------------------------------------
    # Current Step-11 and Step-13 runs on identical xstar support.
    # -------------------------------------------------------------------------
    step11_dirs: Dict[str, Path] = {}
    step13_dirs: Dict[str, Path] = {}
    step11_command_rows: List[Dict[str, Any]] = []
    step13_command_rows: List[Dict[str, Any]] = []

    for mode in FLUCTUATION_MODES:
        mode_tag = mode.lower()
        s11_dir = ensure_dir(step11_root / mode_tag)
        s13_dir = ensure_dir(step13_root / mode_tag)
        step11_dirs[mode] = s11_dir
        step13_dirs[mode] = s13_dir

        cmd11 = build_step11_command(
            args=args,
            mode=mode,
            output_dir=s11_dir,
            x_min=xstar_min,
            x_max=xstar_max,
            dt_values=args.dt_values_parsed,
        )
        step11_command_rows.append(
            {"transition_mode": mode, "output_dir": str(s11_dir), "command": quote_command(cmd11)}
        )

        expected11 = find_output_table(s11_dir, "10_binned_dynamics_long")
        should_run11 = bool(args.run_step11) and not (args.reuse_existing and expected11 is not None)
        if should_run11:
            run_command(cmd11, logs_dir / "step11_mode_{}.log".format(mode_tag))

        cmd13 = build_step13_command(
            args=args,
            mode=mode,
            step11_dir=s11_dir,
            output_dir=s13_dir,
            dt_values=args.dt_values_parsed,
        )
        step13_command_rows.append(
            {"transition_mode": mode, "output_dir": str(s13_dir), "command": quote_command(cmd13)}
        )

        expected13 = find_output_table(s13_dir, "05_model_comparison_by_bin")
        should_run13 = bool(args.run_step13) and not (args.reuse_existing and expected13 is not None)
        if should_run13:
            if find_output_table(s11_dir, "10_binned_dynamics_long") is None:
                raise RuntimeError(
                    "Step 13 requested for {} but Step-11 output is unavailable in {}".format(
                        mode, s11_dir
                    )
                )
            run_command(cmd13, logs_dir / "step13_mode_{}.log".format(mode_tag))

    pl.DataFrame(step11_command_rows).write_csv(outdir / "11_step11_commands.csv")
    pl.DataFrame(step13_command_rows).write_csv(outdir / "12_step13_commands.csv")

    combined_step11 = combine_step11_curves(step11_dirs)
    fluct_diff = None
    fluct_summary = None
    if combined_step11 is not None:
        write_csv_parquet(
            combined_step11,
            outdir / "13_combined_step11_fluctuation_curves",
            args.compression,
        )
        fluct_diff, fluct_summary = build_fluctuation_mode_difference(combined_step11)
        if fluct_diff.height:
            write_csv_parquet(
                fluct_diff,
                outdir / "14_fluctuation_mode_difference_vs_TT",
                args.compression,
            )
        if fluct_summary.height:
            fluct_summary.write_csv(outdir / "15_fluctuation_mode_robustness_summary.csv")

    combined_step13 = combine_step13_model_tables(step13_dirs)
    temporal_core = None
    temporal_core_table_df = None
    temporal_summary = None
    if combined_step13 is not None:
        write_csv_parquet(
            combined_step13,
            outdir / "16_combined_step13_model_comparison_by_bin",
            args.compression,
        )
        temporal_core = cross_mode_temporal_core(combined_step13)
        if temporal_core.height:
            temporal_core.write_csv(outdir / "17_cross_mode_common_temporal_core.csv")
            temporal_core_table_df = temporal_core_table(combined_step13, temporal_core)
            write_csv_parquet(
                temporal_core_table_df,
                outdir / "18_temporal_scaling_mode_comparison_core",
                args.compression,
            )
            temporal_summary = summarize_temporal_sensitivity(temporal_core_table_df)
            if temporal_summary.height:
                temporal_summary.write_csv(outdir / "19_temporal_scaling_sensitivity_summary.csv")

    # -------------------------------------------------------------------------
    # Compact Step-15 summary for manuscript triage.
    # -------------------------------------------------------------------------
    summary_rows: List[Dict[str, Any]] = []
    fwd_row = fwd_summary.row(0, named=True)
    for key in [
        "pearson_r_curve",
        "mean_abs_difference",
        "rms_difference",
        "tt_signed_slope",
        "t0plus_signed_slope",
        "mean_abs_drift_TT",
        "mean_abs_drift_T0PLUS",
    ]:
        summary_rows.append(
            {"component": "forward_T0PLUS_vs_TT", "metric": key, "value": fwd_row.get(key)}
        )

    decomp_row = fwd_decomposition_summary.row(0, named=True)
    for key in [
        "mean_p_TF_given_T0_detectable",
        "median_p_TF_given_T0_detectable",
        "mean_dx_TF_across_bins",
        "tf_signed_slope",
        "mean_abs_TT_contribution",
        "mean_abs_TF_contribution",
        "mean_abs_reconstruction_error",
        "max_abs_reconstruction_error",
        "pearson_direct_vs_reconstructed",
        "max_abs_count_identity_residual",
    ]:
        summary_rows.append(
            {
                "component": "forward_TT_TF_decomposition",
                "metric": key,
                "value": decomp_row.get(key),
            }
        )

    # Overall crossing descriptors by dt, weighted by transition counts.
    for dt in args.dt_values_parsed:
        g = transition_data.filter(pl.col("dt") == int(dt))
        n_all = int(g.height)
        n_tt = int((g["obs_class"] == "TT").sum())
        n_tf = int((g["obs_class"] == "TF").sum())
        n_ft = int((g["obs_class"] == "FT").sum())
        n_ff = int((g["obs_class"] == "FF").sum())
        n_cross = n_tf + n_ft
        n_nonff = n_tt + n_tf + n_ft
        summary_rows.extend(
            [
                {
                    "component": "detectability_crossing_dt{}".format(dt),
                    "metric": "crossing_burden_all",
                    "value": float(n_cross / n_all) if n_all else float("nan"),
                },
                {
                    "component": "detectability_crossing_dt{}".format(dt),
                    "metric": "crossing_burden_non_ff",
                    "value": float(n_cross / n_nonff) if n_nonff else float("nan"),
                },
                {
                    "component": "detectability_crossing_dt{}".format(dt),
                    "metric": "crossing_imbalance",
                    "value": float((n_ft - n_tf) / n_cross) if n_cross else float("nan"),
                },
                {
                    "component": "detectability_crossing_dt{}".format(dt),
                    "metric": "tt_fraction_all",
                    "value": float(n_tt / n_all) if n_all else float("nan"),
                },
                {
                    "component": "detectability_crossing_dt{}".format(dt),
                    "metric": "ff_fraction_all",
                    "value": float(n_ff / n_all) if n_all else float("nan"),
                },
            ]
        )

    if fluct_summary is not None and fluct_summary.height:
        for row in fluct_summary.iter_rows(named=True):
            if row["metric"] not in {"var_dx", "msd"}:
                continue
            summary_rows.append(
                {
                    "component": "fluctuation_NON_FF_vs_TT_dt{}".format(row["dt"]),
                    "metric": "{}_pearson_r".format(row["metric"]),
                    "value": row["pearson_r"],
                }
            )
            summary_rows.append(
                {
                    "component": "fluctuation_NON_FF_vs_TT_dt{}".format(row["dt"]),
                    "metric": "{}_mean_ratio_NON_FF_over_TT".format(row["metric"]),
                    "value": row["mean_ratio_NON_FF_over_TT"],
                }
            )

    if temporal_summary is not None and temporal_summary.height:
        primary = temporal_summary.filter(pl.col("estimator") == "transition_weighted")
        for row in primary.iter_rows(named=True):
            summary_rows.extend(
                [
                    {
                        "component": "temporal_scaling_{}_{}".format(
                            row["metric"], row["transition_mode"]
                        ),
                        "metric": "median_constant_preference_fraction",
                        "value": row["median_constant_preference_fraction"],
                    },
                    {
                        "component": "temporal_scaling_{}_{}".format(
                            row["metric"], row["transition_mode"]
                        ),
                        "metric": "median_signed_linear_slope",
                        "value": row["median_signed_linear_slope"],
                    },
                    {
                        "component": "temporal_scaling_{}_{}".format(
                            row["metric"], row["transition_mode"]
                        ),
                        "metric": "n_slope_ci_below_zero",
                        "value": row["n_slope_ci_below_zero"],
                    },
                    {
                        "component": "temporal_scaling_{}_{}".format(
                            row["metric"], row["transition_mode"]
                        ),
                        "metric": "n_slope_ci_above_zero",
                        "value": row["n_slope_ci_above_zero"],
                    },
                ]
            )

    step15_summary = pl.DataFrame(summary_rows)
    step15_summary.write_csv(outdir / "20_step15_summary.csv")

    config: Dict[str, Any] = {
        "step": 15,
        "script": "15-detectability_boundary_sensitivity.py",
        "script_version": SCRIPT_VERSION,
        "upstream_steps": [4, 5, 9, 11, 13],
        "timestamp": datetime.now().isoformat(),
        "dataset_label": args.dataset_label,
        "trajectories": str(args.trajectories.resolve()),
        "transitions": str(args.transitions.resolve()),
        "detectability_alpha_assumed": float(args.detectability_alpha),
        "transition_labels_reclassified": False,
        "alpha_sweep_performed": False,
        "forward": {
            "dt": int(args.forward_dt),
            "primary_modes": list(FORWARD_MODES),
            "T0PLUS_definition": "observable at t0; t1 unrestricted = TT+TF",
            "TF_definition": "observable at t0; non-observable at t1",
            "conditioning": "replicate-decoupled initial-state x0",
            "AB": "xcond=rep1 posterior median at t0; dx=rep2 posterior mean(t1)-mean(t0)",
            "BA": "xcond=rep2 posterior median at t0; dx=rep1 posterior mean(t1)-mean(t0)",
            "fold_combination": "equal AB/BA weight; folds are not independent N",
            "support_rule": "intersection of TT and T0PLUS x0 quantile ranges",
            "q_low": float(args.forward_q_low),
            "q_high": float(args.forward_q_high),
            "bin_edges": "TT reference quantiles within common support",
            "requested_n_bins": int(args.forward_n_bins),
            "decomposition": (
                "TF is derived algebraically as the exact T0PLUS-minus-TT remainder within "
                "each fold and bin; AB/BA are then combined with equal fold weight"
            ),
        },
        "fluctuations": {
            "modes": list(FLUCTUATION_MODES),
            "NON_FF_definition": "TT+TF+FT; FF excluded",
            "representation": "dx_latent",
            "conditioning": "xstar_latent",
            "primary_metric": "var_dx",
            "complementary_metric": "msd",
            "dt_values": list(args.dt_values_parsed),
            "common_xstar_min": float(xstar_min),
            "common_xstar_max": float(xstar_max),
            "n_bins": int(args.fluctuation_n_bins),
        },
        "crossing_descriptors": {
            "C_all": "(n_TF+n_FT)/n_all",
            "C_nonFF": "(n_TF+n_FT)/(n_TT+n_TF+n_FT)",
            "J": "(n_FT-n_TF)/(n_TF+n_FT)",
        },
        "temporal_scaling": {
            "current_step13_reused": True,
            "cross_mode_core_rule": "largest contiguous intersection of Step-13 tested bins",
            "paired_cross_mode_parameter_test": False,
        },
        "n_bootstrap": int(args.n_bootstrap),
        "seed": int(args.seed),
        "ci": float(args.ci),
        "run_step11": bool(args.run_step11),
        "run_step13": bool(args.run_step13),
        "reuse_existing": bool(args.reuse_existing),
        "step9_script": str(args.step9_script.resolve()),
        "step11_script": str(args.step11_script.resolve()),
        "step13_script": str(args.step13_script.resolve()),
    }
    with open(outdir / "00_run_config.json", "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)

    with open(outdir / "DEBUG_REPORT.txt", "w", encoding="utf-8") as handle:
        handle.write("STEP 15 — DETECTABILITY-BOUNDARY SENSITIVITY\n")
        handle.write("===========================================\n\n")
        handle.write("Script version: {}\n".format(SCRIPT_VERSION))
        handle.write("Dataset: {}\n".format(args.dataset_label))
        handle.write("Assumed upstream detectability alpha: {}\n".format(args.detectability_alpha))
        handle.write("Endpoint labels reclassified: NO\n\n")
        handle.write("FORWARD\n")
        handle.write("  TT vs T0PLUS=TT+TF\n")
        handle.write("  dt={}\n".format(args.forward_dt))
        handle.write("  replicate-decoupled AB/BA geometry imported from current Step 9\n")
        handle.write("  decomposition: TF derived as exact T0PLUS-minus-TT remainder within each fold\n")
        handle.write("  common x0 support: [{:.8g}, {:.8g}]\n".format(
            float(fwd_support["common_x0_min"][0]), float(fwd_support["common_x0_max"][0])
        ))
        handle.write("  actual bins: {}\n\n".format(fwd_edges.height - 1))
        handle.write("FLUCTUATIONS\n")
        handle.write("  TT vs NON_FF=TT+TF+FT\n")
        handle.write("  common xstar support: [{:.8g}, {:.8g}]\n".format(xstar_min, xstar_max))
        handle.write("  bins: {}\n".format(args.fluctuation_n_bins))
        handle.write("  dt: {}\n".format(args.dt_values_parsed))
        handle.write("  Step 11 run requested: {}\n".format(args.run_step11))
        handle.write("  Step 13 run requested: {}\n\n".format(args.run_step13))
        handle.write("CROSSING DESCRIPTORS\n")
        handle.write("  C_all=(TF+FT)/ALL\n")
        handle.write("  C_nonFF=(TF+FT)/(TT+TF+FT)\n")
        handle.write("  J=(FT-TF)/(TF+FT), computed after signed integer casting\n\n")
        if temporal_core is not None and temporal_core.height:
            handle.write("CROSS-MODE TEMPORAL CORE\n")
            for metric in sorted(str(v) for v in temporal_core["metric"].unique().to_list()):
                p = temporal_core.filter(
                    (pl.col("metric") == metric)
                    & (pl.col("estimator") == "transition_weighted")
                )
                if p.height:
                    handle.write(
                        "  {}: bins={} xstar=[{:.8g}, {:.8g}]\n".format(
                            metric,
                            [int(v) for v in p["bin_id"].to_list()],
                            float(p["bin_left"].min()),
                            float(p["bin_right"].max()),
                        )
                    )

    print("[DONE]", outdir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("ERROR:", exc, file=sys.stderr)
        raise
