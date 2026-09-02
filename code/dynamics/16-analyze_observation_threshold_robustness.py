#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
16-analyze_observation_threshold_robustness.py
===============================================

Cross-alpha sensitivity analysis for the operational observation-state threshold
used in ClonoDynamics.

Scientific purpose
------------------
The primary operational observation state is defined from the empirical endpoint
p-value as:

    T(alpha) = 1[p < alpha]
    F(alpha) = 1[p >= alpha]

The primary manuscript threshold is alpha = 0.05. This analyzer evaluates:

    alpha = 0.10, 0.05, 0.025, 0.01

without refitting latent abundance. Only the operational T/F state is changed.

Crucially, this script DOES NOT reimplement the mathematics of the principal
dynamical estimators. It creates threshold-specific trajectory/transition tables
and invokes the current canonical ClonoDynamics analyses:

    P09  9-replicate_decoupled_forward_drift.py
         forward drift from threshold-specific Step-4 trajectories

    P11  11-longitudinal_fluctuation_dynamics.py
         fluctuation magnitude from threshold-specific Step-5 transitions

    P13  13-temporal_fluctuation_scaling.py
         finite-lag temporal scaling from the corresponding P11 outputs

Thus alpha changes the operational observation domain while the estimators,
bootstrap logic and temporal-model mathematics remain those of P09/P11/P13.

Primary design
--------------
1. Obtain one empirical p-value per subject x clonotype x time endpoint from
   Step 2, using minimum p when multiple rows contribute to the endpoint.
2. Reclassify both:
      - the Step-4 trajectory `observable` field;
      - the Step-5 transition obs0/obs1/obs_class fields.
3. Keep all latent states/posterior summaries unchanged.
4. Define one common xstar interval across all alpha x dt TT transition strata
   for P11.
5. Use one common P09 forward-drift bin definition across alpha:
      - reference alpha 0.05 is run first;
      - its generated forward_drift_bin_edges.csv is reused by all other alpha.
6. Invoke P11 for each alpha on identical xstar limits and absolute bins.
7. Invoke P13 for each alpha using the canonical six-model implementation.
8. Combine outputs and resolve the largest contiguous cross-alpha temporal core
   from bins actually tested by P13 across all alpha/metric/estimator strata.

Interpretation boundary
-----------------------
This is sensitivity to an OPERATIONAL OBSERVATION-STATE CLASSIFICATION threshold.
It is not a sensitivity analysis of posterior-predictive detectability itself,
and T/F must not be interpreted as biological presence/absence.

Python >= 3.9
Required: polars, numpy
"""

from __future__ import annotations

import argparse
import json
import math
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl


VERSION = "2.0.0-p09-p11-p13"
DEFAULT_ALPHAS = (0.10, 0.05, 0.025, 0.01)
DEFAULT_DT = (1, 2, 3, 4, 5)
REFERENCE_ALPHA = 0.05

P11_CURVE_METRICS = ("mean_dx", "median_dx", "ppos", "msd", "var_dx", "mad_dx")
P09_CURVE_METRICS = ("mean_dx", "median_dx", "p_dx_gt0")


# =============================================================================
# CLI
# =============================================================================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Reclassify operational observation states across alpha and reuse "
            "canonical P09/P11/P13 analyses without reimplementing their mathematics."
        ),
    )

    p.add_argument(
        "--trajectories",
        required=True,
        type=Path,
        help="Current P04 latent_trajectories_long.parquet used by P09.",
    )
    p.add_argument(
        "--transitions",
        required=True,
        type=Path,
        help="Current P05 latent transition table used by P11.",
    )

    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--observability",
        type=Path,
        default=None,
        help=(
            "Preferred Step-2 endpoint observability table containing "
            "subject, aaSeqCDR3, time and p_value."
        ),
    )
    source.add_argument(
        "--step2-dir",
        type=Path,
        default=None,
        help=(
            "Step-2 output directory. The analyzer searches "
            "latent_trajectory_observability_subject/global/clone and then "
            "per_clone_latent_subject/global/clone."
        ),
    )
    source.add_argument(
        "--per-clone",
        type=Path,
        default=None,
        help="Explicit Step-2 table containing endpoint empirical p-values.",
    )
    source.add_argument(
        "--per-clone-dir",
        type=Path,
        default=None,
        help="Directory containing Step-2 table parts with empirical p-values.",
    )

    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--dataset-label", default="healthy")

    p.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=list(DEFAULT_ALPHAS),
        help="Operational observation-state alpha values.",
    )
    p.add_argument(
        "--reference-alpha",
        type=float,
        default=REFERENCE_ALPHA,
        help="Reference alpha for agreement and cross-alpha comparisons.",
    )
    p.add_argument(
        "--pvalue-col",
        default="auto",
        help="'auto' prefers p_value, then p_emp_low.",
    )
    p.add_argument("--subject-col", default="subject")
    p.add_argument("--clone-col", default="aaSeqCDR3")
    p.add_argument("--time-col", default="time")
    p.add_argument(
        "--allow-missing-pvalues",
        action="store_true",
        help=(
            "Permit missing endpoint p-values and classify them as F. "
            "Default behavior is to fail because reclassification would be ambiguous."
        ),
    )

    # Shared longitudinal settings.
    p.add_argument("--dt", type=int, nargs="*", default=list(DEFAULT_DT))
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--ci", type=float, default=95.0)
    p.add_argument("--streaming", action="store_true")
    p.add_argument("--compression", default="zstd")
    p.add_argument("--python-executable", default=sys.executable)

    # P09 settings.
    p.add_argument(
        "--p09-script",
        type=Path,
        default=Path(__file__).resolve().parent / "9-replicate_decoupled_forward_drift.py",
    )
    p.add_argument("--p09-dt", type=int, default=1)
    p.add_argument("--p09-n-bins", type=int, default=20)
    p.add_argument(
        "--p09-binning",
        choices=["quantile", "equal_width"],
        default="quantile",
    )
    p.add_argument("--p09-q-low", type=float, default=0.005)
    p.add_argument("--p09-q-high", type=float, default=0.995)
    p.add_argument("--p09-n-bootstrap", type=int, default=2000)
    p.add_argument("--p09-row-group-size", type=int, default=100000)
    p.add_argument(
        "--p09-bin-edges-file",
        type=Path,
        default=None,
        help=(
            "Optional externally fixed P09 edges. If omitted, the reference-alpha "
            "P09 run defines the common edges and all other alphas reuse them."
        ),
    )

    # P11 settings.
    p.add_argument(
        "--p11-script",
        type=Path,
        default=Path(__file__).resolve().parent / "11-longitudinal_fluctuation_dynamics.py",
    )
    p.add_argument("--p11-n-bins", type=int, default=30)
    p.add_argument("--p11-min-n", type=int, default=30)
    p.add_argument("--p11-min-subjects", type=int, default=2)
    p.add_argument("--p11-n-bootstrap", type=int, default=2000)

    # P13 settings.
    p.add_argument(
        "--p13-script",
        type=Path,
        default=Path(__file__).resolve().parent / "13-temporal_fluctuation_scaling.py",
    )
    p.add_argument("--p13-n-bootstrap", type=int, default=2000)
    p.add_argument(
        "--p13-core-min-subject-fraction",
        type=float,
        default=1.0,
    )
    p.add_argument("--p13-core-min-subjects", type=int, default=None)
    p.add_argument(
        "--p13-bin-scope",
        choices=["core", "all_complete"],
        default="core",
    )

    # Execution control.
    p.add_argument("--run-p09", action="store_true")
    p.add_argument("--run-p11", action="store_true")
    p.add_argument(
        "--run-p13",
        action="store_true",
        help="Run P13 for every alpha; implies --run-p11.",
    )
    p.add_argument(
        "--run-all",
        action="store_true",
        help="Run P09, P11 and P13 for every alpha.",
    )
    p.add_argument(
        "--restart-p09",
        action="store_true",
        help="Pass --restart to every P09 run.",
    )
    return p


def validate_args(args: argparse.Namespace) -> List[float]:
    alphas = sorted(set(float(v) for v in args.alphas), reverse=True)
    if not alphas:
        raise ValueError("At least one alpha is required.")
    if any((not math.isfinite(v) or v <= 0.0 or v >= 1.0) for v in alphas):
        raise ValueError("Every alpha must be finite and in (0,1).")
    if not any(abs(v - float(args.reference_alpha)) < 1e-12 for v in alphas):
        raise ValueError("--reference-alpha must be included in --alphas.")

    if args.run_all:
        args.run_p09 = True
        args.run_p11 = True
        args.run_p13 = True
    if args.run_p13:
        args.run_p11 = True

    if not args.dt:
        raise ValueError("At least one --dt value is required.")
    args.dt = sorted(set(int(v) for v in args.dt))
    if any(v <= 0 for v in args.dt):
        raise ValueError("All dt values must be positive.")

    if args.p09_dt <= 0:
        raise ValueError("--p09-dt must be positive.")
    if args.p09_n_bins < 3:
        raise ValueError("--p09-n-bins must be >=3.")
    if not (0 <= args.p09_q_low < args.p09_q_high <= 1):
        raise ValueError("Require 0 <= p09-q-low < p09-q-high <= 1.")
    if args.p09_n_bootstrap <= 0:
        raise ValueError("--p09-n-bootstrap must be >0.")

    if args.p11_n_bins < 2:
        raise ValueError("--p11-n-bins must be >=2.")
    if args.p11_min_n < 1 or args.p11_min_subjects < 1:
        raise ValueError("P11 support thresholds must be >=1.")
    if args.p11_n_bootstrap <= 0:
        raise ValueError("--p11-n-bootstrap must be >0.")

    if args.p13_n_bootstrap <= 0:
        raise ValueError("--p13-n-bootstrap must be >0.")
    if not (0 < args.p13_core_min_subject_fraction <= 1):
        raise ValueError("--p13-core-min-subject-fraction must be in (0,1].")

    return alphas


# =============================================================================
# Generic utilities
# =============================================================================

def alpha_tag(alpha: float) -> str:
    text = "{:.8g}".format(float(alpha)).rstrip("0").rstrip(".")
    return text.replace("-", "m").replace(".", "p")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def quote_command(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in command)


def run_command(command: Sequence[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("[RUN] {}".format(quote_command(command)))
    with open(log_path, "w", encoding="utf-8") as handle:
        result = subprocess.run(
            list(command),
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code {}. See {}".format(
                result.returncode, log_path
            )
        )


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
    raise ValueError("Unsupported table format: {}".format(path))


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
    raise ValueError("Unsupported table format: {}".format(path))


def schema_names(lf: pl.LazyFrame) -> List[str]:
    return list(lf.collect_schema().names())


def write_csv(df: pl.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path)


def write_csv_parquet(df: pl.DataFrame, stem: Path, compression: str) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(stem.with_suffix(".csv"))
    df.write_parquet(stem.with_suffix(".parquet"), compression=compression)


def sink_parquet(
    lf: pl.LazyFrame,
    path: Path,
    compression: str,
    streaming: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        lf.sink_parquet(str(path), compression=compression)
    except Exception:
        collect_frame(lf, streaming).write_parquet(
            str(path), compression=compression
        )


def find_output_table(folder: Path, stem: str) -> Optional[Path]:
    for suffix in (".parquet", ".pq", ".csv", ".tsv", ".txt"):
        candidate = folder / (stem + suffix)
        if candidate.exists():
            return candidate
    return None


def largest_contiguous_run(bin_ids: Iterable[int]) -> List[int]:
    values = sorted(set(int(v) for v in bin_ids))
    if not values:
        return []
    runs: List[List[int]] = []
    current = [values[0]]
    for value in values[1:]:
        if value == current[-1] + 1:
            current.append(value)
        else:
            runs.append(current)
            current = [value]
    runs.append(current)
    runs.sort(key=lambda run: (len(run), -run[0]), reverse=True)
    return runs[0]


# =============================================================================
# Step-2 p-value source
# =============================================================================

def choose_pvalue_column(columns: Sequence[str], requested: str) -> str:
    available = set(columns)
    if requested != "auto":
        if requested not in available:
            raise ValueError(
                "Requested p-value column {!r} is absent.".format(requested)
            )
        return requested
    for candidate in ("p_value", "p_emp_low"):
        if candidate in available:
            return candidate
    raise ValueError(
        "No empirical endpoint p-value column found; expected p_value or p_emp_low."
    )


def scan_table_source(
    file_path: Optional[Path],
    directory: Optional[Path],
) -> Tuple[pl.LazyFrame, str]:
    if file_path is not None:
        if not file_path.exists():
            raise FileNotFoundError(file_path)
        return scan_table(file_path), str(file_path)

    if directory is None or not directory.is_dir():
        raise FileNotFoundError(directory)

    parquet = sorted(list(directory.glob("*.parquet")) + list(directory.glob("*.pq")))
    if parquet:
        return pl.scan_parquet([str(p) for p in parquet]), "{} ({} parquet parts)".format(directory, len(parquet))

    text = sorted(list(directory.glob("*.csv")) + list(directory.glob("*.tsv")) + list(directory.glob("*.txt")))
    if text:
        frames = [scan_table(p) for p in text]
        return pl.concat(frames, how="vertical_relaxed"), "{} ({} text parts)".format(directory, len(text))

    ipc = sorted(list(directory.glob("*.feather")) + list(directory.glob("*.arrow")) + list(directory.glob("*.ipc")))
    if ipc:
        frames = [scan_table(p) for p in ipc]
        return pl.concat(frames, how="vertical_relaxed"), "{} ({} IPC parts)".format(directory, len(ipc))

    raise FileNotFoundError("No supported tables found in {}".format(directory))


def resolve_pvalue_source(
    args: argparse.Namespace,
) -> Tuple[pl.LazyFrame, str, str, List[str]]:
    candidates: List[Tuple[Optional[Path], Optional[Path], str]] = []

    if args.observability is not None:
        candidates.append((args.observability, None, "explicit observability"))
    if args.step2_dir is not None:
        root = args.step2_dir
        for null_name in ("subject", "global", "clone"):
            candidates.extend(
                [
                    (
                        root / "latent_trajectory_observability_{}.csv".format(null_name),
                        None,
                        "Step-2 observability",
                    ),
                    (
                        root / "per_clone_latent_{}.parquet".format(null_name),
                        None,
                        "Step-2 aggregate latent",
                    ),
                    (
                        root / "per_clone_latent_{}.csv".format(null_name),
                        None,
                        "Step-2 aggregate latent",
                    ),
                ]
            )
    if args.per_clone is not None or args.per_clone_dir is not None:
        candidates.append((args.per_clone, args.per_clone_dir, "explicit fallback"))

    inspected = []
    for file_path, directory, label in candidates:
        if file_path is not None and not file_path.exists():
            continue
        if directory is not None and not directory.exists():
            continue
        try:
            lf, source_label = scan_table_source(file_path, directory)
            cols = schema_names(lf)
            required = {args.subject_col, args.clone_col, args.time_col}
            if not required.issubset(set(cols)):
                inspected.append((source_label, cols))
                continue
            try:
                pcol = choose_pvalue_column(cols, args.pvalue_col)
            except ValueError:
                inspected.append((source_label, cols))
                continue
            print("[PVALUES] {}: {}".format(label, source_label))
            print("[PVALUES] column={}".format(pcol))
            return lf, source_label, pcol, cols
        except Exception:
            continue

    detail = "\n".join(
        "  {} columns={}".format(label, cols) for label, cols in inspected[:5]
    )
    raise ValueError(
        "Could not resolve a Step-2 empirical endpoint p-value source."
        + (("\nInspected:\n" + detail) if detail else "")
    )


def endpoint_pvalues(
    source_lf: pl.LazyFrame,
    pvalue_col: str,
    args: argparse.Namespace,
) -> pl.LazyFrame:
    return (
        source_lf.select(
            [
                pl.col(args.subject_col).cast(pl.String, strict=False).alias("subject"),
                pl.col(args.clone_col).cast(pl.String, strict=False).alias("aaSeqCDR3"),
                pl.col(args.time_col).cast(pl.Int64, strict=False).alias("time"),
                pl.col(pvalue_col).cast(pl.Float64, strict=False).alias("_p_raw"),
            ]
        )
        .filter(
            pl.col("subject").is_not_null()
            & pl.col("aaSeqCDR3").is_not_null()
            & pl.col("time").is_not_null()
            & pl.col("_p_raw").is_not_null()
            & pl.col("_p_raw").is_finite()
        )
        .group_by(["subject", "aaSeqCDR3", "time"])
        .agg(
            [
                pl.col("_p_raw").min().alias("observation_p_value"),
                pl.len().alias("n_step2_rows_endpoint"),
            ]
        )
    )


# =============================================================================
# Threshold-specific trajectories
# =============================================================================

def build_enriched_trajectories(
    trajectories: Path,
    endpoint_p: pl.LazyFrame,
    out_path: Path,
    compression: str,
    streaming: bool,
) -> Tuple[pl.LazyFrame, pl.DataFrame]:
    lf = scan_table(trajectories)
    names = schema_names(lf)
    required = {"subject", "aaSeqCDR3", "time", "observable"}
    missing = sorted(required - set(names))
    if missing:
        raise ValueError("P04 trajectory table missing: {}".format(missing))

    base = (
        lf.with_columns(
            [
                pl.col("subject").cast(pl.String, strict=False).alias("subject"),
                pl.col("aaSeqCDR3").cast(pl.String, strict=False).alias("aaSeqCDR3"),
                pl.col("time").cast(pl.Int64, strict=False).alias("time"),
                pl.col("observable").alias("observable_original"),
            ]
        )
        .join(endpoint_p, on=["subject", "aaSeqCDR3", "time"], how="left")
    )

    audit = collect_frame(
        base.select(
            [
                pl.len().alias("n_states"),
                pl.col("observation_p_value").is_null().sum().alias("n_missing_p"),
            ]
        ).with_columns(
            (pl.col("n_missing_p") / pl.col("n_states")).alias("fraction_missing_p")
        ),
        streaming,
    )

    sink_parquet(base, out_path, compression, streaming)
    return pl.scan_parquet(str(out_path)), audit


def trajectory_classified(
    enriched_trajectories: pl.LazyFrame,
    alpha: float,
) -> pl.LazyFrame:
    return enriched_trajectories.with_columns(
        [
            (pl.col("observation_p_value") < float(alpha))
            .fill_null(False)
            .alias("observable"),
            pl.lit(float(alpha)).alias("observation_alpha"),
        ]
    )


# =============================================================================
# Threshold-specific transitions
# =============================================================================

def build_enriched_transitions(
    transitions: Path,
    endpoint_p: pl.LazyFrame,
    out_path: Path,
    compression: str,
    streaming: bool,
) -> Tuple[pl.LazyFrame, pl.DataFrame]:
    lf = scan_table(transitions)
    names = schema_names(lf)
    required = {"subject", "aaSeqCDR3", "t0", "t1", "xstar_latent"}
    missing = sorted(required - set(names))
    if missing:
        raise ValueError("P05 transition table missing: {}".format(missing))

    base = lf.with_columns(
        [
            pl.col("subject").cast(pl.String, strict=False).alias("subject"),
            pl.col("aaSeqCDR3").cast(pl.String, strict=False).alias("aaSeqCDR3"),
            pl.col("t0").cast(pl.Int64, strict=False).alias("t0"),
            pl.col("t1").cast(pl.Int64, strict=False).alias("t1"),
        ]
    )
    if "obs0" in names:
        base = base.with_columns(pl.col("obs0").alias("obs0_original"))
    if "obs1" in names:
        base = base.with_columns(pl.col("obs1").alias("obs1_original"))
    if "obs_class" in names:
        base = base.with_columns(pl.col("obs_class").alias("obs_class_original"))

    p0 = endpoint_p.rename(
        {
            "time": "t0",
            "observation_p_value": "observation_p_value_t0",
            "n_step2_rows_endpoint": "n_step2_rows_t0",
        }
    )
    p1 = endpoint_p.rename(
        {
            "time": "t1",
            "observation_p_value": "observation_p_value_t1",
            "n_step2_rows_endpoint": "n_step2_rows_t1",
        }
    )
    base = (
        base.join(p0, on=["subject", "aaSeqCDR3", "t0"], how="left")
        .join(p1, on=["subject", "aaSeqCDR3", "t1"], how="left")
    )

    audit = collect_frame(
        base.select(
            [
                pl.len().alias("n_transitions"),
                pl.col("observation_p_value_t0").is_null().sum().alias("n_missing_p_t0"),
                pl.col("observation_p_value_t1").is_null().sum().alias("n_missing_p_t1"),
                (
                    pl.col("observation_p_value_t0").is_null()
                    | pl.col("observation_p_value_t1").is_null()
                ).sum().alias("n_missing_either_endpoint"),
            ]
        ).with_columns(
            (
                pl.col("n_missing_either_endpoint") / pl.col("n_transitions")
            ).alias("fraction_missing_either_endpoint")
        ),
        streaming,
    )

    sink_parquet(base, out_path, compression, streaming)
    return pl.scan_parquet(str(out_path)), audit


def transition_classified(
    enriched_transitions: pl.LazyFrame,
    alpha: float,
) -> pl.LazyFrame:
    t0 = (pl.col("observation_p_value_t0") < float(alpha)).fill_null(False)
    t1 = (pl.col("observation_p_value_t1") < float(alpha)).fill_null(False)
    return enriched_transitions.with_columns(
        [
            pl.when(t0).then(pl.lit("T")).otherwise(pl.lit("F")).alias("obs0"),
            pl.when(t1).then(pl.lit("T")).otherwise(pl.lit("F")).alias("obs1"),
            (
                pl.when(t0 & t1).then(pl.lit("TT"))
                .when(t0 & ~t1).then(pl.lit("TF"))
                .when(~t0 & t1).then(pl.lit("FT"))
                .otherwise(pl.lit("FF"))
                .alias("obs_class")
            ),
            pl.lit(float(alpha)).alias("observation_alpha"),
        ]
    )


def class_composition(
    classified: pl.LazyFrame,
    alpha: float,
    selected_dt: Sequence[int],
    streaming: bool,
) -> pl.DataFrame:
    names = schema_names(classified)
    dt_expr = (
        pl.col("dt").cast(pl.Int64, strict=False)
        if "dt" in names
        else pl.col("t1").cast(pl.Int64) - pl.col("t0").cast(pl.Int64)
    )
    data = classified.with_columns(dt_expr.alias("_dt")).filter(
        pl.col("_dt").is_in(list(selected_dt))
    )
    counts = collect_frame(
        data.group_by(["_dt", "obs_class"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.struct(["subject", "aaSeqCDR3"]).n_unique().alias("n_subject_clonotype_pairs"),
            ]
        ),
        streaming,
    ).rename({"_dt": "dt"})
    totals = counts.group_by("dt").agg(
        pl.col("n_transitions").sum().alias("n_dt_total")
    )
    return (
        counts.join(totals, on="dt", how="left")
        .with_columns(
            [
                (pl.col("n_transitions") / pl.col("n_dt_total")).alias("fraction_within_dt"),
                pl.lit(float(alpha)).alias("alpha"),
                pl.lit(alpha_tag(alpha)).alias("alpha_tag"),
            ]
        )
        .sort(["dt", "obs_class"])
    )


def tt_support(
    classified: pl.LazyFrame,
    alpha: float,
    selected_dt: Sequence[int],
    streaming: bool,
) -> pl.DataFrame:
    names = schema_names(classified)
    dt_expr = (
        pl.col("dt").cast(pl.Int64, strict=False)
        if "dt" in names
        else pl.col("t1").cast(pl.Int64) - pl.col("t0").cast(pl.Int64)
    )
    tt = (
        classified.with_columns(dt_expr.alias("_dt"))
        .filter(
            pl.col("_dt").is_in(list(selected_dt))
            & (pl.col("obs_class") == "TT")
            & pl.col("xstar_latent").cast(pl.Float64, strict=False).is_finite()
        )
    )
    out = collect_frame(
        tt.group_by("_dt").agg(
            [
                pl.col("xstar_latent").cast(pl.Float64).min().alias("xstar_min"),
                pl.col("xstar_latent").cast(pl.Float64).max().alias("xstar_max"),
                pl.len().alias("n_tt"),
                pl.col("subject").n_unique().alias("n_subjects"),
            ]
        ),
        streaming,
    ).rename({"_dt": "dt"})
    return out.with_columns(
        [
            pl.lit(float(alpha)).alias("alpha"),
            pl.lit(alpha_tag(alpha)).alias("alpha_tag"),
        ]
    ).select(
        ["alpha", "alpha_tag", "dt", "xstar_min", "xstar_max", "n_tt", "n_subjects"]
    ).sort("dt")


def compute_common_xstar_support(
    support: pl.DataFrame,
    alphas: Sequence[float],
    selected_dt: Sequence[int],
) -> Tuple[float, float]:
    expected = len(alphas) * len(selected_dt)
    if support.height != expected:
        present = {
            (round(float(r["alpha"]), 12), int(r["dt"]))
            for r in support.iter_rows(named=True)
        }
        missing = [
            (a, d)
            for a in alphas
            for d in selected_dt
            if (round(float(a), 12), int(d)) not in present
        ]
        raise ValueError("Some alpha x dt strata have no TT support: {}".format(missing))

    x_min = float(support["xstar_min"].max())
    x_max = float(support["xstar_max"].min())
    if not math.isfinite(x_min) or not math.isfinite(x_max) or x_max <= x_min:
        raise ValueError(
            "No common TT xstar support across alpha x dt: [{}, {}]".format(
                x_min, x_max
            )
        )
    return x_min, x_max


def build_tt_retention(
    class_table: pl.DataFrame,
    reference_alpha: float,
) -> pl.DataFrame:
    tt = (
        class_table.filter(pl.col("obs_class") == "TT")
        .group_by(["alpha", "alpha_tag"])
        .agg(pl.col("n_transitions").sum().alias("n_tt"))
        .sort("alpha", descending=True)
    )
    ref = tt.filter((pl.col("alpha") - float(reference_alpha)).abs() < 1e-12)
    if ref.height != 1:
        raise ValueError("Reference alpha TT count is not uniquely available.")
    nref = int(ref["n_tt"][0])
    return tt.with_columns(
        [
            pl.lit(nref).alias("reference_n_tt"),
            (pl.col("n_tt") / float(nref)).alias("fraction_of_reference_tt"),
        ]
    )


# =============================================================================
# Agreement at reference alpha
# =============================================================================

def transition_reference_agreement(
    classified: pl.LazyFrame,
    reference_alpha: float,
    streaming: bool,
) -> pl.DataFrame:
    names = schema_names(classified)
    if "obs_class_original" not in names:
        return pl.DataFrame(
            {
                "reference_alpha": [float(reference_alpha)],
                "note": ["obs_class_original unavailable"],
            }
        )
    return collect_frame(
        classified.group_by(["obs_class_original", "obs_class"])
        .agg(pl.len().alias("n"))
        .with_columns(pl.lit(float(reference_alpha)).alias("reference_alpha"))
        .sort(["obs_class_original", "obs_class"]),
        streaming,
    )


def trajectory_reference_agreement(
    classified: pl.LazyFrame,
    reference_alpha: float,
    streaming: bool,
) -> pl.DataFrame:
    names = schema_names(classified)
    if "observable_original" not in names:
        return pl.DataFrame(
            {
                "reference_alpha": [float(reference_alpha)],
                "note": ["observable_original unavailable"],
            }
        )

    original_text = (
        pl.col("observable_original")
        .cast(pl.String, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    original_bool = (
        pl.when(original_text.is_in(["true", "t", "1", "yes", "y", "observable", "detected", "detectable"]))
        .then(pl.lit(True))
        .when(original_text.is_in(["false", "f", "0", "no", "n", "non-observable", "undetected"]))
        .then(pl.lit(False))
        .otherwise(pl.lit(None, dtype=pl.Boolean))
        .alias("_original_bool")
    )

    return collect_frame(
        classified.with_columns(original_bool)
        .group_by(["_original_bool", "observable"])
        .agg(pl.len().alias("n"))
        .with_columns(pl.lit(float(reference_alpha)).alias("reference_alpha"))
        .sort(["_original_bool", "observable"]),
        streaming,
    )


# =============================================================================
# P09/P11/P13 command builders
# =============================================================================

def build_p09_command(
    args: argparse.Namespace,
    trajectory_table: Path,
    alpha: float,
    output_dir: Path,
    bin_edges_file: Optional[Path],
) -> List[str]:
    cmd = [
        str(args.python_executable),
        str(args.p09_script),
        "--trajectories", str(trajectory_table),
        "--dataset-label", "{}_alpha_{}".format(args.dataset_label, alpha_tag(alpha)),
        "--outdir", str(output_dir),
        "--dt", str(args.p09_dt),
        "--n-bins", str(args.p09_n_bins),
        "--binning", str(args.p09_binning),
        "--q-low", str(args.p09_q_low),
        "--q-high", str(args.p09_q_high),
        "--n-bootstrap", str(args.p09_n_bootstrap),
        "--bootstrap-seed", str(args.seed),
        "--parquet-compression", str(args.compression),
        "--parquet-row-group-size", str(args.p09_row_group_size),
    ]
    if bin_edges_file is not None:
        cmd.extend(["--bin-edges-file", str(bin_edges_file)])
    if args.restart_p09:
        cmd.append("--restart")
    return cmd


def build_p11_command(
    args: argparse.Namespace,
    transition_table: Path,
    alpha: float,
    output_dir: Path,
    x_min: float,
    x_max: float,
) -> List[str]:
    cmd = [
        str(args.python_executable),
        str(args.p11_script),
        "--transitions", str(transition_table),
        "--dataset-label", "{}_alpha_{}".format(args.dataset_label, alpha_tag(alpha)),
        "--outdir", str(output_dir),
        "--representation", "latent",
        "--conditioning", "xstar",
        "--mode", "TT",
        "--shared-range", "intersection",
        "--x-min", "{:.12g}".format(x_min),
        "--x-max", "{:.12g}".format(x_max),
        "--n-bins", str(args.p11_n_bins),
        "--min-n", str(args.p11_min_n),
        "--min-subjects", str(args.p11_min_subjects),
        "--n-bootstrap", str(args.p11_n_bootstrap),
        "--ci", str(args.ci),
        "--seed", str(args.seed),
    ]
    cmd.extend(["--dt"] + [str(v) for v in args.dt])
    if args.streaming:
        cmd.append("--streaming")
    return cmd


def build_p13_command(
    args: argparse.Namespace,
    p11_dir: Path,
    alpha: float,
    output_dir: Path,
) -> List[str]:
    cmd = [
        str(args.python_executable),
        str(args.p13_script),
        "--input-dir", str(p11_dir),
        "--outdir", str(output_dir),
        "--dataset-label", "{}_alpha_{}".format(args.dataset_label, alpha_tag(alpha)),
        "--required-mode", "TT",
        "--dt-values", ",".join(str(v) for v in args.dt),
        "--metrics", "var_dx,msd",
        "--estimators", "transition_weighted,equal_subject_weighted",
        "--bin-scope", str(args.p13_bin_scope),
        "--core-min-subject-fraction", str(args.p13_core_min_subject_fraction),
        "--min-dt-points", str(len(args.dt)),
        "--n-bootstrap", str(args.p13_n_bootstrap),
        "--seed", str(args.seed),
        "--ci", str(args.ci),
        "--compression", str(args.compression),
    ]
    if args.p13_core_min_subjects is not None:
        cmd.extend(["--core-min-subjects", str(args.p13_core_min_subjects)])
    return cmd


# =============================================================================
# Cross-alpha combinations
# =============================================================================

def combine_p11_curves(
    alpha_dirs: Dict[float, Path],
    out_dir: Path,
    compression: str,
) -> Optional[pl.DataFrame]:
    frames = []
    for alpha, folder in alpha_dirs.items():
        path = find_output_table(folder, "10_binned_dynamics_long")
        if path is None:
            continue
        frames.append(
            read_table(path).with_columns(
                [
                    pl.lit(float(alpha)).alias("alpha"),
                    pl.lit(alpha_tag(alpha)).alias("alpha_tag"),
                ]
            )
        )
    if not frames:
        return None
    combined = pl.concat(frames, how="diagonal_relaxed").sort(
        ["alpha", "dt", "bin_id"], descending=[True, False, False]
    )
    write_csv_parquet(combined, out_dir / "09_combined_p11_binned_dynamics", compression)
    # Compatibility alias for the existing threshold plotter.
    write_csv_parquet(combined, out_dir / "09_combined_step5_binned_dynamics", compression)
    return combined


def build_curve_robustness(
    combined: pl.DataFrame,
    reference_alpha: float,
    metrics: Sequence[str],
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    reference = combined.filter(
        (pl.col("alpha") - float(reference_alpha)).abs() < 1e-12
    )
    long_frames = []
    summary_rows = []
    for alpha in sorted(
        [float(v) for v in combined["alpha"].unique().to_list()],
        reverse=True,
    ):
        current = combined.filter((pl.col("alpha") - alpha).abs() < 1e-12)
        for metric in [m for m in metrics if m in combined.columns]:
            ref = reference.select(
                ["dt", "bin_id", pl.col(metric).cast(pl.Float64).alias("reference_value")]
            )
            cur = current.select(
                ["dt", "bin_id", pl.col(metric).cast(pl.Float64).alias("value")]
            )
            j = cur.join(ref, on=["dt", "bin_id"], how="inner").filter(
                pl.col("value").is_finite() & pl.col("reference_value").is_finite()
            )
            if j.height == 0:
                continue
            j = j.with_columns(
                [
                    pl.lit(alpha).alias("alpha"),
                    pl.lit(alpha_tag(alpha)).alias("alpha_tag"),
                    pl.lit(float(reference_alpha)).alias("reference_alpha"),
                    pl.lit(metric).alias("metric"),
                    (pl.col("value") - pl.col("reference_value")).alias("difference"),
                    (pl.col("value") - pl.col("reference_value")).abs().alias("absolute_difference"),
                    pl.when(pl.col("reference_value").abs() > 1e-12)
                    .then((pl.col("value") - pl.col("reference_value")) / pl.col("reference_value").abs())
                    .otherwise(None)
                    .alias("relative_difference"),
                ]
            )
            long_frames.append(
                j.select(
                    [
                        "alpha", "alpha_tag", "reference_alpha", "metric",
                        "dt", "bin_id", "value", "reference_value",
                        "difference", "absolute_difference", "relative_difference",
                    ]
                )
            )
            values = j["value"].to_numpy().astype(float)
            refs = j["reference_value"].to_numpy().astype(float)
            diff = values - refs
            if values.size >= 2 and np.std(values) > 0 and np.std(refs) > 0:
                r = float(np.corrcoef(values, refs)[0, 1])
            elif np.allclose(values, refs, equal_nan=True):
                r = 1.0
            else:
                r = np.nan
            summary_rows.append(
                {
                    "alpha": alpha,
                    "alpha_tag": alpha_tag(alpha),
                    "reference_alpha": float(reference_alpha),
                    "metric": metric,
                    "n_points": int(values.size),
                    "pearson_r": r,
                    "mae": float(np.mean(np.abs(diff))),
                    "rmse": float(np.sqrt(np.mean(diff ** 2))),
                    "median_abs_difference": float(np.median(np.abs(diff))),
                    "max_abs_difference": float(np.max(np.abs(diff))),
                }
            )
    long_table = pl.concat(long_frames, how="vertical") if long_frames else pl.DataFrame()
    summary = pl.DataFrame(summary_rows) if summary_rows else pl.DataFrame()
    return long_table, summary


def combine_p09_curves(
    alpha_dirs: Dict[float, Path],
    out_dir: Path,
    compression: str,
    reference_alpha: float,
) -> Optional[pl.DataFrame]:
    frames = []
    for alpha, folder in alpha_dirs.items():
        path = folder / "forward_drift_by_bin.csv"
        if not path.exists():
            continue
        df = pl.read_csv(path, infer_schema_length=10000)
        if "representation" in df.columns:
            df = df.filter(pl.col("representation") == "cross_combined")
        frames.append(
            df.with_columns(
                [
                    pl.lit(float(alpha)).alias("alpha"),
                    pl.lit(alpha_tag(alpha)).alias("alpha_tag"),
                ]
            )
        )
    if not frames:
        return None
    combined = pl.concat(frames, how="diagonal_relaxed").sort(
        ["alpha", "bin"], descending=[True, False]
    )
    write_csv_parquet(combined, out_dir / "15_combined_p09_forward_drift_by_bin", compression)

    # P09 has one primary lag, so temporarily map bin -> bin_id and add dt.
    temp = combined.with_columns(
        [
            pl.col("bin").cast(pl.Int64).alias("bin_id"),
            pl.lit(1).alias("dt"),
        ]
    )
    long, summary = build_curve_robustness(
        temp,
        reference_alpha,
        [m for m in P09_CURVE_METRICS if m in temp.columns],
    )
    if long.height:
        write_csv_parquet(long, out_dir / "16_p09_forward_robustness_vs_reference_long", compression)
    if summary.height:
        write_csv_parquet(summary, out_dir / "17_p09_forward_robustness_summary", compression)
    return combined


def combine_p13_outputs(
    alpha_dirs: Dict[float, Path],
    stem: str,
) -> Optional[pl.DataFrame]:
    frames = []
    for alpha, folder in alpha_dirs.items():
        path = find_output_table(folder, stem)
        if path is None:
            continue
        frames.append(
            read_table(path).with_columns(
                [
                    pl.lit(float(alpha)).alias("observation_alpha"),
                    pl.lit(alpha_tag(alpha)).alias("observation_alpha_tag"),
                    # Compatibility aliases used by the existing plotter.
                    pl.lit(float(alpha)).alias("detectability_alpha"),
                    pl.lit(alpha_tag(alpha)).alias("detectability_alpha_tag"),
                ]
            )
        )
    return pl.concat(frames, how="diagonal_relaxed") if frames else None


def resolve_cross_alpha_temporal_core(
    model_by_bin: pl.DataFrame,
    alphas: Sequence[float],
) -> pl.DataFrame:
    required = {
        "observation_alpha", "metric", "estimator",
        "bin_id", "bin_left", "bin_mid", "bin_right",
    }
    if not required.issubset(set(model_by_bin.columns)):
        raise ValueError(
            "Combined P13 table lacks common-core columns: {}".format(
                sorted(required - set(model_by_bin.columns))
            )
        )

    metrics = sorted(str(v) for v in model_by_bin["metric"].unique().to_list())
    estimators = sorted(str(v) for v in model_by_bin["estimator"].unique().to_list())
    strata: List[set] = []

    for alpha in alphas:
        for metric in metrics:
            for estimator in estimators:
                d = model_by_bin.filter(
                    (pl.col("observation_alpha") - float(alpha)).abs() < 1e-12
                ).filter(
                    (pl.col("metric") == metric)
                    & (pl.col("estimator") == estimator)
                )
                bins = set(int(v) for v in d["bin_id"].to_list())
                if not bins:
                    raise ValueError(
                        "No P13 bins for alpha={} metric={} estimator={}".format(
                            alpha, metric, estimator
                        )
                    )
                strata.append(bins)

    common = set.intersection(*strata)
    selected = largest_contiguous_run(common)
    if not selected:
        raise ValueError("No contiguous P13 cross-alpha common core exists.")

    geometry = (
        model_by_bin.filter(pl.col("bin_id").is_in(selected))
        .select(["bin_id", "bin_left", "bin_mid", "bin_right"])
        .unique(subset=["bin_id"])
        .sort("bin_id")
    )
    return geometry.with_columns(
        [
            pl.lit(",".join("{:g}".format(v) for v in alphas)).alias("alpha_values"),
            pl.lit(",".join(metrics)).alias("metrics"),
            pl.lit(",".join(estimators)).alias("estimators"),
            pl.lit(True).alias("selected_for_common_core"),
        ]
    )


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    args = build_argparser().parse_args()
    alphas = validate_args(args)

    args.trajectories = args.trajectories.expanduser()
    args.transitions = args.transitions.expanduser()
    args.out_dir = args.out_dir.expanduser()
    args.p09_script = args.p09_script.expanduser()
    args.p11_script = args.p11_script.expanduser()
    args.p13_script = args.p13_script.expanduser()

    if not args.trajectories.exists():
        raise FileNotFoundError(args.trajectories)
    if not args.transitions.exists():
        raise FileNotFoundError(args.transitions)
    if args.run_p09 and not args.p09_script.exists():
        raise FileNotFoundError(args.p09_script)
    if args.run_p11 and not args.p11_script.exists():
        raise FileNotFoundError(args.p11_script)
    if args.run_p13 and not args.p13_script.exists():
        raise FileNotFoundError(args.p13_script)

    out_dir = ensure_dir(args.out_dir)
    threshold_dir = ensure_dir(out_dir / "threshold_tables")
    p09_root = ensure_dir(out_dir / "p09_by_alpha")
    p11_root = ensure_dir(out_dir / "p11_by_alpha")
    p13_root = ensure_dir(out_dir / "p13_by_alpha")
    logs_dir = ensure_dir(out_dir / "logs")

    # Step-2 endpoint empirical p-values.
    source_lf, source_label, pvalue_col, source_columns = resolve_pvalue_source(args)
    endpoint_p = endpoint_pvalues(source_lf, pvalue_col, args)

    # Enrich both P04 and P05 once.
    traj_master_path = threshold_dir / "00_endpoint_pvalue_enriched_trajectories.parquet"
    trans_master_path = threshold_dir / "01_endpoint_pvalue_enriched_transitions.parquet"

    traj_master, traj_audit = build_enriched_trajectories(
        args.trajectories,
        endpoint_p,
        traj_master_path,
        args.compression,
        args.streaming,
    )
    trans_master, trans_audit = build_enriched_transitions(
        args.transitions,
        endpoint_p,
        trans_master_path,
        args.compression,
        args.streaming,
    )
    write_csv(traj_audit, out_dir / "01_trajectory_pvalue_join_audit.csv")
    write_csv(trans_audit, out_dir / "02_transition_pvalue_join_audit.csv")

    n_missing_traj = int(traj_audit["n_missing_p"][0])
    n_missing_trans = int(trans_audit["n_missing_either_endpoint"][0])
    if (n_missing_traj > 0 or n_missing_trans > 0) and not args.allow_missing_pvalues:
        raise ValueError(
            "Missing empirical endpoint p-values detected "
            "(trajectory states={}, transitions with missing endpoint={}). "
            "Inspect join audits. Use --allow-missing-pvalues only if treating "
            "missing endpoints as F is explicitly intended.".format(
                n_missing_traj, n_missing_trans
            )
        )

    class_tables: List[pl.DataFrame] = []
    support_tables: List[pl.DataFrame] = []
    trajectory_paths: Dict[float, Path] = {}
    transition_paths: Dict[float, Path] = {}
    transition_agreement = None
    trajectory_agreement = None

    for alpha in alphas:
        tag = alpha_tag(alpha)

        c_traj = trajectory_classified(traj_master, alpha)
        traj_path = threshold_dir / "trajectories_alpha_{}.parquet".format(tag)
        sink_parquet(c_traj, traj_path, args.compression, args.streaming)
        trajectory_paths[alpha] = traj_path

        c_trans = transition_classified(trans_master, alpha)
        trans_path = threshold_dir / "transitions_alpha_{}.parquet".format(tag)
        sink_parquet(c_trans, trans_path, args.compression, args.streaming)
        transition_paths[alpha] = trans_path

        class_tables.append(
            class_composition(c_trans, alpha, args.dt, args.streaming)
        )
        support_tables.append(
            tt_support(c_trans, alpha, args.dt, args.streaming)
        )

        if abs(alpha - float(args.reference_alpha)) < 1e-12:
            transition_agreement = transition_reference_agreement(
                c_trans, args.reference_alpha, args.streaming
            )
            trajectory_agreement = trajectory_reference_agreement(
                c_traj, args.reference_alpha, args.streaming
            )

        print("[WRITE] alpha={} trajectories={}".format(alpha, traj_path.name))
        print("[WRITE] alpha={} transitions={}".format(alpha, trans_path.name))

    class_table = pl.concat(class_tables, how="vertical").sort(
        ["alpha", "dt", "obs_class"], descending=[True, False, False]
    )
    support_table = pl.concat(support_tables, how="vertical").sort(
        ["alpha", "dt"], descending=[True, False]
    )

    write_csv(class_table, out_dir / "03_class_composition_by_alpha_dt.csv")
    write_csv(support_table, out_dir / "04_tt_support_by_alpha_dt.csv")

    if transition_agreement is not None:
        write_csv(
            transition_agreement,
            out_dir / "02_reference_alpha_classification_agreement.csv",
        )
    if trajectory_agreement is not None:
        write_csv(
            trajectory_agreement,
            out_dir / "02b_reference_alpha_trajectory_agreement.csv",
        )

    x_min, x_max = compute_common_xstar_support(
        support_table, alphas, args.dt
    )
    common_support = pl.DataFrame(
        [
            {
                "xstar_common_min": x_min,
                "xstar_common_max": x_max,
                "width": x_max - x_min,
                "n_alphas": len(alphas),
                "n_dt": len(args.dt),
                "n_alpha_dt_strata": len(alphas) * len(args.dt),
                "n_bins": args.p11_n_bins,
                "alpha_values": ",".join("{:g}".format(v) for v in alphas),
                "reference_alpha": float(args.reference_alpha),
            }
        ]
    )
    write_csv(common_support, out_dir / "05_common_xstar_support.csv")
    write_csv(
        build_tt_retention(class_table, args.reference_alpha),
        out_dir / "06_tt_retention_by_alpha.csv",
    )
    print("[SUPPORT] common P11 xstar=[{:.6f}, {:.6f}]".format(x_min, x_max))

    # ------------------------------------------------------------------
    # Build/run P09. Reference alpha first so it defines common x0 bins.
    # ------------------------------------------------------------------
    p09_dirs: Dict[float, Path] = {}
    p09_commands: List[Dict[str, object]] = []

    p09_order = [float(args.reference_alpha)] + [
        a for a in alphas if abs(a - float(args.reference_alpha)) >= 1e-12
    ]
    common_p09_edges: Optional[Path] = (
        args.p09_bin_edges_file.expanduser()
        if args.p09_bin_edges_file is not None
        else None
    )

    for alpha in p09_order:
        tag = alpha_tag(alpha)
        folder = ensure_dir(p09_root / "alpha_{}".format(tag))
        p09_dirs[alpha] = folder

        edge_arg = common_p09_edges
        cmd = build_p09_command(
            args,
            trajectory_paths[alpha],
            alpha,
            folder,
            edge_arg,
        )
        p09_commands.append(
            {
                "alpha": alpha,
                "alpha_tag": tag,
                "output_dir": str(folder),
                "bin_edges_input": str(edge_arg) if edge_arg is not None else "",
                "command": quote_command(cmd),
            }
        )
        if args.run_p09:
            run_command(cmd, logs_dir / "p09_alpha_{}.log".format(tag))

            if (
                common_p09_edges is None
                and abs(alpha - float(args.reference_alpha)) < 1e-12
            ):
                generated = folder / "forward_drift_bin_edges.csv"
                if not generated.exists():
                    raise RuntimeError(
                        "Reference P09 completed but forward_drift_bin_edges.csv is missing."
                    )
                common_p09_edges = out_dir / "14_p09_common_forward_bin_edges.csv"
                common_p09_edges.write_bytes(generated.read_bytes())
                print("[P09 BINS] reference alpha defined {}".format(common_p09_edges))

    write_csv(pl.DataFrame(p09_commands), out_dir / "20_p09_commands.csv")

    # ------------------------------------------------------------------
    # Build/run P11 and then P13.
    # ------------------------------------------------------------------
    p11_dirs: Dict[float, Path] = {}
    p13_dirs: Dict[float, Path] = {}
    p11_commands: List[Dict[str, object]] = []
    p13_commands: List[Dict[str, object]] = []

    for alpha in alphas:
        tag = alpha_tag(alpha)

        p11_dir = ensure_dir(p11_root / "alpha_{}".format(tag))
        p11_dirs[alpha] = p11_dir
        cmd11 = build_p11_command(
            args,
            transition_paths[alpha],
            alpha,
            p11_dir,
            x_min,
            x_max,
        )
        p11_commands.append(
            {
                "alpha": alpha,
                "alpha_tag": tag,
                "output_dir": str(p11_dir),
                "command": quote_command(cmd11),
            }
        )
        if args.run_p11:
            run_command(cmd11, logs_dir / "p11_alpha_{}.log".format(tag))

        p13_dir = ensure_dir(p13_root / "alpha_{}".format(tag))
        p13_dirs[alpha] = p13_dir
        cmd13 = build_p13_command(args, p11_dir, alpha, p13_dir)
        p13_commands.append(
            {
                "alpha": alpha,
                "alpha_tag": tag,
                "output_dir": str(p13_dir),
                "command": quote_command(cmd13),
            }
        )
        if args.run_p13:
            run_command(cmd13, logs_dir / "p13_alpha_{}.log".format(tag))

    write_csv(pl.DataFrame(p11_commands), out_dir / "21_p11_commands.csv")
    write_csv(pl.DataFrame(p13_commands), out_dir / "22_p13_commands.csv")

    # ------------------------------------------------------------------
    # Combine P09 results.
    # ------------------------------------------------------------------
    if args.run_p09:
        combine_p09_curves(
            p09_dirs,
            out_dir,
            args.compression,
            args.reference_alpha,
        )

    # ------------------------------------------------------------------
    # Combine P11 and write compatibility tables for existing plotter.
    # ------------------------------------------------------------------
    if args.run_p11:
        combined11 = combine_p11_curves(
            p11_dirs,
            out_dir,
            args.compression,
        )
        if combined11 is None:
            raise RuntimeError("P11 completed but no 10_binned_dynamics_long outputs were found.")
        robust_long, robust_summary = build_curve_robustness(
            combined11,
            args.reference_alpha,
            P11_CURVE_METRICS,
        )
        write_csv_parquet(
            robust_long,
            out_dir / "10_curve_robustness_vs_reference_long",
            args.compression,
        )
        write_csv_parquet(
            robust_summary,
            out_dir / "11_curve_robustness_summary",
            args.compression,
        )

    # ------------------------------------------------------------------
    # Combine P13, preserve old plotter aliases, resolve common temporal core.
    # ------------------------------------------------------------------
    cross_alpha_core = None
    if args.run_p13:
        model_by_bin = combine_p13_outputs(
            p13_dirs, "05_model_comparison_by_bin"
        )
        model_summary = combine_p13_outputs(
            p13_dirs, "06_model_comparison_summary"
        )
        if model_by_bin is None:
            raise RuntimeError("P13 completed but no 05_model_comparison_by_bin outputs were found.")

        write_csv_parquet(
            model_by_bin,
            out_dir / "12_combined_p13_model_comparison_by_bin",
            args.compression,
        )
        # Compatibility with the existing threshold-robustness plotter.
        write_csv_parquet(
            model_by_bin,
            out_dir / "12_combined_step51_model_comparison_by_bin",
            args.compression,
        )

        if model_summary is not None:
            write_csv_parquet(
                model_summary,
                out_dir / "13_combined_p13_model_summary",
                args.compression,
            )
            write_csv_parquet(
                model_summary,
                out_dir / "13_combined_step51_model_summary",
                args.compression,
            )

        cross_alpha_core = resolve_cross_alpha_temporal_core(
            model_by_bin, alphas
        )
        write_csv(
            cross_alpha_core,
            out_dir / "14_cross_alpha_common_core.csv",
        )
        print(
            "[P13 CORE] bins={}".format(
                [int(v) for v in cross_alpha_core["bin_id"].to_list()]
            )
        )

    # ------------------------------------------------------------------
    # Run metadata.
    # ------------------------------------------------------------------
    metadata = {
        "script": Path(__file__).name,
        "version": VERSION,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "dataset_label": args.dataset_label,
        "scientific_scope": "operational_observation_state_threshold_sensitivity",
        "latent_states_refitted": False,
        "posterior_predictive_detectability_thresholded": False,
        "classification_rule": "T(alpha)=1[p_endpoint<alpha]",
        "alphas": alphas,
        "reference_alpha": float(args.reference_alpha),
        "selected_dt": list(args.dt),
        "pvalue_source": source_label,
        "pvalue_column": pvalue_col,
        "source_columns": source_columns,
        "trajectories": str(args.trajectories),
        "transitions": str(args.transitions),
        "common_p11_xstar_min": x_min,
        "common_p11_xstar_max": x_max,
        "p09_common_bin_edges": str(common_p09_edges) if common_p09_edges is not None else None,
        "canonical_estimators_reused": {
            "P09": str(args.p09_script),
            "P11": str(args.p11_script),
            "P13": str(args.p13_script),
        },
        "run_p09": bool(args.run_p09),
        "run_p11": bool(args.run_p11),
        "run_p13": bool(args.run_p13),
        "bootstrap": {
            "P09": int(args.p09_n_bootstrap),
            "P11": int(args.p11_n_bootstrap),
            "P13": int(args.p13_n_bootstrap),
            "seed": int(args.seed),
        },
        "threshold_trajectory_tables": {
            alpha_tag(a): str(trajectory_paths[a]) for a in alphas
        },
        "threshold_transition_tables": {
            alpha_tag(a): str(transition_paths[a]) for a in alphas
        },
        "cross_alpha_temporal_core_bins": (
            [int(v) for v in cross_alpha_core["bin_id"].to_list()]
            if cross_alpha_core is not None else None
        ),
    }
    (out_dir / "00_run_config.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    debug_lines = [
        "OPERATIONAL OBSERVATION-THRESHOLD ROBUSTNESS",
        "==========================================",
        "",
        "Version: {}".format(VERSION),
        "Alphas: {}".format(alphas),
        "Reference alpha: {}".format(args.reference_alpha),
        "Rule: T(alpha) = 1[p_endpoint < alpha]",
        "Latent states refitted: NO",
        "Posterior-predictive detectability thresholded: NO",
        "",
        "P09 reused: {}".format(args.p09_script),
        "P11 reused: {}".format(args.p11_script),
        "P13 reused: {}".format(args.p13_script),
        "",
        "P11 common xstar support: [{:.12g}, {:.12g}]".format(x_min, x_max),
        "P09 common bin edges: {}".format(common_p09_edges),
        "P09 run: {}".format(args.run_p09),
        "P11 run: {}".format(args.run_p11),
        "P13 run: {}".format(args.run_p13),
        "",
        "Interpretation: sensitivity to the operational observation-state "
        "classification threshold, not biological presence/absence and not "
        "posterior-predictive detectability.",
    ]
    if cross_alpha_core is not None:
        debug_lines.append(
            "P13 cross-alpha common core bins: {}".format(
                [int(v) for v in cross_alpha_core["bin_id"].to_list()]
            )
        )
    (out_dir / "DEBUG_REPORT.txt").write_text(
        "\n".join(debug_lines) + "\n", encoding="utf-8"
    )

    print("[DONE] {}".format(out_dir))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise
