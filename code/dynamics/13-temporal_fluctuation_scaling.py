#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13-temporal_fluctuation_scaling.py
==================================

ClonoDynamics Step 13: finite-lag temporal scaling of genuine longitudinal
clonotype fluctuations downstream of Step 11.

Scientific role
---------------
The manuscript-primary input is the current Step-11 analysis:

    representation = latent
    conditioning   = xstar
    mode           = TT
    dt             = 1,2,3,4,5

Step 13 tests whether abundance-resolved longitudinal fluctuation magnitude
shows systematic accumulation with temporal lag. The primary metric is:

    Var(dx_latent | xstar_latent, dt)

and MSD is retained as a complementary second-moment measure because genuine
forward drift need not be zero.

The analysis is restricted to an internally comparable absolute-xstar abundance
core. Under the default strict core rule, each retained abundance bin must
contain every subject structurally available at each requested lag; the largest
contiguous qualifying abundance interval is used.

Six finite-lag models are compared in every retained abundance bin, metric and
estimator:

    constant
        M(dt) = K

    signed_linear
        M(dt) = K + D * (dt - 1), D signed

    incremental_linear
        M(dt) = K + D * (dt - 1), D >= 0

    linear
        M(dt) = D * dt

    power_law
        M(dt) = C * dt**alpha

    saturating
        M(dt) = A * (1 - exp(-dt / tau))

The constant model is empirical over the observed lag range only. The
zero-origin linear model is a Brownian comparator. A saturating second-moment
fit is not interpreted as a direct Ornstein-Uhlenbeck test.

Joint subject bootstrap
-----------------------
Subjects are sampled once with replacement per bootstrap replicate and the same
multiplicities are propagated jointly across all lags, abundance bins, metrics
and estimators. This preserves subject-level covariance across dt.

Interpretation combines:

    1. signed-slope evidence from the joint subject-bootstrap distribution;
    2. descriptive/AICc support for the finite-lag model set.

With only five lags, AICc preference is not treated as a stand-alone
mechanistic winner criterion.

Input
-----
--input-dir points to the Step-11 output directory and uses:

    00_analysis_metadata.json
    08_subject_availability_by_dt.csv/parquet
    09_subject_bin_dynamics.csv/parquet
    10_binned_dynamics_long.csv/parquet

The Step-5 transition table is not reloaded.

Outputs
-------
    00_run_config.json
    01_input_manifest.csv/parquet
    02_tested_abundance_bins.csv/parquet
    03_observed_six_model_fits_long.csv/parquet
    04_joint_subject_bootstrap_six_model_comparison.csv/parquet
    05_model_comparison_by_bin.csv/parquet
    06_model_comparison_summary.csv/parquet
    07_interpretation_aid.csv/parquet
    08_subject_availability_and_requirements_by_dt.csv/parquet
    DEBUG_REPORT.txt

Pipeline boundary
-----------------
Step 12 does not feed Step 13 computationally and the pseudo technical baseline
is not subtracted from Step-11 curves before temporal-model fitting.

Step 14 evaluates interval/calendar-position and subject-composition structure.
Step 15 evaluates sensitivity of the Step-13 conclusion to inclusion of
detectability-boundary-crossing transitions.

Preferred interpretation
------------------------
When positive accumulation is not supported, use wording such as:

    "no systematic positive accumulation over the observed temporal lags"

rather than claiming complete temporal independence.

Recommended primary run
-----------------------
python3 ./code/analysis/13-temporal_fluctuation_scaling.py \
    --input-dir ./dataset_longitudinal_clonodynamics_results/longitudinal/11-longitudinal_fluctuation_dynamics \
    --outdir ./dataset_longitudinal_clonodynamics_results/longitudinal/13-temporal-fluctuation-scaling \
    --dt-values 1,2,3,4,5 \
    --metrics var_dx,msd \
    --estimators transition_weighted,equal_subject_weighted \
    --bin-scope core \
    --core-min-subject-fraction 1.0 \
    --min-dt-points 5 \
    --n-bootstrap 2000 \
    --seed 123 \
    --ci 95

Python >= 3.9
Required: numpy, pandas
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v2-step13-temporal-fluctuation-scaling-2026-08-25"
PRIMARY_REPRESENTATION = "latent"
PRIMARY_CONDITIONING = "xstar"
PRIMARY_MODE = "TT"
PRIMARY_METRIC = "var_dx"
COMPLEMENTARY_METRIC = "msd"
MODELS = ("constant", "signed_linear", "incremental_linear", "linear", "power_law", "saturating")
SUPPORTED_METRICS = ("var_dx", "msd")
SUPPORTED_ESTIMATORS = ("transition_weighted", "equal_subject_weighted")


# =============================================================================
# CLI and I/O
# =============================================================================


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 13: analyze finite-lag temporal scaling of genuine longitudinal "
            "Step-11 fluctuations using a joint subject-level bootstrap."
        ),
    )
    parser.add_argument(
        "--input-dir",
        dest="step11_dir",
        required=True,
        help=("Current Step-11 output directory from "
              "11-longitudinal_fluctuation_dynamics.py."),
    )
    parser.add_argument("--outdir", "--out-dir", dest="out_dir", required=True)
    parser.add_argument(
        "--dataset-label",
        default="",
        help="Dataset label. Empty uses Step-11 metadata or the input-directory name.",
    )
    parser.add_argument(
        "--required-mode",
        choices=["TT", "TF", "FT", "FF", "NON_TT", "NON_FF", "ALL"],
        default=PRIMARY_MODE,
        help=(
            "Required upstream Step-11 transition mode. The manuscript-primary run "
            "uses TT; other modes require --allow-nonprimary."
        ),
    )
    parser.add_argument(
        "--allow-nonprimary",
        action="store_true",
        help="Explicitly permit a non-TT upstream sensitivity analysis.",
    )
    parser.add_argument(
        "--subject-stats",
        default="",
        help="Optional explicit 09_subject_bin_dynamics table.",
    )
    parser.add_argument(
        "--curves",
        default="",
        help="Optional explicit 10_binned_dynamics_long table.",
    )
    parser.add_argument(
        "--subject-availability",
        default="",
        help=(
            "Optional explicit 08_subject_availability_by_dt table from Step 11. "
            "If omitted, the script uses that table when present and otherwise "
            "derives per-dt availability from 09_subject_bin_dynamics."
        ),
    )
    parser.add_argument(
        "--dt-values",
        default="1,2,3,4,5",
        help="Comma-separated temporal lags fitted jointly.",
    )
    parser.add_argument(
        "--metrics",
        default="var_dx,msd",
        help="Metrics to fit. Primary: var_dx; complementary: msd.",
    )
    parser.add_argument(
        "--estimators",
        default="transition_weighted,equal_subject_weighted",
        help=(
            "Estimators to summarize. transition_weighted is primary; "
            "equal_subject_weighted is a sensitivity analysis."
        ),
    )
    parser.add_argument(
        "--bin-scope",
        choices=["core", "all_complete"],
        default="core",
        help=(
            "core retains the largest contiguous run satisfying the per-dt subject "
            "requirement at every selected lag; all_complete retains all qualifying bins."
        ),
    )
    parser.add_argument(
        "--bin-ids",
        default="",
        help="Optional explicit comma-separated Step-11 bin IDs; overrides --bin-scope.",
    )
    parser.add_argument(
        "--core-min-subjects",
        type=int,
        default=None,
        help=(
            "Optional absolute subject requirement per dt, capped at the number "
            "structurally available at that lag."
        ),
    )
    parser.add_argument(
        "--core-min-subject-fraction",
        type=float,
        default=1.0,
        help=(
            "Fraction of structurally available subjects required in each bin at "
            "each dt. Default 1.0 defines the strict longitudinal core."
        ),
    )
    parser.add_argument("--min-dt-points", type=int, default=5)
    parser.add_argument("--n-bootstrap", "--n-boot", dest="n_boot", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--ci", type=float, default=95.0)

    parser.add_argument("--alpha-min", type=float, default=0.01)
    parser.add_argument("--alpha-max", type=float, default=2.0)
    parser.add_argument("--tau-min", type=float, default=0.001)
    parser.add_argument("--tau-max", type=float, default=100.0)

    parser.add_argument(
        "--model-preference-threshold",
        type=float,
        default=0.80,
        help=(
            "Operational bootstrap support threshold used only in the interpretation aid; "
            "it is not a mechanistic winner criterion."
        ),
    )
    parser.add_argument(
        "--signed-slope-probability-threshold",
        type=float,
        default=0.95,
        help="Bootstrap sign-probability threshold for positive/negative signed slope evidence.",
    )
    parser.add_argument(
        "--subdiffusion-probability-threshold",
        type=float,
        default=0.95,
        help="Operational bootstrap fraction alpha<1 used for a subdiffusion descriptor.",
    )
    parser.add_argument(
        "--delta-aicc-threshold",
        type=float,
        default=2.0,
        help="Operational pairwise AICc separation threshold.",
    )
    parser.add_argument(
        "--parameter-bound-warning-fraction",
        type=float,
        default=0.25,
        help=(
            "Maximum tolerated bootstrap fraction at a parameter bound when "
            "describing nonlinear parameters as identifiable."
        ),
    )

    parser.add_argument("--write-csv", action="store_true", default=True)
    parser.add_argument("--no-write-csv", action="store_false", dest="write_csv")
    parser.add_argument("--write-parquet", action="store_true", default=True)
    parser.add_argument("--no-write-parquet", action="store_false", dest="write_parquet")
    parser.add_argument("--compression", default="zstd")
    return parser

def parse_csv_strings(text: str) -> List[str]:
    return [value.strip() for value in str(text or "").split(",") if value.strip()]


def parse_csv_ints(text: str) -> List[int]:
    values = []
    for raw in parse_csv_strings(text):
        value = float(raw)
        if not value.is_integer():
            raise ValueError("Expected integer value, got {!r}".format(raw))
        values.append(int(value))
    return values


def find_table(folder: Path, stem: str) -> Path:
    for suffix in (".parquet", ".pq", ".csv", ".tsv", ".txt"):
        candidate = folder / (stem + suffix)
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find '{}' in {}".format(stem, folder))


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".parquet", ".pq"}:
        try:
            return pd.read_parquet(path)
        except Exception as pandas_error:
            try:
                import polars as pl  # type: ignore

                return pl.read_parquet(path).to_pandas()
            except Exception as polars_error:
                # Step 11 normally writes CSV and parquet copies together.
                # Parquet support is optional for Step 13, so if no parquet
                # engine is available, transparently use the canonical
                # same-stem CSV when it exists.
                csv_fallback = path.with_suffix(".csv")
                if csv_fallback.exists():
                    warnings.warn(
                        "Could not read parquet {} (pandas: {}; Polars: {}). "
                        "Falling back to {}.".format(
                            path, pandas_error, polars_error, csv_fallback
                        )
                    )
                    return pd.read_csv(csv_fallback)

                raise RuntimeError(
                    "Could not read parquet {}. pandas error: {}; Polars error: {}".format(
                        path, pandas_error, polars_error
                    )
                )
    raise ValueError("Unsupported table format: {}".format(path))


def write_table(
    frame: pd.DataFrame,
    stem: Path,
    write_csv: bool,
    write_parquet: bool,
    compression: str,
) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    if write_csv:
        frame.to_csv(stem.with_suffix(".csv"), index=False)
    if write_parquet:
        try:
            frame.to_parquet(
                stem.with_suffix(".parquet"), index=False, compression=compression
            )
        except Exception as exc:
            warnings.warn(
                "Could not write parquet {}: {}. CSV output is unaffected.".format(
                    stem.name, exc
                )
            )


def load_metadata(step11_dir: Path) -> Dict[str, object]:
    path = step11_dir / "00_analysis_metadata.json"
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_upstream_metadata(
    metadata: Dict[str, object],
    selected_dts: Sequence[int],
    required_mode: str,
) -> None:
    """
    Fail fast unless Step 11 used latent relative-frequency displacement,
    xstar conditioning, the explicitly required transition mode, and the
    requested temporal lags. The default required mode is TT.
    """
    if not metadata:
        raise ValueError(
            "Step-11 metadata are required for temporal scaling. "
            "00_analysis_metadata.json was not found or is empty."
        )

    expected = {
        "representation": "latent",
        "conditioning": "xstar",
        "mode": str(required_mode).upper(),
    }
    mismatches = []
    for key, wanted in expected.items():
        observed = str(metadata.get(key, "")).strip()
        if observed != wanted:
            mismatches.append("{}={!r} (expected {!r})".format(key, observed, wanted))

    upstream_dts = metadata.get("selected_dt", [])
    try:
        upstream_dt_set = set(int(v) for v in upstream_dts)
    except Exception:
        upstream_dt_set = set()
    missing = sorted(set(int(v) for v in selected_dts) - upstream_dt_set)
    if missing:
        mismatches.append(
            "requested dt values absent from Step-11 metadata: {}".format(missing)
        )

    resolved = metadata.get("resolved_columns", {})
    if isinstance(resolved, dict):
        dx_source = str(resolved.get("dx", ""))
        x_source = str(resolved.get("x_condition", ""))
        if dx_source not in {"dx_latent", "__computed_x1_minus_x0"}:
            mismatches.append(
                "resolved dx source={!r} (expected dx_latent)".format(dx_source)
            )
        if x_source != "xstar_latent":
            mismatches.append(
                "resolved conditioning source={!r} (expected xstar_latent)".format(
                    x_source
                )
            )

    if mismatches:
        raise ValueError(
            "Upstream Step-11 representation does not match the frozen temporal-scaling "
            "analysis:\n  - " + "\n  - ".join(mismatches)
        )


# =============================================================================
# Validation and tested-bin selection
# =============================================================================


def require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(
            "{} is missing required columns {}. Available columns: {}".format(
                label, missing, list(frame.columns)
            )
        )


def normalize_inputs(
    subject_stats: pd.DataFrame,
    curves: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    require_columns(
        subject_stats,
        [
            "dt",
            "bin_id",
            "subject",
            "n_subject_bin",
            "sum_dx",
            "sum_dx2",
            "subject_mean_dx",
            "subject_msd",
        ],
        "09_subject_bin_dynamics",
    )
    require_columns(
        curves,
        ["dt", "bin_id", "bin_left", "bin_mid", "bin_right", "n_subjects"],
        "10_binned_dynamics_long",
    )

    subject_stats = subject_stats.copy()
    curves = curves.copy()
    for column in [
        "dt",
        "bin_id",
        "n_subject_bin",
        "sum_dx",
        "sum_dx2",
        "subject_mean_dx",
        "subject_msd",
    ]:
        subject_stats[column] = pd.to_numeric(subject_stats[column], errors="coerce")
    subject_stats["subject"] = subject_stats["subject"].astype(str)
    subject_stats = subject_stats.replace([np.inf, -np.inf], np.nan)
    subject_stats = subject_stats.dropna(
        subset=["dt", "bin_id", "subject", "n_subject_bin", "sum_dx", "sum_dx2"]
    )
    subject_stats["dt"] = subject_stats["dt"].astype(int)
    subject_stats["bin_id"] = subject_stats["bin_id"].astype(int)

    # Collapse accidental duplicate subject x dt x bin rows using sufficient statistics.
    grouped = (
        subject_stats.groupby(["dt", "bin_id", "subject"], as_index=False)
        .agg(
            n_subject_bin=("n_subject_bin", "sum"),
            sum_dx=("sum_dx", "sum"),
            sum_dx2=("sum_dx2", "sum"),
        )
        .reset_index(drop=True)
    )
    grouped["subject_mean_dx"] = grouped["sum_dx"] / grouped["n_subject_bin"]
    grouped["subject_msd"] = grouped["sum_dx2"] / grouped["n_subject_bin"]

    for column in ["dt", "bin_id", "bin_left", "bin_mid", "bin_right", "n_subjects"]:
        curves[column] = pd.to_numeric(curves[column], errors="coerce")
    curves = curves.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["dt", "bin_id", "bin_left", "bin_mid", "bin_right", "n_subjects"]
    )
    curves["dt"] = curves["dt"].astype(int)
    curves["bin_id"] = curves["bin_id"].astype(int)
    return grouped, curves


def largest_contiguous_run(bin_ids: Iterable[int]) -> List[int]:
    values = sorted(set(int(value) for value in bin_ids))
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


def normalize_subject_availability(
    availability: pd.DataFrame,
    selected_dts: Sequence[int],
) -> pd.DataFrame:
    """Normalize a Step-11 subject-availability table to one row per dt."""
    require_columns(availability, ["dt", "subject"], "08_subject_availability_by_dt")
    d = availability.copy()
    d["dt"] = pd.to_numeric(d["dt"], errors="coerce")
    d = d.dropna(subset=["dt", "subject"]).copy()
    d["dt"] = d["dt"].astype(int)
    d["subject"] = d["subject"].astype(str)
    d = d[d["dt"].isin([int(v) for v in selected_dts])].copy()
    summary = (
        d.groupby("dt", as_index=False)
        .agg(n_subjects_available=("subject", "nunique"))
        .sort_values("dt")
        .reset_index(drop=True)
    )
    missing = sorted(set(int(v) for v in selected_dts) - set(summary["dt"].astype(int)))
    if missing:
        raise ValueError(
            "08_subject_availability_by_dt is missing requested dt values: {}".format(missing)
        )
    return summary


def build_subject_requirements_by_dt(
    subject_stats: pd.DataFrame,
    selected_dts: Sequence[int],
    absolute_min_subjects: Optional[int],
    subject_fraction: float,
    availability_summary: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Define the subject requirement relative to structural availability at each dt.

    Preferred source of structural availability is Step 11
    08_subject_availability_by_dt. If unavailable, availability is derived from
    09_subject_bin_dynamics. With --core-min-subject-fraction 1.0, a qualifying
    bin must therefore contain every subject structurally available at that lag,
    not every subject in the full cohort.
    """
    if availability_summary is not None:
        avail_map = dict(
            zip(
                availability_summary["dt"].astype(int),
                availability_summary["n_subjects_available"].astype(int),
            )
        )
        availability_source = "08_subject_availability_by_dt"
    else:
        avail_map = {}
        for dt in selected_dts:
            d = subject_stats[subject_stats["dt"].eq(int(dt))]
            avail_map[int(dt)] = int(d["subject"].astype(str).nunique())
        availability_source = "09_subject_bin_dynamics_fallback"

    rows = []
    for dt in selected_dts:
        d = subject_stats[subject_stats["dt"].eq(int(dt))]
        n_available = int(avail_map.get(int(dt), 0))
        if n_available < 2:
            raise ValueError(
                "dt={} has fewer than two structurally available subjects.".format(dt)
            )

        if absolute_min_subjects is None:
            required = int(math.ceil(float(subject_fraction) * n_available))
        else:
            required = int(absolute_min_subjects)

        required = max(2, min(n_available, required))
        rows.append(
            {
                "dt": int(dt),
                "n_subjects_available": n_available,
                "n_subjects_required": required,
                "required_fraction_of_available": required / float(n_available),
                "availability_source": availability_source,
                "n_subject_bin_rows": int(len(d)),
                "n_transitions_in_subject_stats": int(
                    pd.to_numeric(d["n_subject_bin"], errors="coerce").fillna(0).sum()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("dt").reset_index(drop=True)

def select_tested_bins(
    curves: pd.DataFrame,
    selected_dts: Sequence[int],
    subject_requirements: pd.DataFrame,
    scope: str,
    explicit_bin_ids: Sequence[int],
) -> Tuple[List[int], pd.DataFrame]:
    """Select bins with complete lag coverage and availability-aware subject coverage."""
    d = curves[curves["dt"].isin(selected_dts)].copy()
    req_map = dict(
        zip(
            subject_requirements["dt"].astype(int),
            subject_requirements["n_subjects_required"].astype(int),
        )
    )
    avail_map = dict(
        zip(
            subject_requirements["dt"].astype(int),
            subject_requirements["n_subjects_available"].astype(int),
        )
    )

    rows = []
    for bin_id, group in d.groupby("bin_id", sort=True):
        present_dts = set(group["dt"].astype(int))
        present_all = all(int(dt) in present_dts for dt in selected_dts)
        dt_ok = {}
        availability_fractions = []

        first = group.iloc[0]
        row: Dict[str, object] = {
            "bin_id": int(bin_id),
            "n_dt_present": int(group["dt"].nunique()),
            "bin_left": float(first["bin_left"]),
            "bin_mid": float(first["bin_mid"]),
            "bin_right": float(first["bin_right"]),
            "present_in_all_selected_dt": bool(present_all),
        }

        observed_subject_counts = []
        for dt in selected_dts:
            dt_int = int(dt)
            gdt = group[group["dt"].eq(dt_int)]
            required = int(req_map[dt_int])
            available = int(avail_map[dt_int])
            if gdt.empty:
                n_subjects = 0
                ok = False
            else:
                n_subjects = int(pd.to_numeric(gdt["n_subjects"], errors="coerce").max())
                ok = n_subjects >= required

            observed_subject_counts.append(n_subjects)
            dt_ok[dt_int] = ok
            availability_fraction = n_subjects / float(available) if available > 0 else np.nan
            availability_fractions.append(availability_fraction)

            row["n_subjects_dt{}".format(dt_int)] = n_subjects
            row["available_subjects_dt{}".format(dt_int)] = available
            row["required_subjects_dt{}".format(dt_int)] = required
            row["coverage_fraction_of_available_dt{}".format(dt_int)] = availability_fraction
            row["meets_requirement_dt{}".format(dt_int)] = bool(ok)

        qualifies = present_all and all(dt_ok.values())
        row["min_subjects_across_dt"] = int(min(observed_subject_counts)) if observed_subject_counts else 0
        row["min_subject_coverage_fraction_across_dt"] = float(np.nanmin(availability_fractions))
        row["meets_subject_requirement"] = bool(all(dt_ok.values()))
        row["qualifies"] = bool(qualifies)
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values("bin_id").reset_index(drop=True)
    qualifying = summary.loc[summary["qualifies"], "bin_id"].astype(int).tolist()

    if explicit_bin_ids:
        missing = sorted(set(explicit_bin_ids) - set(qualifying))
        if missing:
            raise ValueError(
                "Explicit bins do not satisfy complete-dt/per-dt subject requirements: {}".format(
                    missing
                )
            )
        selected = sorted(set(int(value) for value in explicit_bin_ids))
    elif scope == "core":
        selected = largest_contiguous_run(qualifying)
    else:
        selected = sorted(qualifying)

    summary["selected_for_analysis"] = summary["bin_id"].isin(selected)
    return selected, summary


# =============================================================================
# Subject arrays and temporal metrics
# =============================================================================


def prepare_subject_arrays(
    subject_stats: pd.DataFrame,
    curves: pd.DataFrame,
    all_subjects: Sequence[str],
    selected_dts: Sequence[int],
    bin_ids: Sequence[int],
) -> List[Dict[str, object]]:
    subject_to_index = {subject: index for index, subject in enumerate(all_subjects)}
    dt_to_index = {int(dt): index for index, dt in enumerate(selected_dts)}
    shape = (len(selected_dts), len(all_subjects))

    bin_reference = (
        curves[curves["bin_id"].isin(bin_ids)][
            ["bin_id", "bin_left", "bin_mid", "bin_right"]
        ]
        .drop_duplicates("bin_id")
        .set_index("bin_id")
    )

    arrays: Dict[int, Dict[str, object]] = {}
    for bin_id in bin_ids:
        arrays[int(bin_id)] = {
            "n": np.zeros(shape, dtype=float),
            "sum_dx": np.zeros(shape, dtype=float),
            "sum_dx2": np.zeros(shape, dtype=float),
            "subject_mean_dx": np.full(shape, np.nan, dtype=float),
            "subject_msd": np.full(shape, np.nan, dtype=float),
        }

    selected = subject_stats[
        subject_stats["dt"].isin(selected_dts) & subject_stats["bin_id"].isin(bin_ids)
    ]
    for row in selected.itertuples(index=False):
        dt = int(row.dt)
        bin_id = int(row.bin_id)
        subject = str(row.subject)
        if dt not in dt_to_index or subject not in subject_to_index:
            continue
        i = dt_to_index[dt]
        j = subject_to_index[subject]
        target = arrays[bin_id]
        target["n"][i, j] = float(row.n_subject_bin)  # type: ignore[index]
        target["sum_dx"][i, j] = float(row.sum_dx)  # type: ignore[index]
        target["sum_dx2"][i, j] = float(row.sum_dx2)  # type: ignore[index]
        target["subject_mean_dx"][i, j] = float(row.subject_mean_dx)  # type: ignore[index]
        target["subject_msd"][i, j] = float(row.subject_msd)  # type: ignore[index]

    output: List[Dict[str, object]] = []
    for bin_id in bin_ids:
        info = bin_reference.loc[int(bin_id)]
        item = arrays[int(bin_id)]
        item.update(
            {
                "bin_id": int(bin_id),
                "bin_left": float(info["bin_left"]),
                "bin_mid": float(info["bin_mid"]),
                "bin_right": float(info["bin_right"]),
            }
        )
        output.append(item)
    return output


def temporal_metric_series(
    item: Dict[str, object],
    subject_multiplicity: np.ndarray,
    estimator: str,
    metric: str,
) -> np.ndarray:
    if estimator == "transition_weighted":
        n = np.asarray(item["n"], dtype=float)
        sum_dx = np.asarray(item["sum_dx"], dtype=float)
        sum_dx2 = np.asarray(item["sum_dx2"], dtype=float)
        total_n = n @ subject_multiplicity
        total_sum = sum_dx @ subject_multiplicity
        total_sum2 = sum_dx2 @ subject_multiplicity
        mean = np.divide(
            total_sum,
            total_n,
            out=np.full(total_n.shape, np.nan, dtype=float),
            where=total_n > 0,
        )
        msd = np.divide(
            total_sum2,
            total_n,
            out=np.full(total_n.shape, np.nan, dtype=float),
            where=total_n > 0,
        )
        variance = np.divide(
            total_sum2 - total_n * mean ** 2,
            total_n - 1.0,
            out=np.full(total_n.shape, np.nan, dtype=float),
            where=total_n > 1.0,
        )
    elif estimator == "equal_subject_weighted":
        subject_mean = np.asarray(item["subject_mean_dx"], dtype=float)
        subject_msd = np.asarray(item["subject_msd"], dtype=float)
        mean_valid = np.isfinite(subject_mean)
        msd_valid = np.isfinite(subject_msd)
        mean_denominator = mean_valid.astype(float) @ subject_multiplicity
        msd_denominator = msd_valid.astype(float) @ subject_multiplicity
        mean = np.divide(
            np.nan_to_num(subject_mean, nan=0.0) @ subject_multiplicity,
            mean_denominator,
            out=np.full(mean_denominator.shape, np.nan, dtype=float),
            where=mean_denominator > 0,
        )
        msd = np.divide(
            np.nan_to_num(subject_msd, nan=0.0) @ subject_multiplicity,
            msd_denominator,
            out=np.full(msd_denominator.shape, np.nan, dtype=float),
            where=msd_denominator > 0,
        )
        variance = msd - mean ** 2
    else:
        raise ValueError("Unsupported estimator: {}".format(estimator))

    # Numerical round-off can produce tiny negative equal-subject variances.
    variance = np.where((variance < 0) & (variance > -1e-12), 0.0, variance)
    if metric == "msd":
        return msd
    if metric == "var_dx":
        return variance
    raise ValueError("Unsupported metric: {}".format(metric))


# =============================================================================
# Model fitting
# =============================================================================


def fit_statistics(
    y: np.ndarray,
    predicted: np.ndarray,
    n_structural_parameters: int,
) -> Dict[str, float]:
    """
    Compute Gaussian-residual AIC/AICc.

    AIC parameter count includes the residual-variance parameter in addition to
    the structural model parameters. This is consequential with only five lags.
    """
    valid = np.isfinite(y) & np.isfinite(predicted)
    y_valid = y[valid]
    predicted_valid = predicted[valid]
    n = int(y_valid.size)
    if n == 0:
        return {
            "rss": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "aic": np.nan,
            "aicc": np.nan,
            "n_parameters_structural": int(n_structural_parameters),
            "n_parameters_aic": int(n_structural_parameters + 1),
        }

    residual = y_valid - predicted_valid
    rss = float(np.sum(residual ** 2))
    rmse = float(np.sqrt(rss / n))
    tss = float(np.sum((y_valid - np.mean(y_valid)) ** 2))
    r2 = float(1.0 - rss / tss) if tss > 0 else np.nan
    safe_rss = max(rss, np.finfo(float).tiny)

    k = int(n_structural_parameters) + 1  # + residual variance
    aic = float(n * np.log(safe_rss / n) + 2.0 * k)
    if n > k + 1:
        aicc = float(
            aic + (2.0 * k * (k + 1.0)) / (n - k - 1.0)
        )
    else:
        aicc = np.nan

    return {
        "rss": rss,
        "rmse": rmse,
        "r2": r2,
        "aic": aic,
        "aicc": aicc,
        "n_parameters_structural": int(n_structural_parameters),
        "n_parameters_aic": k,
    }


def fit_constant(dt: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    """Fit an empirical constant level over the observed dt range only."""
    valid = np.isfinite(dt) & np.isfinite(y) & (dt > 0) & (y >= 0)
    x = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x.size < 2:
        return {"fit_success": False, "fit_error": "fewer than 2 valid dt values"}
    level = max(0.0, float(np.mean(values)))
    predicted = np.full(values.shape, level, dtype=float)
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "level": level,
        "predicted": predicted,
    }
    output.update(fit_statistics(values, predicted, 1))
    return output


def fit_linear_zero(dt: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    """Fit the theoretical zero-origin Brownian comparator M(dt)=D*dt."""
    valid = np.isfinite(dt) & np.isfinite(y) & (dt > 0) & (y >= 0)
    x = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x.size < 2:
        return {"fit_success": False, "fit_error": "fewer than 2 valid dt values"}
    denominator = float(np.sum(x ** 2))
    if denominator <= 0:
        return {"fit_success": False, "fit_error": "non-positive denominator"}
    slope = max(0.0, float(np.sum(x * values) / denominator))
    predicted = slope * x
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "slope": slope,
        "predicted": predicted,
    }
    output.update(fit_statistics(values, predicted, 1))
    return output


def fit_incremental_linear(dt: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    """
    Fit M(dt)=K + D*(dt-1), constrained to K>=0 and D>=0.

    This is the direct nested alternative to the constant model for detecting
    additional accumulation beyond the minimum observed one-week lag.
    """
    valid = np.isfinite(dt) & np.isfinite(y) & (dt >= 1) & (y >= 0)
    x_dt = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x_dt.size < 3:
        return {"fit_success": False, "fit_error": "fewer than 3 valid dt values"}

    x = x_dt - 1.0
    X = np.column_stack([np.ones_like(x), x])

    candidates = []
    try:
        beta, _, _, _ = np.linalg.lstsq(X, values, rcond=None)
        if beta[0] >= 0 and beta[1] >= 0:
            pred = X @ beta
            candidates.append((float(np.sum((values - pred) ** 2)), float(beta[0]), float(beta[1]), pred))
    except np.linalg.LinAlgError:
        pass

    # Boundary D=0 -> constant.
    k_only = max(0.0, float(np.mean(values)))
    pred_k = np.full(values.shape, k_only, dtype=float)
    candidates.append((float(np.sum((values - pred_k) ** 2)), k_only, 0.0, pred_k))

    # Boundary K=0 -> slope on (dt-1).
    denom = float(np.sum(x ** 2))
    d_only = max(0.0, float(np.sum(x * values) / denom)) if denom > 0 else 0.0
    pred_d = d_only * x
    candidates.append((float(np.sum((values - pred_d) ** 2)), 0.0, d_only, pred_d))

    pred_zero = np.zeros(values.shape, dtype=float)
    candidates.append((float(np.sum(values ** 2)), 0.0, 0.0, pred_zero))

    rss, intercept, slope, predicted = min(candidates, key=lambda item: item[0])
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "intercept_at_dt1": intercept,
        "slope": slope,
        "predicted": predicted,
    }
    output.update(fit_statistics(values, predicted, 2))
    output["rss"] = rss
    return output


def fit_signed_linear(dt: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    """
    Fit M(dt)=K + D*(dt-1) with signed D and non-negative fitted values.

    D is allowed to be positive or negative. Because MSD and variance are
    non-negative quantities, the fitted line is constrained to remain >= 0
    over the observed lag range. For a linear function this is equivalent to
    requiring non-negative fitted values at the two endpoints.
    """
    valid = np.isfinite(dt) & np.isfinite(y) & (dt >= 1) & (y >= 0)
    x_dt = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x_dt.size < 3:
        return {"fit_success": False, "fit_error": "fewer than 3 valid dt values"}

    x = x_dt - 1.0
    X = np.column_stack([np.ones_like(x), x])
    x_max = float(np.max(x))
    candidates = []

    # Interior ordinary least-squares solution when physically admissible.
    try:
        beta, _, _, _ = np.linalg.lstsq(X, values, rcond=None)
        pred = X @ beta
        if beta[0] >= -1e-12 and np.all(pred >= -1e-12):
            intercept = max(0.0, float(beta[0]))
            slope = float(beta[1])
            pred = intercept + slope * x
            candidates.append((float(np.sum((values - pred) ** 2)), intercept, slope, pred, "interior"))
    except np.linalg.LinAlgError:
        pass

    # Boundary K=0. Non-negativity then requires D>=0 for x>=0.
    denom = float(np.sum(x ** 2))
    slope_k0 = max(0.0, float(np.sum(x * values) / denom)) if denom > 0 else 0.0
    pred_k0 = slope_k0 * x
    candidates.append((float(np.sum((values - pred_k0) ** 2)), 0.0, slope_k0, pred_k0, "K_at_zero"))

    # Boundary at the longest observed lag: M(dt_max)=0. This permits D<0.
    if x_max > 0:
        basis = 1.0 - x / x_max
        denom_end = float(np.sum(basis ** 2))
        intercept_end = (
            max(0.0, float(np.sum(basis * values) / denom_end))
            if denom_end > 0 else 0.0
        )
        slope_end = -intercept_end / x_max
        pred_end = intercept_end + slope_end * x
        candidates.append((float(np.sum((values - pred_end) ** 2)), intercept_end, slope_end, pred_end, "dtmax_at_zero"))

    # Degenerate zero solution.
    pred_zero = np.zeros(values.shape, dtype=float)
    candidates.append((float(np.sum(values ** 2)), 0.0, 0.0, pred_zero, "all_zero"))

    rss, intercept, slope, predicted, boundary = min(candidates, key=lambda item: item[0])
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "intercept_at_dt1": float(intercept),
        "slope": float(slope),
        "predicted": predicted,
        "predicted_min": float(np.min(predicted)),
        "constraint_status": boundary,
    }
    output.update(fit_statistics(values, predicted, 2))
    output["rss"] = rss
    return output


def power_profile(dt: np.ndarray, y: np.ndarray, alpha: float) -> Tuple[float, float, np.ndarray]:
    basis = np.power(dt, alpha)
    denominator = float(np.sum(basis ** 2))
    coefficient = (
        max(0.0, float(np.sum(basis * y) / denominator)) if denominator > 0 else 0.0
    )
    predicted = coefficient * basis
    rss = float(np.sum((y - predicted) ** 2))
    return rss, coefficient, predicted


def golden_section(
    objective,
    low: float,
    high: float,
    max_iter: int = 100,
    tolerance: float = 1e-10,
) -> float:
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    a = float(low)
    b = float(high)
    c = b - golden * (b - a)
    d = a + golden * (b - a)
    fc = float(objective(c))
    fd = float(objective(d))
    for _ in range(max_iter):
        if abs(b - a) <= tolerance:
            break
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - golden * (b - a)
            fc = float(objective(c))
        else:
            a, c, fc = c, d, fd
            d = a + golden * (b - a)
            fd = float(objective(d))
    return float(0.5 * (a + b))


def fit_power_law(
    dt: np.ndarray,
    y: np.ndarray,
    alpha_min: float,
    alpha_max: float,
) -> Dict[str, object]:
    valid = np.isfinite(dt) & np.isfinite(y) & (dt > 0) & (y >= 0)
    x = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x.size < 3:
        return {"fit_success": False, "fit_error": "fewer than 3 valid dt values"}
    if not (math.isfinite(alpha_min) and math.isfinite(alpha_max)):
        return {"fit_success": False, "fit_error": "non-finite alpha bounds"}
    if alpha_min <= 0 or alpha_max <= alpha_min:
        return {"fit_success": False, "fit_error": "invalid alpha bounds"}

    grid = np.linspace(alpha_min, alpha_max, 257)
    rss_grid = np.array([power_profile(x, values, float(a))[0] for a in grid])
    best_index = int(np.nanargmin(rss_grid))
    left_index = max(0, best_index - 1)
    right_index = min(len(grid) - 1, best_index + 1)
    if left_index == right_index:
        alpha = float(grid[best_index])
    else:
        alpha = golden_section(
            lambda a: power_profile(x, values, float(a))[0],
            float(grid[left_index]),
            float(grid[right_index]),
        )

    rss, coefficient, predicted = power_profile(x, values, alpha)
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "coefficient": float(coefficient),
        "alpha": float(alpha),
        "alpha_at_lower_bound": bool(alpha <= alpha_min + 0.01 * (alpha_max - alpha_min)),
        "alpha_at_upper_bound": bool(alpha >= alpha_max - 0.01 * (alpha_max - alpha_min)),
        "predicted": predicted,
    }
    output.update(fit_statistics(values, predicted, 2))
    output["rss"] = rss
    return output


def saturating_profile(dt: np.ndarray, y: np.ndarray, tau: float) -> Tuple[float, float, np.ndarray]:
    basis = 1.0 - np.exp(-dt / tau)
    denominator = float(np.sum(basis ** 2))
    amplitude = (
        max(0.0, float(np.sum(basis * y) / denominator)) if denominator > 0 else 0.0
    )
    predicted = amplitude * basis
    rss = float(np.sum((y - predicted) ** 2))
    return rss, amplitude, predicted


def fit_saturating_zero(
    dt: np.ndarray,
    y: np.ndarray,
    tau_min: float,
    tau_max: float,
) -> Dict[str, object]:
    valid = np.isfinite(dt) & np.isfinite(y) & (dt > 0) & (y >= 0)
    x = dt[valid].astype(float)
    values = y[valid].astype(float)
    if x.size < 3:
        return {"fit_success": False, "fit_error": "fewer than 3 valid dt values"}
    if not (math.isfinite(tau_min) and math.isfinite(tau_max)):
        return {"fit_success": False, "fit_error": "non-finite tau bounds"}
    if tau_min <= 0 or tau_max <= tau_min:
        return {"fit_success": False, "fit_error": "invalid tau bounds"}

    log_grid = np.linspace(math.log(tau_min), math.log(tau_max), 257)
    rss_grid = np.array(
        [saturating_profile(x, values, float(math.exp(v)))[0] for v in log_grid]
    )
    best_index = int(np.nanargmin(rss_grid))
    left_index = max(0, best_index - 1)
    right_index = min(len(log_grid) - 1, best_index + 1)
    if left_index == right_index:
        log_tau = float(log_grid[best_index])
    else:
        log_tau = golden_section(
            lambda value: saturating_profile(x, values, float(math.exp(value)))[0],
            float(log_grid[left_index]),
            float(log_grid[right_index]),
        )
    tau = float(math.exp(log_tau))
    rss, amplitude, predicted = saturating_profile(x, values, tau)
    output: Dict[str, object] = {
        "fit_success": True,
        "fit_error": "",
        "amplitude": float(amplitude),
        "tau": tau,
        "tau_at_lower_bound": bool(tau <= tau_min * 1.01),
        "tau_at_upper_bound": bool(tau >= tau_max * 0.99),
        "predicted": predicted,
    }
    output.update(fit_statistics(values, predicted, 2))
    output["rss"] = rss
    return output


def fit_six_models(
    dt: np.ndarray,
    y: np.ndarray,
    min_dt_points: int,
    alpha_min: float,
    alpha_max: float,
    tau_min: float,
    tau_max: float,
) -> Dict[str, object]:
    n_valid = int(np.sum(np.isfinite(dt) & np.isfinite(y) & (dt > 0) & (y >= 0)))
    if n_valid < min_dt_points:
        return {
            "n_dt": n_valid,
            "all_models_success": False,
            "preferred_model": "unresolved",
            "fit_error": "fewer than {} valid dt values".format(min_dt_points),
        }

    fits = {
        "constant": fit_constant(dt, y),
        "signed_linear": fit_signed_linear(dt, y),
        "incremental_linear": fit_incremental_linear(dt, y),
        "linear": fit_linear_zero(dt, y),
        "power_law": fit_power_law(dt, y, alpha_min, alpha_max),
        "saturating": fit_saturating_zero(dt, y, tau_min, tau_max),
    }
    output: Dict[str, object] = {"n_dt": n_valid, "fit_error": ""}
    for model, fit in fits.items():
        output["{}_fit_success".format(model)] = bool(fit.get("fit_success", False))
        output["{}_fit_error".format(model)] = str(fit.get("fit_error", ""))
        for statistic in (
            "rss", "rmse", "r2", "aic", "aicc",
            "n_parameters_structural", "n_parameters_aic",
        ):
            value = fit.get(statistic, np.nan)
            output["{}_{}".format(model, statistic)] = (
                float(value) if np.isscalar(value) and value is not None else np.nan
            )

    output["constant_level"] = float(fits["constant"].get("level", np.nan))
    output["linear_slope"] = float(fits["linear"].get("slope", np.nan))
    output["signed_linear_intercept_at_dt1"] = float(
        fits["signed_linear"].get("intercept_at_dt1", np.nan)
    )
    output["signed_linear_slope"] = float(
        fits["signed_linear"].get("slope", np.nan)
    )
    output["signed_linear_constraint_status"] = str(
        fits["signed_linear"].get("constraint_status", "")
    )
    output["incremental_linear_intercept_at_dt1"] = float(
        fits["incremental_linear"].get("intercept_at_dt1", np.nan)
    )
    output["incremental_linear_slope"] = float(
        fits["incremental_linear"].get("slope", np.nan)
    )
    output["power_coefficient"] = float(fits["power_law"].get("coefficient", np.nan))
    output["power_alpha"] = float(fits["power_law"].get("alpha", np.nan))
    output["power_alpha_at_lower_bound"] = bool(
        fits["power_law"].get("alpha_at_lower_bound", False)
    )
    output["power_alpha_at_upper_bound"] = bool(
        fits["power_law"].get("alpha_at_upper_bound", False)
    )
    output["saturating_amplitude"] = float(fits["saturating"].get("amplitude", np.nan))
    output["saturating_tau"] = float(fits["saturating"].get("tau", np.nan))
    output["saturating_tau_at_lower_bound"] = bool(
        fits["saturating"].get("tau_at_lower_bound", False)
    )
    output["saturating_tau_at_upper_bound"] = bool(
        fits["saturating"].get("tau_at_upper_bound", False)
    )

    all_success = all(bool(fits[m].get("fit_success", False)) for m in MODELS)
    output["all_models_success"] = all_success
    if not all_success:
        errors = [
            "{}: {}".format(model, fits[model].get("fit_error", "failed"))
            for model in MODELS
            if not bool(fits[model].get("fit_success", False))
        ]
        output["fit_error"] = "; ".join(errors)
        output["preferred_model"] = "unresolved"
    else:
        scores = {model: float(output["{}_aicc".format(model)]) for model in MODELS}
        finite_scores = {
            model: score for model, score in scores.items() if np.isfinite(score)
        }
        if len(finite_scores) != len(MODELS):
            output["preferred_model"] = "unresolved"
        else:
            minimum = min(finite_scores.values())
            winners = [
                model for model, score in finite_scores.items()
                if abs(score - minimum) <= 1e-10
            ]
            output["preferred_model"] = (
                winners[0] if len(winners) == 1 else "unresolved_tie"
            )

    # Positive delta means the model named first has larger (worse) AICc.
    comparison_pairs = [
        ("signed_linear", "constant"),
        ("incremental_linear", "constant"),
        ("linear", "constant"),
        ("power_law", "constant"),
        ("saturating", "constant"),
        ("incremental_linear", "signed_linear"),
        ("linear", "signed_linear"),
        ("power_law", "signed_linear"),
        ("saturating", "signed_linear"),
        ("linear", "incremental_linear"),
        ("power_law", "incremental_linear"),
        ("saturating", "incremental_linear"),
        ("linear", "power_law"),
        ("linear", "saturating"),
        ("power_law", "saturating"),
    ]
    for first, second in comparison_pairs:
        output["delta_aicc_{}_minus_{}".format(first, second)] = float(
            output.get("{}_aicc".format(first), np.nan)
            - output.get("{}_aicc".format(second), np.nan)
        )

    for model in MODELS:
        output["_predicted_{}".format(model)] = fits[model].get("predicted")
    return output


# =============================================================================
# Observed fits and joint bootstrap
# =============================================================================


def common_base(
    dataset_label: str,
    estimator: str,
    metric: str,
    item: Dict[str, object],
    bin_scope: str,
    core_min_subjects: int,
) -> Dict[str, object]:
    return {
        "dataset_label": dataset_label,
        "estimator": estimator,
        "metric": metric,
        "bin_scope": bin_scope,
        "core_min_subjects": int(core_min_subjects),
        "bin_id": int(item["bin_id"]),
        "bin_left": float(item["bin_left"]),
        "bin_mid": float(item["bin_mid"]),
        "bin_right": float(item["bin_right"]),
    }


def observed_fits(
    arrays_by_bin: Sequence[Dict[str, object]],
    all_subjects: Sequence[str],
    selected_dts: Sequence[int],
    estimators: Sequence[str],
    metrics: Sequence[str],
    dataset_label: str,
    args: argparse.Namespace,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    dt = np.asarray(selected_dts, dtype=float)
    multiplicity = np.ones(len(all_subjects), dtype=float)
    pair_rows: List[Dict[str, object]] = []
    long_rows: List[Dict[str, object]] = []

    for item in arrays_by_bin:
        for estimator in estimators:
            for metric in metrics:
                values = temporal_metric_series(item, multiplicity, estimator, metric)
                fit = fit_six_models(
                    dt,
                    values,
                    args.min_dt_points,
                    args.alpha_min,
                    args.alpha_max,
                    args.tau_min,
                    args.tau_max,
                )
                base = common_base(
                    dataset_label,
                    estimator,
                    metric,
                    item,
                    args.bin_scope,
                    args.core_min_subjects_resolved,
                )
                pair = dict(base)
                pair.update({key: value for key, value in fit.items() if not key.startswith("_predicted")})
                pair["dt_values"] = ",".join(str(value) for value in selected_dts)
                pair["metric_values"] = ",".join(
                    "{:.12g}".format(value) if np.isfinite(value) else "nan"
                    for value in values
                )
                pair_rows.append(pair)

                for model in MODELS:
                    row = dict(base)
                    row.update(
                        {
                            "model": model,
                            "n_dt": int(fit.get("n_dt", 0)),
                            "fit_success": bool(fit.get("{}_fit_success".format(model), False)),
                            "fit_error": str(fit.get("{}_fit_error".format(model), "")),
                            "rss": fit.get("{}_rss".format(model), np.nan),
                            "rmse": fit.get("{}_rmse".format(model), np.nan),
                            "r2": fit.get("{}_r2".format(model), np.nan),
                            "aic": fit.get("{}_aic".format(model), np.nan),
                            "aicc": fit.get("{}_aicc".format(model), np.nan),
                            "best_model": bool(fit.get("preferred_model") == model),
                            "preferred_model": fit.get("preferred_model", "unresolved"),
                            "level": fit.get("constant_level", np.nan) if model == "constant" else np.nan,
                            "slope": (
                                fit.get("linear_slope", np.nan)
                                if model == "linear"
                                else fit.get("signed_linear_slope", np.nan)
                                if model == "signed_linear"
                                else fit.get("incremental_linear_slope", np.nan)
                                if model == "incremental_linear"
                                else np.nan
                            ),
                            "intercept_at_dt1": (
                                fit.get("signed_linear_intercept_at_dt1", np.nan)
                                if model == "signed_linear"
                                else fit.get("incremental_linear_intercept_at_dt1", np.nan)
                                if model == "incremental_linear"
                                else np.nan
                            ),
                            "coefficient": fit.get("power_coefficient", np.nan) if model == "power_law" else np.nan,
                            "alpha": fit.get("power_alpha", np.nan) if model == "power_law" else np.nan,
                            "alpha_at_lower_bound": fit.get("power_alpha_at_lower_bound", False) if model == "power_law" else False,
                            "alpha_at_upper_bound": fit.get("power_alpha_at_upper_bound", False) if model == "power_law" else False,
                            "amplitude": fit.get("saturating_amplitude", np.nan) if model == "saturating" else np.nan,
                            "tau": fit.get("saturating_tau", np.nan) if model == "saturating" else np.nan,
                            "tau_at_lower_bound": fit.get("saturating_tau_at_lower_bound", False) if model == "saturating" else False,
                            "tau_at_upper_bound": fit.get("saturating_tau_at_upper_bound", False) if model == "saturating" else False,
                            **{
                                key: fit.get(key, np.nan)
                                for key in fit
                                if str(key).startswith("delta_aicc_")
                            },
                            "dt_values": ",".join(str(value) for value in selected_dts),
                            "metric_values": pair["metric_values"],
                            "predicted_values": ",".join(
                                "{:.12g}".format(value) if np.isfinite(value) else "nan"
                                for value in np.asarray(
                                    fit.get("_predicted_{}".format(model), np.full(len(selected_dts), np.nan)),
                                    dtype=float,
                                )
                            ),
                        }
                    )
                    long_rows.append(row)
    return pd.DataFrame(long_rows), pd.DataFrame(pair_rows)


def joint_subject_bootstrap(
    arrays_by_bin: Sequence[Dict[str, object]],
    all_subjects: Sequence[str],
    selected_dts: Sequence[int],
    estimators: Sequence[str],
    metrics: Sequence[str],
    dataset_label: str,
    args: argparse.Namespace,
) -> pd.DataFrame:
    if args.n_boot <= 0:
        return pd.DataFrame()
    dt = np.asarray(selected_dts, dtype=float)
    n_subjects = len(all_subjects)
    rng = np.random.default_rng(args.seed)
    rows: List[Dict[str, object]] = []

    for bootstrap_id in range(args.n_boot):
        sampled_indices = rng.integers(0, n_subjects, size=n_subjects)
        multiplicity = np.bincount(sampled_indices, minlength=n_subjects).astype(float)
        n_unique = int(np.count_nonzero(multiplicity))

        for item in arrays_by_bin:
            for estimator in estimators:
                for metric in metrics:
                    values = temporal_metric_series(item, multiplicity, estimator, metric)
                    fit = fit_six_models(
                        dt,
                        values,
                        args.min_dt_points,
                        args.alpha_min,
                        args.alpha_max,
                        args.tau_min,
                        args.tau_max,
                    )
                    row = common_base(
                        dataset_label,
                        estimator,
                        metric,
                        item,
                        args.bin_scope,
                        args.core_min_subjects_resolved,
                    )
                    row.update(
                        {
                            "bootstrap_id": int(bootstrap_id),
                            "n_subject_draws": int(n_subjects),
                            "n_unique_subjects_sampled": n_unique,
                        }
                    )
                    row.update(
                        {key: value for key, value in fit.items() if not key.startswith("_predicted")}
                    )
                    rows.append(row)

        if (bootstrap_id + 1) % max(1, min(100, args.n_boot // 10)) == 0:
            print("[BOOTSTRAP] {}/{}".format(bootstrap_id + 1, args.n_boot))

    return pd.DataFrame(rows).sort_values(
        ["bootstrap_id", "estimator", "metric", "bin_id"]
    ).reset_index(drop=True)


# =============================================================================
# Summaries and decision aid
# =============================================================================


def quantile(values: pd.Series, probability: float) -> float:
    x = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.quantile(probability)) if len(x) else np.nan


def summarize_by_bin(
    observed_pairs: pd.DataFrame,
    bootstrap: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    keys = [
        "dataset_label",
        "estimator",
        "metric",
        "bin_scope",
        "core_min_subjects",
        "bin_id",
        "bin_left",
        "bin_mid",
        "bin_right",
    ]
    alpha_low = (1.0 - args.ci / 100.0) / 2.0
    alpha_high = 1.0 - alpha_low
    rows: List[Dict[str, object]] = []

    for group_key, group in bootstrap.groupby(keys, observed=True, sort=True):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        base = dict(zip(keys, group_key))
        observed = observed_pairs
        for key, value in base.items():
            observed = observed[observed[key] == value]
        obs = observed.iloc[0] if len(observed) else pd.Series(dtype=object)

        valid_mask = group["all_models_success"].astype(bool)
        for model in MODELS:
            valid_mask &= pd.to_numeric(
                group["{}_aicc".format(model)], errors="coerce"
            ).notna()
        valid = group[valid_mask].copy()
        n_valid = int(len(valid))

        fractions = {
            model: (
                float((valid["preferred_model"].astype(str) == model).mean())
                if n_valid else np.nan
            )
            for model in MODELS
        }

        row: Dict[str, object] = dict(base)
        row.update(
            {
                "n_bootstrap_requested": int(args.n_boot),
                "n_valid_bootstrap": n_valid,
                "valid_bootstrap_fraction": (
                    n_valid / float(args.n_boot) if args.n_boot > 0 else np.nan
                ),
                "observed_preferred_model": obs.get("preferred_model", "unresolved"),
            }
        )

        for model in MODELS:
            row["observed_aicc_{}".format(model)] = obs.get(
                "{}_aicc".format(model), np.nan
            )
            row["{}_preference_fraction".format(model)] = fractions[model]

        # Constant-referenced effect sizes are the primary model-comparison summaries.
        for model in ("signed_linear", "incremental_linear", "linear", "power_law", "saturating"):
            column = "delta_aicc_{}_minus_constant".format(model)
            values = pd.to_numeric(valid.get(column), errors="coerce")
            row["observed_{}".format(column)] = obs.get(column, np.nan)
            row["{}_median".format(column)] = quantile(values, 0.5)
            row["{}_q_low".format(column)] = quantile(values, alpha_low)
            row["{}_q_high".format(column)] = quantile(values, alpha_high)
            row["constant_beats_{}_fraction".format(model)] = (
                float((values > 0).mean()) if len(values) else np.nan
            )
            row["constant_beats_{}_by_2_fraction".format(model)] = (
                float((values > args.delta_aicc_threshold).mean())
                if len(values) else np.nan
            )

        # Signed-linear trend diagnostic: D may be positive or negative.
        signed_delta = pd.to_numeric(
            valid.get("delta_aicc_signed_linear_minus_constant"),
            errors="coerce",
        )
        row["signed_linear_beats_constant_fraction"] = (
            float((signed_delta < 0).mean()) if len(signed_delta) else np.nan
        )
        row["signed_linear_beats_constant_by_2_fraction"] = (
            float((signed_delta < -args.delta_aicc_threshold).mean())
            if len(signed_delta) else np.nan
        )
        signed_slopes = pd.to_numeric(
            valid.get("signed_linear_slope"), errors="coerce"
        )
        row["observed_signed_linear_slope"] = obs.get(
            "signed_linear_slope", np.nan
        )
        row["signed_linear_slope_median"] = quantile(signed_slopes, 0.5)
        row["signed_linear_slope_q_low"] = quantile(signed_slopes, alpha_low)
        row["signed_linear_slope_q_high"] = quantile(signed_slopes, alpha_high)
        row["signed_linear_positive_slope_fraction"] = (
            float((signed_slopes > 0).mean()) if signed_slopes.notna().any() else np.nan
        )
        row["signed_linear_negative_slope_fraction"] = (
            float((signed_slopes < 0).mean()) if signed_slopes.notna().any() else np.nan
        )
        row["signed_linear_slope_ci_below_zero"] = bool(
            np.isfinite(row["signed_linear_slope_q_high"])
            and row["signed_linear_slope_q_high"] < 0
        )
        row["signed_linear_slope_ci_above_zero"] = bool(
            np.isfinite(row["signed_linear_slope_q_low"])
            and row["signed_linear_slope_q_low"] > 0
        )

        # Direct nested comparison: constant vs positive incremental accumulation.
        inc_delta = pd.to_numeric(
            valid.get("delta_aicc_incremental_linear_minus_constant"),
            errors="coerce",
        )
        row["incremental_linear_beats_constant_fraction"] = (
            float((inc_delta < 0).mean()) if len(inc_delta) else np.nan
        )
        row["incremental_linear_beats_constant_by_2_fraction"] = (
            float((inc_delta < -args.delta_aicc_threshold).mean())
            if len(inc_delta) else np.nan
        )

        inc_slopes = pd.to_numeric(
            valid.get("incremental_linear_slope"), errors="coerce"
        )
        row["observed_incremental_linear_slope"] = obs.get(
            "incremental_linear_slope", np.nan
        )
        row["incremental_linear_slope_median"] = quantile(inc_slopes, 0.5)
        row["incremental_linear_slope_q_low"] = quantile(inc_slopes, alpha_low)
        row["incremental_linear_slope_q_high"] = quantile(inc_slopes, alpha_high)
        row["incremental_linear_positive_slope_fraction"] = (
            float((inc_slopes > 0).mean()) if inc_slopes.notna().any() else np.nan
        )

        alpha_values = pd.to_numeric(valid.get("power_alpha"), errors="coerce")
        tau_values = pd.to_numeric(valid.get("saturating_tau"), errors="coerce")
        power_preferred = valid["preferred_model"].astype(str).eq("power_law")
        alpha_when_power = alpha_values[power_preferred]

        row.update(
            {
                "observed_power_alpha": obs.get("power_alpha", np.nan),
                "power_alpha_median": quantile(alpha_values, 0.5),
                "power_alpha_q_low": quantile(alpha_values, alpha_low),
                "power_alpha_q_high": quantile(alpha_values, alpha_high),
                "power_alpha_lt_1_fraction": (
                    float((alpha_values < 1.0).mean())
                    if alpha_values.notna().any() else np.nan
                ),
                "power_alpha_lt_1_when_power_preferred_fraction": (
                    float((alpha_when_power < 1.0).mean())
                    if alpha_when_power.notna().any() else np.nan
                ),
                "power_alpha_ci_below_1": (
                    bool(quantile(alpha_values, alpha_high) < 1.0)
                    if alpha_values.notna().any() else False
                ),
                "power_alpha_at_lower_bound_fraction": (
                    float(valid["power_alpha_at_lower_bound"].astype(bool).mean())
                    if n_valid else np.nan
                ),
                "power_alpha_at_upper_bound_fraction": (
                    float(valid["power_alpha_at_upper_bound"].astype(bool).mean())
                    if n_valid else np.nan
                ),
                "observed_saturating_tau": obs.get("saturating_tau", np.nan),
                "saturating_tau_median": quantile(tau_values, 0.5),
                "saturating_tau_q_low": quantile(tau_values, alpha_low),
                "saturating_tau_q_high": quantile(tau_values, alpha_high),
                "saturating_tau_at_lower_bound_fraction": (
                    float(valid["saturating_tau_at_lower_bound"].astype(bool).mean())
                    if n_valid else np.nan
                ),
                "saturating_tau_at_upper_bound_fraction": (
                    float(valid["saturating_tau_at_upper_bound"].astype(bool).mean())
                    if n_valid else np.nan
                ),
                "ou_second_moment_compatibility_only": True,
                "ou_restoring_drift_tested": False,
            }
        )

        power_parameter_identifiable = bool(
            np.isfinite(row["power_alpha_at_lower_bound_fraction"])
            and np.isfinite(row["power_alpha_at_upper_bound_fraction"])
            and float(row["power_alpha_at_lower_bound_fraction"])
            < args.parameter_bound_warning_fraction
            and float(row["power_alpha_at_upper_bound_fraction"])
            < args.parameter_bound_warning_fraction
        )
        saturating_parameter_identifiable = bool(
            np.isfinite(row["saturating_tau_at_lower_bound_fraction"])
            and np.isfinite(row["saturating_tau_at_upper_bound_fraction"])
            and float(row["saturating_tau_at_lower_bound_fraction"])
            < args.parameter_bound_warning_fraction
            and float(row["saturating_tau_at_upper_bound_fraction"])
            < args.parameter_bound_warning_fraction
        )
        row["power_parameter_identifiable"] = power_parameter_identifiable
        row["saturating_parameter_identifiable"] = saturating_parameter_identifiable

        valid_fraction = row["valid_bootstrap_fraction"]
        signed_negative = row["signed_linear_negative_slope_fraction"]
        signed_positive = row["signed_linear_positive_slope_fraction"]
        if not np.isfinite(valid_fraction) or valid_fraction < 0.80:
            classification = "fit_unstable"
        elif (
            fractions["signed_linear"] >= args.model_preference_threshold
            and np.isfinite(signed_negative)
            and signed_negative >= args.signed_slope_probability_threshold
        ):
            classification = "signed_linear_negative_trend_supported"
        elif (
            fractions["signed_linear"] >= args.model_preference_threshold
            and np.isfinite(signed_positive)
            and signed_positive >= args.signed_slope_probability_threshold
        ):
            classification = "signed_linear_positive_trend_supported"
        elif fractions["constant"] >= args.model_preference_threshold:
            classification = "constant_level_preferred"
        elif fractions["incremental_linear"] >= args.model_preference_threshold:
            classification = "positive_incremental_accumulation_supported"
        elif fractions["linear"] >= args.model_preference_threshold:
            classification = "zero_origin_linear_scaling_supported"
        elif (
            fractions["power_law"] >= args.model_preference_threshold
            and row["power_alpha_lt_1_fraction"]
            >= args.subdiffusion_probability_threshold
            and power_parameter_identifiable
        ):
            classification = "subdiffusive_power_law_supported"
        elif (
            fractions["saturating"] >= args.model_preference_threshold
            and saturating_parameter_identifiable
        ):
            classification = "exponential_saturation_supported"
        else:
            classification = "inconclusive"

        row["operational_classification"] = classification

        # Slope-sign evidence is reported separately from AICc model preference.
        # With only five lags, AICc can prefer the constant model even when the
        # subject-bootstrap distribution of the signed slope is directionally
        # concentrated. Keeping these two inferential layers separate prevents
        # the small-sample AICc penalty from being mistaken for a zero slope.
        if (
            np.isfinite(signed_negative)
            and signed_negative >= args.signed_slope_probability_threshold
        ):
            signed_trend_classification = "negative_slope_supported"
        elif (
            np.isfinite(signed_positive)
            and signed_positive >= args.signed_slope_probability_threshold
        ):
            signed_trend_classification = "positive_slope_supported"
        else:
            signed_trend_classification = "signed_slope_inconclusive"
        row["signed_trend_classification"] = signed_trend_classification
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["estimator", "metric", "bin_id"]
    ).reset_index(drop=True)


def summarize_overall(by_bin: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for (dataset, estimator, metric), group in by_bin.groupby(
        ["dataset_label", "estimator", "metric"], observed=True, sort=True
    ):
        classifications = group["operational_classification"].astype(str)
        signed_trends = group["signed_trend_classification"].astype(str)
        n_bins = int(len(group))
        row: Dict[str, object] = {
            "dataset_label": dataset,
            "estimator": estimator,
            "metric": metric,
            "n_bins_tested": n_bins,
            "tested_x_min": float(group["bin_left"].min()),
            "tested_x_max": float(group["bin_right"].max()),
            "n_observed_constant_best": int(
                group["observed_preferred_model"].astype(str).eq("constant").sum()
            ),
            "n_observed_linear_best": int(
                group["observed_preferred_model"].astype(str).eq("linear").sum()
            ),
            "n_observed_signed_linear_best": int(
                group["observed_preferred_model"].astype(str).eq("signed_linear").sum()
            ),
            "n_observed_incremental_linear_best": int(
                group["observed_preferred_model"].astype(str)
                .eq("incremental_linear").sum()
            ),
            "n_observed_power_law_best": int(
                group["observed_preferred_model"].astype(str).eq("power_law").sum()
            ),
            "n_observed_saturating_best": int(
                group["observed_preferred_model"].astype(str).eq("saturating").sum()
            ),
            "median_constant_preference_fraction": float(
                group["constant_preference_fraction"].median()
            ),
            "median_linear_preference_fraction": float(
                group["linear_preference_fraction"].median()
            ),
            "median_signed_linear_preference_fraction": float(
                group["signed_linear_preference_fraction"].median()
            ),
            "median_incremental_linear_preference_fraction": float(
                group["incremental_linear_preference_fraction"].median()
            ),
            "median_power_law_preference_fraction": float(
                group["power_law_preference_fraction"].median()
            ),
            "median_saturating_preference_fraction": float(
                group["saturating_preference_fraction"].median()
            ),
            "median_signed_linear_slope": float(
                group["signed_linear_slope_median"].median()
            ),
            "median_signed_linear_negative_slope_fraction": float(
                group["signed_linear_negative_slope_fraction"].median()
            ),
            "median_signed_linear_positive_slope_fraction": float(
                group["signed_linear_positive_slope_fraction"].median()
            ),
            "n_bins_signed_linear_slope_ci_below_zero": int(
                group["signed_linear_slope_ci_below_zero"].astype(bool).sum()
            ),
            "n_bins_signed_linear_slope_ci_above_zero": int(
                group["signed_linear_slope_ci_above_zero"].astype(bool).sum()
            ),
            "median_incremental_linear_slope": float(
                group["incremental_linear_slope_median"].median()
            ),
            "median_power_alpha": float(group["power_alpha_median"].median()),
            "median_saturating_tau": float(group["saturating_tau_median"].median()),
            "n_bins_signed_linear_negative_trend_supported": int(
                classifications.eq("signed_linear_negative_trend_supported").sum()
            ),
            "n_bins_signed_linear_positive_trend_supported": int(
                classifications.eq("signed_linear_positive_trend_supported").sum()
            ),
            "n_bins_negative_slope_supported": int(
                signed_trends.eq("negative_slope_supported").sum()
            ),
            "n_bins_positive_slope_supported": int(
                signed_trends.eq("positive_slope_supported").sum()
            ),
            "n_bins_signed_slope_inconclusive": int(
                signed_trends.eq("signed_slope_inconclusive").sum()
            ),
            "n_bins_constant_level_preferred": int(
                classifications.eq("constant_level_preferred").sum()
            ),
            "n_bins_positive_incremental_accumulation_supported": int(
                classifications.eq("positive_incremental_accumulation_supported").sum()
            ),
            "n_bins_zero_origin_linear_scaling_supported": int(
                classifications.eq("zero_origin_linear_scaling_supported").sum()
            ),
            "n_bins_subdiffusive_power_law_supported": int(
                classifications.eq("subdiffusive_power_law_supported").sum()
            ),
            "n_bins_exponential_saturation_supported": int(
                classifications.eq("exponential_saturation_supported").sum()
            ),
            "n_bins_inconclusive_or_unstable": int(
                classifications.isin(["inconclusive", "fit_unstable"]).sum()
            ),
            "model_preference_threshold": float(args.model_preference_threshold),
            "signed_slope_probability_threshold": float(
                args.signed_slope_probability_threshold
            ),
            "subdiffusion_probability_threshold": float(
                args.subdiffusion_probability_threshold
            ),
            "parameter_bound_warning_fraction": float(
                args.parameter_bound_warning_fraction
            ),
        }

        half = max(1, int(math.ceil(n_bins / 2.0)))
        if row["n_bins_signed_linear_negative_trend_supported"] >= half:
            recommendation = "retain_signed_negative_temporal_trend"
            wording = (
                "A signed linear model with negative slope is robustly preferred across "
                "at least half of tested bins, supporting a decrease in dispersion over "
                "the observed 1--5 week lag range rather than positive accumulation."
            )
        elif row["n_bins_signed_linear_positive_trend_supported"] >= half:
            recommendation = "retain_signed_positive_temporal_trend"
            wording = (
                "A signed linear model with positive slope is robustly preferred across "
                "at least half of tested bins, supporting increasing dispersion over the "
                "observed 1--5 week lag range."
            )
        elif (
            row["n_bins_negative_slope_supported"] >= half
            and row["n_bins_constant_level_preferred"] >= half
        ):
            recommendation = "retain_negative_temporal_tendency_with_constant_aicc_preference"
            wording = (
                "The constant model remains AICc-preferred across at least half of tested "
                "bins, but the subject-bootstrap signed slope is negative with the predefined "
                "sign probability in at least half of bins. Report this as a reproducible "
                "negative temporal tendency without claiming that a two-parameter linear "
                "model is preferred over the constant model."
            )
        elif (
            row["n_bins_positive_slope_supported"] >= half
            and row["n_bins_constant_level_preferred"] >= half
        ):
            recommendation = "retain_positive_temporal_tendency_with_constant_aicc_preference"
            wording = (
                "The constant model remains AICc-preferred across at least half of tested "
                "bins, but the subject-bootstrap signed slope is positive with the predefined "
                "sign probability in at least half of bins. Report this as a positive temporal "
                "tendency without claiming signed-linear model preference."
            )
        elif row["n_bins_constant_level_preferred"] >= half:
            recommendation = "retain_no_detectable_temporal_trend_result"
            wording = (
                "The constant model is preferred across at least half of tested bins after "
                "including a signed-linear alternative, and no directional signed-slope "
                "pattern reaches the predefined criterion across at least half of bins. "
                "Interpret this as no detectable systematic temporal trend over the observed "
                "1--5 week lag range."
            )
        elif row["n_bins_positive_incremental_accumulation_supported"] >= half:
            recommendation = "retain_positive_incremental_accumulation_result"
            wording = (
                "Additional approximately linear accumulation above the one-week level "
                "is supported across at least half of tested bins."
            )
        elif row["n_bins_subdiffusive_power_law_supported"] >= half:
            recommendation = "retain_formal_subdiffusive_power_law_result"
            wording = (
                "Power-law subdiffusion is supported across at least half of tested "
                "bins and alpha is not boundary-dominated."
            )
        elif row["n_bins_exponential_saturation_supported"] >= half:
            recommendation = "retain_exponential_saturation_as_second_moment_result"
            wording = (
                "Exponential saturation is supported across at least half of tested "
                "bins and tau is not boundary-dominated; this remains a second-moment "
                "description rather than proof of an OU mechanism."
            )
        else:
            recommendation = "do_not_promote_mechanistic_temporal_model"
            wording = (
                "Evidence is heterogeneous; avoid assigning a unique temporal law."
            )

        row["decision_aid"] = recommendation
        row["recommended_wording"] = wording
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["estimator", "metric"]
    ).reset_index(drop=True)


# =============================================================================
# Main
# =============================================================================


def validate_args(args: argparse.Namespace) -> None:
    args.dt_values_parsed = parse_csv_ints(args.dt_values)
    args.metrics_parsed = parse_csv_strings(args.metrics)
    args.estimators_parsed = parse_csv_strings(args.estimators)
    args.bin_ids_parsed = parse_csv_ints(args.bin_ids) if args.bin_ids else []

    if len(args.dt_values_parsed) < args.min_dt_points:
        raise ValueError(
            "Selected {} dt values but --min-dt-points={}.".format(
                len(args.dt_values_parsed), args.min_dt_points
            )
        )
    unknown_metrics = sorted(set(args.metrics_parsed) - set(SUPPORTED_METRICS))
    unknown_estimators = sorted(set(args.estimators_parsed) - set(SUPPORTED_ESTIMATORS))
    if unknown_metrics:
        raise ValueError("Unsupported metrics: {}".format(unknown_metrics))
    if unknown_estimators:
        raise ValueError("Unsupported estimators: {}".format(unknown_estimators))
    if str(args.required_mode).upper() != PRIMARY_MODE and not args.allow_nonprimary:
        raise ValueError(
            "The primary Step-11.2 analysis requires upstream TT mode. "
            "Use --allow-nonprimary only for an explicit sensitivity run."
        )
    if args.n_boot <= 0:
        raise ValueError("--n-bootstrap must be > 0")
    if not (0 < args.core_min_subject_fraction <= 1):
        raise ValueError("--core-min-subject-fraction must be in (0,1]")
    if args.core_min_subjects is not None and args.core_min_subjects < 2:
        raise ValueError("--core-min-subjects must be >=2")
    if not (0 < args.ci < 100):
        raise ValueError("--ci must be between 0 and 100")
    if args.alpha_min <= 0 or args.alpha_max <= args.alpha_min:
        raise ValueError("Invalid alpha bounds")
    if args.tau_min <= 0 or args.tau_max <= args.tau_min:
        raise ValueError("Invalid tau bounds")
    if not (0 <= args.model_preference_threshold <= 1):
        raise ValueError("--model-preference-threshold must be in [0,1]")
    if not (0 <= args.signed_slope_probability_threshold <= 1):
        raise ValueError("--signed-slope-probability-threshold must be in [0,1]")
    if not (0 <= args.subdiffusion_probability_threshold <= 1):
        raise ValueError("--subdiffusion-probability-threshold must be in [0,1]")
    if not (0 <= args.parameter_bound_warning_fraction <= 1):
        raise ValueError("--parameter-bound-warning-fraction must be in [0,1]")


def main() -> int:
    args = build_argparser().parse_args()
    validate_args(args)

    step11_dir = Path(args.step11_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not step11_dir.exists():
        raise FileNotFoundError("Step-11 analysis directory not found: {}".format(step11_dir))

    metadata = load_metadata(step11_dir)
    validate_upstream_metadata(metadata, args.dt_values_parsed, args.required_mode)
    dataset_label = str(args.dataset_label).strip()
    if not dataset_label:
        dataset_label = str(metadata.get("dataset_label") or step11_dir.name)

    subject_path = (
        Path(args.subject_stats).expanduser()
        if str(args.subject_stats).strip()
        else find_table(step11_dir, "09_subject_bin_dynamics")
    )
    curves_path = (
        Path(args.curves).expanduser()
        if str(args.curves).strip()
        else find_table(step11_dir, "10_binned_dynamics_long")
    )

    availability_path = None
    if str(args.subject_availability).strip():
        availability_path = Path(args.subject_availability).expanduser()
    else:
        try:
            availability_path = find_table(step11_dir, "08_subject_availability_by_dt")
        except FileNotFoundError:
            availability_path = None

    print("[READ] {}".format(subject_path))
    subject_stats = read_table(subject_path)
    print("[READ] {}".format(curves_path))
    curves = read_table(curves_path)
    subject_stats, curves = normalize_inputs(subject_stats, curves)

    availability_summary = None
    if availability_path is not None:
        print("[READ] {}".format(availability_path))
        availability_raw = read_table(availability_path)
        availability_summary = normalize_subject_availability(
            availability_raw, args.dt_values_parsed
        )
    else:
        warnings.warn(
            "08_subject_availability_by_dt was not found; structural availability "
            "will be derived from 09_subject_bin_dynamics."
        )

    available_dts = sorted(set(subject_stats["dt"].astype(int)))
    missing_dts = sorted(set(args.dt_values_parsed) - set(available_dts))
    if missing_dts:
        raise ValueError(
            "Requested dt values are absent from subject-bin statistics: {}. Available: {}".format(
                missing_dts, available_dts
            )
        )

    all_subjects = sorted(subject_stats["subject"].astype(str).unique().tolist())
    if len(all_subjects) < 3:
        raise ValueError("At least 3 subjects are required; found {}".format(len(all_subjects)))

    subject_requirements = build_subject_requirements_by_dt(
        subject_stats,
        args.dt_values_parsed,
        args.core_min_subjects,
        args.core_min_subject_fraction,
        availability_summary=availability_summary,
    )
    args.core_min_subjects_resolved = int(
        subject_requirements["n_subjects_required"].min()
    )

    tested_bins, bin_summary = select_tested_bins(
        curves,
        args.dt_values_parsed,
        subject_requirements,
        args.bin_scope,
        args.bin_ids_parsed,
    )
    if not tested_bins:
        write_table(
            bin_summary,
            out_dir / "02_tested_abundance_bins",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )
        write_table(
            subject_requirements,
            out_dir / "08_subject_availability_and_requirements_by_dt",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )
        print("[DESIGN] Per-dt structural availability and requirements:")
        print(subject_requirements.to_string(index=False))
        if len(bin_summary):
            show_cols = [
                "bin_id", "bin_left", "bin_right",
                "min_subject_coverage_fraction_across_dt",
                "present_in_all_selected_dt", "meets_subject_requirement", "qualifies",
            ]
            print("[DESIGN] Best-covered candidate bins:")
            print(
                bin_summary.sort_values(
                    ["min_subject_coverage_fraction_across_dt", "bin_id"],
                    ascending=[False, True],
                )[show_cols].head(15).to_string(index=False)
            )
        raise ValueError(
            "No abundance bins satisfy the requested joint-dt criteria. "
            "Diagnostic tables were written to 02_tested_abundance_bins and "
            "08_subject_availability_and_requirements_by_dt."
        )

    print("[DESIGN] Per-dt structural availability and requirements:")
    print(subject_requirements.to_string(index=False))
    print(
        "[DESIGN] subjects={} dt={} bins={} estimators={} metrics={} bootstrap={}".format(
            len(all_subjects),
            args.dt_values_parsed,
            tested_bins,
            args.estimators_parsed,
            args.metrics_parsed,
            args.n_boot,
        )
    )

    arrays = prepare_subject_arrays(
        subject_stats,
        curves,
        all_subjects,
        args.dt_values_parsed,
        tested_bins,
    )

    observed_long, observed_pairs = observed_fits(
        arrays,
        all_subjects,
        args.dt_values_parsed,
        args.estimators_parsed,
        args.metrics_parsed,
        dataset_label,
        args,
    )
    bootstrap = joint_subject_bootstrap(
        arrays,
        all_subjects,
        args.dt_values_parsed,
        args.estimators_parsed,
        args.metrics_parsed,
        dataset_label,
        args,
    )
    by_bin = summarize_by_bin(observed_pairs, bootstrap, args)
    summary = summarize_overall(by_bin, args)
    interpretation = summary[
        [
            "dataset_label",
            "estimator",
            "metric",
            "n_bins_tested",
            "decision_aid",
            "recommended_wording",
        ]
    ].copy()

    manifest = pd.DataFrame(
        [
            {
                "dataset_label": dataset_label,
                "step11_dir": str(step11_dir),
                "subject_stats_path": str(subject_path),
                "curves_path": str(curves_path),
                "subject_availability_path": str(availability_path) if availability_path is not None else "",
                "n_subject_rows": int(len(subject_stats)),
                "n_curve_rows": int(len(curves)),
                "n_subjects": int(len(all_subjects)),
                "dt_values": ",".join(str(v) for v in args.dt_values_parsed),
                "tested_bin_ids": ",".join(str(v) for v in tested_bins),
            }
        ]
    )

    write_table(manifest, out_dir / "01_input_manifest", args.write_csv, args.write_parquet, args.compression)
    write_table(bin_summary, out_dir / "02_tested_abundance_bins", args.write_csv, args.write_parquet, args.compression)
    write_table(
        subject_requirements,
        out_dir / "08_subject_availability_and_requirements_by_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )
    write_table(observed_long, out_dir / "03_observed_six_model_fits_long", args.write_csv, args.write_parquet, args.compression)
    write_table(bootstrap, out_dir / "04_joint_subject_bootstrap_six_model_comparison", args.write_csv, args.write_parquet, args.compression)
    write_table(by_bin, out_dir / "05_model_comparison_by_bin", args.write_csv, args.write_parquet, args.compression)
    write_table(summary, out_dir / "06_model_comparison_summary", args.write_csv, args.write_parquet, args.compression)
    write_table(interpretation, out_dir / "07_interpretation_aid", args.write_csv, args.write_parquet, args.compression)

    config = vars(args).copy()
    for key in ["dt_values_parsed", "metrics_parsed", "estimators_parsed", "bin_ids_parsed"]:
        config[key] = list(config[key])
    config.update(
        {
            "script_version": SCRIPT_VERSION,
            "step": 13,
            "script": "13-temporal_fluctuation_scaling.py",
            "analysis_layer": "step13_temporal_fluctuation_scaling",
            "upstream_step": 11,
            "pseudo_baseline_step12_used": False,
            "generated": datetime.now().isoformat(timespec="seconds"),
            "primary_metric": PRIMARY_METRIC,
            "complementary_metric": COMPLEMENTARY_METRIC,
            "primary_estimator": "transition_weighted",
            "sensitivity_estimator": "equal_subject_weighted",
            "abundance_domain": "Step-11 internally comparable longitudinal core",
            "pseudo_overlap_used_for_temporal_scaling": False,
            "dataset_label_resolved": dataset_label,
            "n_subjects": len(all_subjects),
            "all_subjects": all_subjects,
            "tested_bins": tested_bins,
            "bootstrap_unit": "joint_subject_cluster",
            "same_subject_sample_across_all_dt_bins_metrics_estimators": True,
            "constant_model": "K (empirical over observed dt range only)",
            "constant_model_dt0_extrapolation_allowed": False,
            "linear_model": "D * dt (zero-origin Brownian comparator)",
            "signed_linear_model": "K + D * (dt - 1), K>=0, D signed, fitted values constrained >=0 over observed dt",
            "incremental_linear_model": "K + D * (dt - 1), K>=0, D>=0",
            "power_law_model": "C * dt**alpha",
            "saturating_model": "A * (1 - exp(-dt/tau))",
            "ou_mechanism_tested": False,
            "ou_second_moment_compatibility_only": True,
            "aicc_parameter_count_includes_residual_variance": True,
            "aicc_caveat": "Only five temporal lags; model discrimination is interpreted through joint subject-bootstrap stability, parameter identifiability and AICc together.",
            "nonlinear_model_promotion_requires_parameter_identifiability": True,
            "upstream_representation_validated": True,
            "upstream_required_representation": PRIMARY_REPRESENTATION,
            "upstream_required_conditioning": PRIMARY_CONDITIONING,
            "upstream_required_x_condition_source": "xstar_latent",
            "upstream_required_mode": str(args.required_mode).upper(),
        }
    )
    with open(out_dir / "00_run_config.json", "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)

    with open(out_dir / "DEBUG_REPORT.txt", "w", encoding="utf-8") as handle:
        handle.write("STEP 13 — TEMPORAL FLUCTUATION SCALING\n")
        handle.write("==========================================================\n\n")
        handle.write("Dataset: {}\n".format(dataset_label))
        handle.write("Subjects: {}\n".format(len(all_subjects)))
        handle.write("Temporal lags: {}\n".format(args.dt_values_parsed))
        handle.write("Required upstream mode: {}\n".format(str(args.required_mode).upper()))
        handle.write("Tested bins: {}\n".format(tested_bins))
        handle.write("Bootstrap replicates: {}\n".format(args.n_boot))
        handle.write("Same subject resample across all dt/bins/metrics/estimators: yes\n")
        handle.write("Subject requirements by dt:\n{}\n\n".format(subject_requirements.to_string(index=False)))
        handle.write("Compared models:\n")
        handle.write("  constant weekly-resolution level: K\n")
        handle.write("  signed linear: K + D * (dt - 1), D signed\n")
        handle.write("  positive incremental linear: K + D * (dt - 1), D>=0\n")
        handle.write("  zero-origin linear comparator: D * dt\n")
        handle.write("  power law: C * dt**alpha\n")
        handle.write("  exponential saturation: A * (1 - exp(-dt/tau))\n\n")
        handle.write("Interpretation boundary:\n")
        handle.write("  The constant model is empirical over the observed 1--5 week range and must not be extrapolated to dt=0.\n")
        handle.write("  The signed-linear model tests positive or negative trend; the positive incremental model tests accumulation only.\n")
        handle.write("  AIC/AICc parameter counts include the residual-variance parameter.\n")
        handle.write("  Nonlinear-model support requires both bootstrap model preference and non-boundary-dominated alpha/tau estimates.\n")
        handle.write("  The saturating second moment is OU-like compatible, but the restoring drift required for an OU process is not tested here.\n\n")
        handle.write(summary.to_string(index=False))
        handle.write("\n")

    print("[DONE] {}".format(out_dir))
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise