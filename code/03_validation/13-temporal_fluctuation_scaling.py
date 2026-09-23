#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13-temporal_fluctuation_scaling.py
===================================

Canonical ClonoDynamics Step 13.

PURPOSE
-------
Test whether the genuine longitudinal shared fluctuation component identified in
Step 11 shows systematic temporal accumulation over lags of 1–5 weeks.

The primary upstream fluctuation estimand is:

    cross_cov
        = Cov(dx_observed_rep1, dx_observed_rep2
              | xmid_latent, dt, common4)

Step 13 does NOT subtract the pseudo technical-null baseline from Step 12.
Step 12 establishes that the shared covariance is above the technical null;
Step 13 asks a different question:

    does the genuine longitudinal covariance increase with elapsed lag?

PRIMARY DESIGN
--------------
Temporal scaling is evaluated only after fixing both abundance support and
subject composition.

1. A fixed longitudinal abundance core is selected from the Step-11
   xmid_latent grid.

2. A bin qualifies when, at every requested lag, it contains at least

       ceil(core_min_subject_fraction × subjects structurally available at dt)

   biological subjects with a finite subject-level cross covariance.

3. The largest contiguous qualifying run is the primary abundance core.

4. The primary complete-case cohort contains only subjects with a valid
   PRIMARY `cross_cov` value in EVERY selected abundance bin at EVERY requested
   lag.

Technical comparator availability is then checked on this frozen primary
cohort/core and is never allowed to redefine the biological subject set.

Default:

    core_min_subject_fraction = 0.90

Because the six-week longitudinal design contains only nine subjects at dt=5,
the default 90% rule generally requires nine subjects at every lag.

PRIMARY TEMPORAL ESTIMATOR
--------------------------
For each complete-case subject and lag, the primary core summary is the
EQUAL-BIN mean across the fixed abundance core.

The cohort curve is then the EQUAL-SUBJECT mean across complete-case subjects.

Thus the primary temporal curve does not change because of lag-specific
abundance composition or subject composition.

Primary estimator label:

    equal_bin_equal_subject

SENSITIVITY ESTIMATOR
---------------------
A precision-oriented sensitivity weights abundance bins within each
subject × lag by covariance degrees of freedom:

    weight = n_subject_bin - 1

and then averages subjects equally.

Sensitivity estimator label:

    df_weighted_within_subject_equal_subject

This weighting is a sensitivity only. It does not replace the fixed-geometry
equal-bin primary analysis.

METRICS
-------
Primary:
    cross_cov

Technical comparators:
    same_var_mean
    replicate_specific_excess

The default analysis does not use shared_fraction as a primary temporal-scaling
metric because ratios can become unstable when the denominator is small.

SIGNED TEMPORAL SLOPE
---------------------
For each metric and estimator:

    M(dt) = intercept_at_dt1 + slope × (dt - 1)

with dt = 1,2,3,4,5.

The signed slope is the principal temporal descriptor.

Biological uncertainty is estimated by a JOINT subject-cluster bootstrap:
subjects are sampled once with replacement per bootstrap replicate and the same
subject draw is propagated across all lags, metrics and estimators.

Reported quantities include:

    observed slope
    bootstrap median slope
    2.5th–97.5th bootstrap interval
    bootstrap fraction slope > 0
    bootstrap fraction slope < 0

The preferred primary wording when positive accumulation is not supported is:

    "no systematic positive accumulation over the observed temporal lags"

rather than claiming complete temporal independence.

EXACT SUBJECT SIGN-FLIP TEST
----------------------------
For each metric and estimator, a subject-specific temporal slope is calculated.
An exact paired sign-flip test is then applied to the mean subject slope.

The script reports:

    two-sided exact P
    one-sided negative-trend P
    one-sided positive-trend P

This is retained as a conservative subject-level inferential complement to the
bootstrap interval.

FINITE-LAG MODEL DESCRIPTORS
----------------------------
Three deliberately simple finite-lag descriptions are compared on the observed
cohort curve:

    constant
        M(dt) = K

    signed_linear
        M(dt) = K + D × (dt - 1), D signed

    positive_incremental
        M(dt) = K + D × (dt - 1), D >= 0

AICc includes the residual-variance parameter.

With only five temporal lags, AICc model preference is descriptive and is kept
separate from signed-slope evidence. In particular, an AICc preference for the
constant model must not be interpreted automatically as proof of zero slope.

ABUNDANCE-RESOLVED SECONDARY ANALYSIS
-------------------------------------
Within the selected core, Step 13 additionally reports one temporal slope per
abundance bin for the primary cross_cov metric, using the same complete-case
subjects and the same joint subject-bootstrap draws.

This localizes temporal structure without assigning a mechanistic law to
individual abundance bins.

INPUT
-----
--step11-dir must contain final Step-11 outputs:

    00_run_config.json
    04_subject_fluctuation_by_bin_dt.csv
    06_support_by_dt.csv

The Step-5 transition table is not reloaded.

Step 12 is not a computational input.

OUTPUT
------
<outdir>/
    00_run_config.json
    00_analysis_summary.csv
    00_pipeline_manifest.csv
    README_outputs.md

    01_core_bin_selection.csv
    02_complete_case_subjects.csv
    03_subject_core_metric_by_dt.csv
    04_cohort_core_metric_by_dt.csv
    05_temporal_slope_summary.csv
    06_subject_slopes.csv
    07_exact_signflip_tests.csv
    08_binwise_temporal_slopes.csv
    09_model_comparison.csv
    10_joint_subject_bootstrap.csv
    11_step14_contract.csv
    12_core_cell_support_audit.csv

    00_run_signature.json
    00_input_manifest.csv

PRIMARY RUN
-----------
python3 13-temporal_fluctuation_scaling.py \
    --step11-dir \
    ./results_4/11-cross_replicate_fluctuation_dynamics \
    --outdir \
    ./results_4/13-temporal_fluctuation_scaling \
    --dt-values 1 2 3 4 5 \
    --core-min-subject-fraction 0.90 \
    --min-subject-cell-n 2 \
    --n-bootstrap 2000 \
    --bootstrap-seed 123

INTERPRETATION BOUNDARY
-----------------------
Step 13 analyzes genuine longitudinal temporal scaling only.

It does NOT:
    - subtract the pseudo Step-12 covariance;
    - test calendar-position composition (Step 14);
    - impose TT/TF/FT/FF operational classes;
    - test observation-threshold sensitivity (Steps 15–16);
    - infer a specific stochastic mechanism from five lag values.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v3-final-step11v2-fixed-primary-cohort-2026-09-17"

PRIMARY_METRIC = "cross_cov"
DEFAULT_METRICS = (
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
)

PRIMARY_ESTIMATOR = "equal_bin_equal_subject"
SENSITIVITY_ESTIMATOR = "df_weighted_within_subject_equal_subject"
ESTIMATORS = (
    PRIMARY_ESTIMATOR,
    SENSITIVITY_ESTIMATOR,
)

MODELS = (
    "constant",
    "signed_linear",
    "positive_incremental",
)


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


def read_json(path: Path) -> Dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
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
    pd.DataFrame(rows).to_csv(
        outdir / "00_pipeline_manifest.csv",
        index=False,
    )


def largest_contiguous_run(values: Iterable[int]) -> List[int]:
    values = sorted(set(int(x) for x in values))
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

    # Longest run; if tied, choose the lower-abundance-starting run.
    runs.sort(
        key=lambda run: (len(run), -run[0]),
        reverse=True,
    )
    return runs[0]


def parse_int_list(values: Sequence[int]) -> List[int]:
    return sorted(set(int(x) for x in values))


def finite_numeric(values) -> np.ndarray:
    """Return finite numeric values from any array-like input."""
    x = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    ).to_numpy(dtype=float)
    return x[np.isfinite(x)]


def safe_linear_fit(
    dt_values: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    dt_values = np.asarray(dt_values, dtype=float)
    y = np.asarray(y, dtype=float)

    ok = np.isfinite(dt_values) & np.isfinite(y)

    if int(ok.sum()) < 2:
        return {
            "n_dt": int(ok.sum()),
            "intercept_at_dt1": np.nan,
            "slope_per_week": np.nan,
        }

    x = dt_values[ok] - 1.0
    values = y[ok]

    if np.ptp(x) <= 0:
        return {
            "n_dt": int(ok.sum()),
            "intercept_at_dt1": np.nan,
            "slope_per_week": np.nan,
        }

    slope, intercept = np.polyfit(
        x,
        values,
        1,
    )

    return {
        "n_dt": int(ok.sum()),
        "intercept_at_dt1": float(intercept),
        "slope_per_week": float(slope),
    }


def percentile_interval(
    values: np.ndarray,
    low: float = 0.025,
    high: float = 0.975,
) -> Tuple[float, float, float, int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan, np.nan, np.nan, 0

    return (
        float(np.quantile(x, low)),
        float(np.median(x)),
        float(np.quantile(x, high)),
        int(x.size),
    )


# =============================================================================
# Input validation
# =============================================================================

def validate_step11_contract(
    config: Mapping[str, object],
    signature: Mapping[str, object],
) -> None:
    if str(config.get("primary_support")) != "common4":
        raise ValueError(
            "Step 13 requires Step 11 primary_support='common4'."
        )

    if str(config.get("primary_conditioning")) != "xmid_latent":
        raise ValueError(
            "Step 13 requires Step 11 primary_conditioning='xmid_latent'."
        )

    if config.get("operational_class_filter") is not None:
        raise ValueError(
            "Step 13 requires Step 11 without TT/TF/FT/FF primary filtering."
        )

    estimands = config.get("primary_estimands", {})
    if isinstance(estimands, dict):
        trunc = estimands.get("cross_cov_truncated_at_zero")
        if trunc is True:
            raise ValueError(
                "Step 13 requires signed, untruncated Step-11 cross covariance."
            )


    payload = signature.get("signature_payload", {})
    if payload.get("primary_support") != "common4":
        raise ValueError(
            "Step-11 run signature does not certify common4 support."
        )
    if payload.get("primary_conditioning") != "xmid_latent":
        raise ValueError(
            "Step-11 run signature does not certify xmid_latent conditioning."
        )
    if payload.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError(
            "Step-11 run signature does not certify signed cross_cov."
        )

    bootstrap = payload.get("bootstrap", {})
    if (
        isinstance(bootstrap, dict)
        and bootstrap.get(
            "bootstrap_cell_min_n_matches_point_estimate"
        ) is not True
    ):
        raise ValueError(
            "Step 13 requires final Step-11 bootstrap support contract."
        )


# =============================================================================
# Core selection
# =============================================================================

def normalize_subject_table(
    subject_table: pd.DataFrame,
    metrics: Sequence[str],
) -> pd.DataFrame:
    required = [
        "subject",
        "dt",
        "bin",
        "n",
        "x_left",
        "x_right",
        "x_center",
        *metrics,
    ]
    require_columns(
        subject_table,
        required,
        "Step-11 subject fluctuation table",
    )

    d = subject_table.copy()

    numeric_cols = [
        "subject",
        "dt",
        "bin",
        "n",
        "x_left",
        "x_right",
        "x_center",
        *metrics,
    ]

    for c in numeric_cols:
        d[c] = pd.to_numeric(
            d[c],
            errors="coerce",
        )

    d = d.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    d = d.dropna(
        subset=[
            "subject",
            "dt",
            "bin",
            "n",
            "x_center",
        ]
    ).copy()

    d["subject"] = d["subject"].astype(int)
    d["dt"] = d["dt"].astype(int)
    d["bin"] = d["bin"].astype(int)

    # There should already be one row per subject × dt × bin.
    duplicated = d.duplicated(
        subset=["subject", "dt", "bin"],
        keep=False,
    )
    if duplicated.any():
        bad = d.loc[
            duplicated,
            ["subject", "dt", "bin"],
        ].drop_duplicates()
        raise ValueError(
            "Step-11 subject fluctuation table contains duplicate "
            f"subject × dt × bin rows. Examples:\n{bad.head(10)}"
        )

    return d


def build_core_selection(
    subject_table: pd.DataFrame,
    support: pd.DataFrame,
    dt_values: Sequence[int],
    core_min_subject_fraction: float,
    min_subject_cell_n: int,
    explicit_bins: Sequence[int],
) -> Tuple[pd.DataFrame, List[int], Dict[int, int]]:
    require_columns(
        support,
        ["dt", "n_subjects"],
        "Step-11 support table",
    )

    support = support.copy()
    support["dt"] = pd.to_numeric(
        support["dt"],
        errors="coerce",
    )
    support["n_subjects"] = pd.to_numeric(
        support["n_subjects"],
        errors="coerce",
    )

    availability = {}

    for dt in dt_values:
        row = support[
            support["dt"].eq(int(dt))
        ]

        if row.empty:
            raise ValueError(
                f"Step-11 support table is missing dt={dt}."
            )

        n_available = int(
            row.iloc[0]["n_subjects"]
        )
        if n_available < 2:
            raise ValueError(
                f"dt={dt}: fewer than two structurally available subjects."
            )

        availability[int(dt)] = n_available

    requirements = {
        dt: max(
            2,
            int(
                math.ceil(
                    float(core_min_subject_fraction)
                    * availability[dt]
                )
            ),
        )
        for dt in dt_values
    }

    # A cell is structurally valid for temporal scaling when covariance is
    # defined and at least two transitions contribute. min_subject_cell_n can
    # be raised explicitly as a robustness choice.
    valid = subject_table[
        subject_table["dt"].isin(
            [int(x) for x in dt_values]
        )
        & (
            subject_table["n"]
            >= int(min_subject_cell_n)
        )
        & np.isfinite(
            subject_table[PRIMARY_METRIC]
        )
    ].copy()

    bins = sorted(
        int(x)
        for x in subject_table["bin"].unique()
    )

    reference = (
        subject_table[
            ["bin", "x_left", "x_right", "x_center"]
        ]
        .drop_duplicates("bin")
        .set_index("bin")
    )

    rows = []
    qualifying = []

    for bin_id in bins:
        ref = reference.loc[int(bin_id)]

        row: Dict[str, object] = {
            "bin": int(bin_id),
            "x_left": float(ref["x_left"]),
            "x_right": float(ref["x_right"]),
            "x_center": float(ref["x_center"]),
        }

        qualifies = True

        for dt in dt_values:
            g = valid[
                valid["dt"].eq(int(dt))
                & valid["bin"].eq(int(bin_id))
            ]

            n_subjects = int(
                g["subject"].nunique()
            )

            row[
                f"n_subjects_valid_dt{dt}"
            ] = n_subjects
            row[
                f"n_subjects_available_dt{dt}"
            ] = int(
                availability[int(dt)]
            )
            row[
                f"n_subjects_required_dt{dt}"
            ] = int(
                requirements[int(dt)]
            )
            row[
                f"coverage_fraction_dt{dt}"
            ] = (
                n_subjects
                / float(availability[int(dt)])
            )

            meets = (
                n_subjects
                >= requirements[int(dt)]
            )
            row[
                f"meets_requirement_dt{dt}"
            ] = bool(meets)

            qualifies = (
                qualifies and meets
            )

        row["qualifies"] = bool(qualifies)

        if qualifies:
            qualifying.append(
                int(bin_id)
            )

        rows.append(row)

    summary = pd.DataFrame(rows)

    if explicit_bins:
        explicit = sorted(
            set(int(x) for x in explicit_bins)
        )

        if explicit != list(
            range(explicit[0], explicit[-1] + 1)
        ):
            raise ValueError(
                "Explicit --core-bins must form one contiguous abundance core."
            )

        missing = sorted(
            set(explicit)
            - set(qualifying)
        )

        if missing:
            raise ValueError(
                "Explicit bins fail the temporal-core subject-coverage "
                f"requirements: {missing}"
            )

        selected = explicit
    else:
        selected = largest_contiguous_run(
            qualifying
        )

    if not selected:
        raise ValueError(
            "No contiguous abundance core satisfies the requested coverage rule."
        )

    summary[
        "selected_for_primary_core"
    ] = summary["bin"].isin(
        selected
    )

    return (
        summary.sort_values("bin"),
        selected,
        requirements,
    )


def complete_case_subjects(
    subject_table: pd.DataFrame,
    dt_values: Sequence[int],
    core_bins: Sequence[int],
    metrics: Sequence[str],
    min_subject_cell_n: int,
) -> Tuple[List[int], pd.DataFrame]:
    subjects = sorted(
        int(x)
        for x in subject_table["subject"].unique()
    )

    rows = []
    complete = []

    expected_cells = (
        len(dt_values)
        * len(core_bins)
    )

    for subject in subjects:
        d = subject_table[
            subject_table["subject"].eq(int(subject))
            & subject_table["dt"].isin(
                [int(x) for x in dt_values]
            )
            & subject_table["bin"].isin(
                [int(x) for x in core_bins]
            )
        ].copy()

        cell_valid = (
            d["n"].ge(
                int(min_subject_cell_n)
            )
        )

        for metric in metrics:
            cell_valid &= np.isfinite(
                pd.to_numeric(
                    d[metric],
                    errors="coerce",
                )
            )

        n_valid = int(
            cell_valid.sum()
        )

        is_complete = (
            len(d) == expected_cells
            and n_valid == expected_cells
        )

        rows.append(
            {
                "subject": int(subject),
                "n_expected_cells": int(expected_cells),
                "n_present_cells": int(len(d)),
                "n_valid_cells": int(n_valid),
                "complete_case": bool(
                    is_complete
                ),
            }
        )

        if is_complete:
            complete.append(
                int(subject)
            )

    return (
        complete,
        pd.DataFrame(rows),
    )


def validate_metrics_on_fixed_primary_cohort(
    subject_table: pd.DataFrame,
    complete_subjects: Sequence[int],
    dt_values: Sequence[int],
    core_bins: Sequence[int],
    metrics: Sequence[str],
    min_subject_cell_n: int,
) -> None:
    """
    The primary cohort/core is defined by cross_cov only.

    Every requested technical comparator must be finite on that exact frozen
    support. Missing comparator values cause an explicit error rather than
    silently changing the primary subject composition.
    """
    d = subject_table[
        subject_table["subject"].isin(
            [int(x) for x in complete_subjects]
        )
        & subject_table["dt"].isin(
            [int(x) for x in dt_values]
        )
        & subject_table["bin"].isin(
            [int(x) for x in core_bins]
        )
    ].copy()

    expected = (
        len(complete_subjects)
        * len(dt_values)
        * len(core_bins)
    )

    if len(d) != expected:
        raise RuntimeError(
            "Frozen primary cohort/core does not contain the expected number "
            "of subject x dt x bin cells."
        )

    if not (
        pd.to_numeric(
            d["n"], errors="coerce"
        ) >= int(min_subject_cell_n)
    ).all():
        raise RuntimeError(
            "Frozen primary cohort/core contains an under-supported subject cell."
        )

    for metric in metrics:
        values = pd.to_numeric(
            d[metric],
            errors="coerce",
        ).to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(
                f"Metric {metric!r} is incomplete on the frozen primary "
                "cross_cov cohort/core. The primary cohort will not be changed."
            )


def core_cell_support_audit(
    subject_table: pd.DataFrame,
    complete_subjects: Sequence[int],
    dt_values: Sequence[int],
    core_bins: Sequence[int],
) -> pd.DataFrame:
    d = subject_table[
        subject_table["subject"].isin(
            [int(x) for x in complete_subjects]
        )
        & subject_table["dt"].isin(
            [int(x) for x in dt_values]
        )
        & subject_table["bin"].isin(
            [int(x) for x in core_bins]
        )
    ].copy()

    rows = []
    for dt in dt_values:
        g = d[d["dt"].eq(int(dt))].copy()
        n = pd.to_numeric(
            g["n"], errors="coerce"
        ).to_numpy(float)

        rows.append(
            {
                "dt": int(dt),
                "n_complete_case_subjects": int(
                    g["subject"].nunique()
                ),
                "n_core_bins": int(
                    g["bin"].nunique()
                ),
                "n_subject_bin_cells": int(len(g)),
                "min_subject_bin_n": float(np.nanmin(n)),
                "q025_subject_bin_n": float(
                    np.nanquantile(n, 0.025)
                ),
                "median_subject_bin_n": float(
                    np.nanmedian(n)
                ),
                "q975_subject_bin_n": float(
                    np.nanquantile(n, 0.975)
                ),
                "max_subject_bin_n": float(np.nanmax(n)),
            }
        )

    return pd.DataFrame(rows)


# =============================================================================
# Core summaries
# =============================================================================

def build_subject_core_table(
    subject_table: pd.DataFrame,
    complete_subjects: Sequence[int],
    dt_values: Sequence[int],
    core_bins: Sequence[int],
    metrics: Sequence[str],
) -> pd.DataFrame:
    d = subject_table[
        subject_table["subject"].isin(
            [int(x) for x in complete_subjects]
        )
        & subject_table["dt"].isin(
            [int(x) for x in dt_values]
        )
        & subject_table["bin"].isin(
            [int(x) for x in core_bins]
        )
    ].copy()

    rows = []

    for (subject, dt), g in d.groupby(
        ["subject", "dt"],
        sort=True,
    ):
        if len(g) != len(core_bins):
            raise RuntimeError(
                f"Complete-case contract broken for subject={subject}, dt={dt}."
            )

        for estimator in ESTIMATORS:
            if estimator == PRIMARY_ESTIMATOR:
                weights = np.ones(
                    len(g),
                    dtype=float,
                )
            elif estimator == SENSITIVITY_ESTIMATOR:
                weights = np.maximum(
                    pd.to_numeric(
                        g["n"],
                        errors="coerce",
                    ).to_numpy(float)
                    - 1.0,
                    1.0,
                )
            else:
                raise ValueError(
                    estimator
                )

            for metric in metrics:
                values = pd.to_numeric(
                    g[metric],
                    errors="coerce",
                ).to_numpy(float)

                valid = (
                    np.isfinite(values)
                    & np.isfinite(weights)
                    & (weights > 0)
                )

                if not valid.all():
                    raise RuntimeError(
                        "Non-finite complete-case core metric encountered."
                    )

                value = float(
                    np.average(
                        values,
                        weights=weights,
                    )
                )

                rows.append(
                    {
                        "subject": int(subject),
                        "dt": int(dt),
                        "estimator": estimator,
                        "metric": metric,
                        "value": value,
                        "n_core_bins": int(
                            len(g)
                        ),
                        "sum_subject_bin_n": int(
                            pd.to_numeric(
                                g["n"],
                                errors="coerce",
                            ).sum()
                        ),
                        "sum_covariance_df": float(
                            np.maximum(
                                pd.to_numeric(
                                    g["n"],
                                    errors="coerce",
                                ).to_numpy(float)
                                - 1.0,
                                0.0,
                            ).sum()
                        ),
                    }
                )

    return pd.DataFrame(rows).sort_values(
        [
            "estimator",
            "metric",
            "subject",
            "dt",
        ]
    )


def cohort_table_from_subject_core(
    subject_core: pd.DataFrame,
) -> pd.DataFrame:
    return (
        subject_core.groupby(
            ["estimator", "metric", "dt"],
            as_index=False,
            sort=True,
        )
        .agg(
            value=("value", "mean"),
            subject_sd=("value", "std"),
            n_subjects=("subject", "nunique"),
        )
    )


# =============================================================================
# Models / AICc
# =============================================================================

def aicc_from_prediction(
    y: np.ndarray,
    predicted: np.ndarray,
    structural_parameters: int,
) -> Dict[str, float]:
    y = np.asarray(
        y,
        dtype=float,
    )
    predicted = np.asarray(
        predicted,
        dtype=float,
    )

    valid = (
        np.isfinite(y)
        & np.isfinite(predicted)
    )

    y = y[valid]
    predicted = predicted[valid]

    n = len(y)

    if n == 0:
        return {
            "rss": np.nan,
            "aic": np.nan,
            "aicc": np.nan,
        }

    rss = float(
        np.sum(
            (y - predicted) ** 2
        )
    )
    safe_rss = max(
        rss,
        np.finfo(float).tiny,
    )

    # residual variance counts as one fitted parameter
    k = int(structural_parameters) + 1

    aic = float(
        n * np.log(
            safe_rss / n
        )
        + 2.0 * k
    )

    if n > k + 1:
        aicc = float(
            aic
            + (
                2.0 * k * (k + 1)
            )
            / (
                n - k - 1.0
            )
        )
    else:
        aicc = np.nan

    return {
        "rss": rss,
        "aic": aic,
        "aicc": aicc,
    }


def fit_temporal_models(
    dt_values: np.ndarray,
    y: np.ndarray,
) -> Dict[str, Dict[str, float]]:
    dt_values = np.asarray(
        dt_values,
        dtype=float,
    )
    y = np.asarray(
        y,
        dtype=float,
    )

    valid = (
        np.isfinite(dt_values)
        & np.isfinite(y)
    )

    x = (
        dt_values[valid] - 1.0
    )
    values = y[valid]

    if len(values) < 3:
        return {}

    # Constant.
    K = float(
        np.mean(values)
    )
    pred_constant = np.full(
        len(values),
        K,
        dtype=float,
    )

    constant = {
        "intercept_at_dt1": K,
        "slope_per_week": 0.0,
        **aicc_from_prediction(
            values,
            pred_constant,
            structural_parameters=1,
        ),
    }

    # Signed linear.
    slope, intercept = np.polyfit(
        x,
        values,
        1,
    )
    pred_signed = (
        intercept
        + slope * x
    )

    signed = {
        "intercept_at_dt1": float(intercept),
        "slope_per_week": float(slope),
        **aicc_from_prediction(
            values,
            pred_signed,
            structural_parameters=2,
        ),
    }

    # Positive incremental model.
    # If unconstrained slope is negative, the constrained optimum is D=0 and
    # K=mean(values). The model still carries two structural parameters for
    # AICc accounting.
    if slope >= 0:
        inc_slope = float(slope)
        # The declared constraint is D >= 0 only. K is unconstrained.
        inc_intercept = float(intercept)
        pred_incremental = (
            inc_intercept
            + inc_slope * x
        )
    else:
        inc_slope = 0.0
        inc_intercept = K
        pred_incremental = (
            pred_constant.copy()
        )

    incremental = {
        "intercept_at_dt1": float(
            inc_intercept
        ),
        "slope_per_week": float(
            inc_slope
        ),
        **aicc_from_prediction(
            values,
            pred_incremental,
            structural_parameters=2,
        ),
    }

    return {
        "constant": constant,
        "signed_linear": signed,
        "positive_incremental": incremental,
    }


# =============================================================================
# Bootstrap and subject slopes
# =============================================================================

def prepare_subject_arrays(
    subject_core: pd.DataFrame,
    complete_subjects: Sequence[int],
    dt_values: Sequence[int],
    metrics: Sequence[str],
) -> Dict[Tuple[str, str], np.ndarray]:
    arrays: Dict[
        Tuple[str, str],
        np.ndarray
    ] = {}

    for estimator in ESTIMATORS:
        for metric in metrics:
            matrix = np.full(
                (
                    len(complete_subjects),
                    len(dt_values),
                ),
                np.nan,
                dtype=float,
            )

            for i, subject in enumerate(
                complete_subjects
            ):
                g = subject_core[
                    subject_core["subject"].eq(
                        int(subject)
                    )
                    & subject_core["estimator"].eq(
                        estimator
                    )
                    & subject_core["metric"].eq(
                        metric
                    )
                ]

                value_map = dict(
                    zip(
                        g["dt"].astype(int),
                        g["value"].astype(float),
                    )
                )

                for j, dt in enumerate(
                    dt_values
                ):
                    matrix[i, j] = value_map[
                        int(dt)
                    ]

            if not np.isfinite(
                matrix
            ).all():
                raise RuntimeError(
                    f"Non-finite subject core matrix for {estimator}/{metric}."
                )

            arrays[
                (estimator, metric)
            ] = matrix

    return arrays


def joint_subject_bootstrap(
    arrays: Mapping[
        Tuple[str, str],
        np.ndarray
    ],
    complete_subjects: Sequence[int],
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
) -> Tuple[pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(
        int(seed)
    )

    n_subjects = len(
        complete_subjects
    )

    # Store the actual subject indices so the exact same draws can also be
    # propagated into the abundance-resolved binwise sensitivity.
    draws = rng.integers(
        0,
        n_subjects,
        size=(
            int(n_bootstrap),
            n_subjects,
        ),
    )

    rows = []

    dt_array = np.asarray(
        dt_values,
        dtype=float,
    )

    for boot_id in range(
        int(n_bootstrap)
    ):
        idx = draws[boot_id]

        for (
            estimator,
            metric,
        ), matrix in arrays.items():

            curve = matrix[
                idx,
                :
            ].mean(
                axis=0
            )

            linear = safe_linear_fit(
                dt_array,
                curve,
            )

            models = fit_temporal_models(
                dt_array,
                curve,
            )

            preferred_model = ""
            finite_aicc = {
                model: fit["aicc"]
                for model, fit in models.items()
                if np.isfinite(
                    fit["aicc"]
                )
            }

            if finite_aicc:
                preferred_model = min(
                    finite_aicc,
                    key=finite_aicc.get,
                )

            row = {
                "bootstrap": int(
                    boot_id
                ),
                "estimator": estimator,
                "metric": metric,
                "slope_per_week": linear[
                    "slope_per_week"
                ],
                "intercept_at_dt1": linear[
                    "intercept_at_dt1"
                ],
                "preferred_model": preferred_model,
                "n_unique_subjects_sampled": int(
                    len(
                        np.unique(
                            idx
                        )
                    )
                ),
            }

            for j, dt in enumerate(
                dt_values
            ):
                row[
                    f"value_dt{dt}"
                ] = float(
                    curve[j]
                )

            for model, fit in models.items():
                row[
                    f"{model}_aicc"
                ] = fit["aicc"]

            rows.append(row)

    return (
        pd.DataFrame(rows),
        draws,
    )


def build_subject_slopes(
    subject_core: pd.DataFrame,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    rows = []

    dt_array = np.asarray(
        dt_values,
        dtype=float,
    )

    for (
        estimator,
        metric,
        subject,
    ), g in subject_core.groupby(
        [
            "estimator",
            "metric",
            "subject",
        ],
        sort=True,
    ):
        g = g.sort_values("dt")

        if list(
            g["dt"].astype(int)
        ) != [
            int(x)
            for x in dt_values
        ]:
            raise RuntimeError(
                "Subject slope table is not complete across dt."
            )

        fit = safe_linear_fit(
            dt_array,
            g["value"].to_numpy(float),
        )

        rows.append(
            {
                "subject": int(subject),
                "estimator": estimator,
                "metric": metric,
                "n_dt": int(
                    fit["n_dt"]
                ),
                "intercept_at_dt1": fit[
                    "intercept_at_dt1"
                ],
                "slope_per_week": fit[
                    "slope_per_week"
                ],
            }
        )

    return pd.DataFrame(rows).sort_values(
        [
            "estimator",
            "metric",
            "subject",
        ]
    )


def exact_signflip(
    slopes: np.ndarray,
) -> Dict[str, float]:
    slopes = np.asarray(
        slopes,
        dtype=float,
    )
    slopes = slopes[
        np.isfinite(slopes)
    ]

    n = len(slopes)

    if n == 0:
        return {
            "n_subjects": 0,
            "mean_subject_slope": np.nan,
            "two_sided_p": np.nan,
            "one_sided_negative_p": np.nan,
            "one_sided_positive_p": np.nan,
        }

    if n > 20:
        raise ValueError(
            "Exact sign-flip enumeration is restricted to <=20 subjects."
        )

    observed = float(
        np.mean(slopes)
    )

    null_values = np.empty(
        2 ** n,
        dtype=float,
    )

    for mask in range(
        2 ** n
    ):
        signs = np.array(
            [
                1.0
                if (
                    mask >> i
                ) & 1
                else -1.0
                for i in range(n)
            ],
            dtype=float,
        )

        null_values[mask] = float(
            np.mean(
                signs * slopes
            )
        )

    eps = 1e-15

    two_sided = float(
        np.mean(
            np.abs(
                null_values
            )
            >= abs(
                observed
            ) - eps
        )
    )

    negative = float(
        np.mean(
            null_values
            <= observed + eps
        )
    )

    positive = float(
        np.mean(
            null_values
            >= observed - eps
        )
    )

    return {
        "n_subjects": int(n),
        "mean_subject_slope": observed,
        "two_sided_p": two_sided,
        "one_sided_negative_p": negative,
        "one_sided_positive_p": positive,
        "n_negative_subject_slopes": int(
            np.sum(
                slopes < 0
            )
        ),
        "n_positive_subject_slopes": int(
            np.sum(
                slopes > 0
            )
        ),
        "n_zero_subject_slopes": int(
            np.sum(
                slopes == 0
            )
        ),
    }


# =============================================================================
# Binwise temporal slopes
# =============================================================================

def binwise_primary_slopes(
    subject_table: pd.DataFrame,
    complete_subjects: Sequence[int],
    core_bins: Sequence[int],
    dt_values: Sequence[int],
    bootstrap_draws: np.ndarray,
) -> pd.DataFrame:
    d = subject_table[
        subject_table["subject"].isin(
            [int(x) for x in complete_subjects]
        )
        & subject_table["bin"].isin(
            [int(x) for x in core_bins]
        )
        & subject_table["dt"].isin(
            [int(x) for x in dt_values]
        )
    ].copy()

    subject_index = {
        int(subject): i
        for i, subject in enumerate(
            complete_subjects
        )
    }

    dt_array = np.asarray(
        dt_values,
        dtype=float,
    )

    rows = []

    for bin_id in core_bins:
        g = d[
            d["bin"].eq(
                int(bin_id)
            )
        ]

        matrix = np.full(
            (
                len(complete_subjects),
                len(dt_values),
            ),
            np.nan,
            dtype=float,
        )

        for row in g.itertuples(
            index=False
        ):
            i = subject_index[
                int(row.subject)
            ]
            j = list(
                dt_values
            ).index(
                int(row.dt)
            )
            matrix[
                i,
                j,
            ] = float(
                row.cross_cov
            )

        if not np.isfinite(
            matrix
        ).all():
            raise RuntimeError(
                f"Core bin {bin_id} is not complete in the complete-case cohort."
            )

        observed_curve = matrix.mean(
            axis=0
        )
        observed_fit = safe_linear_fit(
            dt_array,
            observed_curve,
        )

        boot_slopes = np.empty(
            bootstrap_draws.shape[0],
            dtype=float,
        )

        for b in range(
            bootstrap_draws.shape[0]
        ):
            curve = matrix[
                bootstrap_draws[b],
                :
            ].mean(
                axis=0
            )

            boot_slopes[b] = safe_linear_fit(
                dt_array,
                curve,
            )[
                "slope_per_week"
            ]

        lo, med, hi, n_valid = (
            percentile_interval(
                boot_slopes
            )
        )

        ref = g.iloc[0]

        rows.append(
            {
                "bin": int(bin_id),
                "x_left": float(
                    ref["x_left"]
                ),
                "x_right": float(
                    ref["x_right"]
                ),
                "x_center": float(
                    ref["x_center"]
                ),
                "observed_slope_per_week": observed_fit[
                    "slope_per_week"
                ],
                "bootstrap_median_slope": med,
                "bootstrap_q025": lo,
                "bootstrap_q975": hi,
                "n_bootstrap_valid": n_valid,
                "bootstrap_positive_fraction": float(
                    np.mean(
                        boot_slopes[
                            np.isfinite(
                                boot_slopes
                            )
                        ] > 0
                    )
                ),
                "bootstrap_negative_fraction": float(
                    np.mean(
                        boot_slopes[
                            np.isfinite(
                                boot_slopes
                            )
                        ] < 0
                    )
                ),
            }
        )

    return pd.DataFrame(rows).sort_values(
        "bin"
    )


# =============================================================================
# Summaries
# =============================================================================

def build_slope_summary(
    cohort: pd.DataFrame,
    bootstrap: pd.DataFrame,
    subject_slopes: pd.DataFrame,
    dt_values: Sequence[int],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    slope_rows = []
    signflip_rows = []
    model_rows = []

    dt_array = np.asarray(
        dt_values,
        dtype=float,
    )

    for estimator in ESTIMATORS:
        for metric in sorted(
            cohort["metric"].unique()
        ):
            d = cohort[
                cohort["estimator"].eq(
                    estimator
                )
                & cohort["metric"].eq(
                    metric
                )
            ].sort_values(
                "dt"
            )

            values = d[
                "value"
            ].to_numpy(float)

            fit = safe_linear_fit(
                dt_array,
                values,
            )

            b = bootstrap[
                bootstrap["estimator"].eq(
                    estimator
                )
                & bootstrap["metric"].eq(
                    metric
                )
            ]

            slopes = finite_numeric(
                b["slope_per_week"]
            )

            lo, med, hi, n_valid = (
                percentile_interval(
                    slopes
                )
            )

            slope_rows.append(
                {
                    "estimator": estimator,
                    "metric": metric,
                    "n_complete_case_subjects": int(
                        d["n_subjects"].min()
                    ),
                    "n_core_bins": None,
                    "observed_intercept_at_dt1": fit[
                        "intercept_at_dt1"
                    ],
                    "observed_slope_per_week": fit[
                        "slope_per_week"
                    ],
                    "bootstrap_median_slope": med,
                    "bootstrap_q025": lo,
                    "bootstrap_q975": hi,
                    "n_bootstrap_valid": n_valid,
                    "bootstrap_positive_slope_fraction": (
                        float(
                            np.mean(
                                slopes > 0
                            )
                        )
                        if slopes.size
                        else np.nan
                    ),
                    "bootstrap_negative_slope_fraction": (
                        float(
                            np.mean(
                                slopes < 0
                            )
                        )
                        if slopes.size
                        else np.nan
                    ),
                    "bootstrap_ci_excludes_positive": bool(
                        np.isfinite(hi)
                        and hi < 0
                    ),
                    "bootstrap_ci_excludes_negative": bool(
                        np.isfinite(lo)
                        and lo > 0
                    ),
                }
            )

            subj = subject_slopes[
                subject_slopes["estimator"].eq(
                    estimator
                )
                & subject_slopes["metric"].eq(
                    metric
                )
            ]

            test = exact_signflip(
                subj[
                    "slope_per_week"
                ].to_numpy(float)
            )

            signflip_rows.append(
                {
                    "estimator": estimator,
                    "metric": metric,
                    **test,
                }
            )

            models = fit_temporal_models(
                dt_array,
                values,
            )

            finite = {
                model: fit_info[
                    "aicc"
                ]
                for model, fit_info in models.items()
                if np.isfinite(
                    fit_info["aicc"]
                )
            }

            preferred = (
                min(
                    finite,
                    key=finite.get,
                )
                if finite
                else ""
            )

            for model, fit_info in models.items():
                bb = b[
                    [
                        "preferred_model",
                        f"{model}_aicc",
                    ]
                ].copy()

                aicc_values = finite_numeric(
                    bb[
                        f"{model}_aicc"
                    ]
                )

                model_rows.append(
                    {
                        "estimator": estimator,
                        "metric": metric,
                        "model": model,
                        "observed_intercept_at_dt1": fit_info[
                            "intercept_at_dt1"
                        ],
                        "observed_slope_per_week": fit_info[
                            "slope_per_week"
                        ],
                        "observed_rss": fit_info[
                            "rss"
                        ],
                        "observed_aic": fit_info[
                            "aic"
                        ],
                        "observed_aicc": fit_info[
                            "aicc"
                        ],
                        "observed_preferred_model": preferred,
                        "bootstrap_preference_fraction": float(
                            np.mean(
                                b[
                                    "preferred_model"
                                ].astype(str)
                                .eq(model)
                            )
                        ),
                        "bootstrap_median_aicc": (
                            float(
                                np.median(
                                    aicc_values
                                )
                            )
                            if aicc_values.size
                            else np.nan
                        ),
                    }
                )

    model_df = pd.DataFrame(model_rows)

    if not model_df.empty:
        min_aicc = (
            model_df.groupby(
                ["estimator", "metric"]
            )["observed_aicc"]
            .transform("min")
        )
        model_df["delta_aicc"] = (
            model_df["observed_aicc"]
            - min_aicc
        )

    return (
        pd.DataFrame(slope_rows),
        pd.DataFrame(signflip_rows),
        model_df,
    )


def add_bootstrap_ci_to_cohort(
    cohort: pd.DataFrame,
    bootstrap: pd.DataFrame,
    dt_values: Sequence[int],
) -> pd.DataFrame:
    rows = []

    for row in cohort.itertuples(
        index=False
    ):
        b = bootstrap[
            bootstrap["estimator"].eq(
                row.estimator
            )
            & bootstrap["metric"].eq(
                row.metric
            )
        ]

        values = finite_numeric(
            b[
                f"value_dt{int(row.dt)}"
            ]
        )

        lo, med, hi, n_valid = (
            percentile_interval(
                values
            )
        )

        record = dict(
            row._asdict()
        )
        record.update(
            {
                "bootstrap_median": med,
                "bootstrap_q025": lo,
                "bootstrap_q975": hi,
                "n_bootstrap_valid": n_valid,
            }
        )
        rows.append(record)

    return pd.DataFrame(rows)


# =============================================================================
# README / CLI
# =============================================================================

def write_readme(outdir: Path) -> None:
    text = """# Step 13 — temporal scaling of genuine cross-replicate fluctuations

Primary metric
--------------
cross_cov = Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt, common4)

Primary abundance domain
------------------------
Largest contiguous fixed xmid_latent core satisfying the lag-specific
subject-coverage rule.

Primary cohort
--------------
Complete-case biological subjects are defined by the PRIMARY cross_cov metric
in every selected core bin at every requested lag. Technical comparators are
evaluated on this frozen biological cohort and cannot redefine it.

Primary estimator
-----------------
equal_bin_equal_subject:
1. equal weight to every abundance bin within subject × lag;
2. equal weight to every complete-case biological subject.

Sensitivity
-----------
df_weighted_within_subject_equal_subject:
abundance bins are weighted by n_subject_bin - 1 within subject × lag, then
subjects are averaged equally.

Primary temporal claim
----------------------
The signed slope across dt = 1–5 is the principal temporal descriptor.

Subject-cluster bootstrap uses the same subject resample across all lags,
metrics and estimators.

Exact subject-level sign-flip tests are reported separately.

Model comparison
----------------
constant, signed-linear and positive-incremental finite-lag descriptions are
compared by AICc. The positive-incremental model constrains only D >= 0; its
intercept K remains unconstrained. With only five lag values, AICc is
descriptive and does not override signed-slope evidence.

Pipeline boundary
-----------------
Step 12 is not subtracted or otherwise used computationally.

Step 14 consumes the selected core bins and complete-case subject cohort from
11_step14_contract.csv together with Step-11 interval-resolved outputs.
"""
    (
        outdir / "README_outputs.md"
    ).write_text(
        text,
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 13: fixed-core complete-case temporal scaling of genuine "
            "cross-replicate fluctuation covariance."
        ),
    )

    p.add_argument(
        "--step11-dir",
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
        "--metrics",
        nargs="+",
        default=list(
            DEFAULT_METRICS
        ),
        choices=[
            "cross_cov",
            "same_var_mean",
            "replicate_specific_excess",
            "shared_fraction",
        ],
    )

    p.add_argument(
        "--core-min-subject-fraction",
        type=float,
        default=0.90,
    )

    p.add_argument(
        "--min-subject-cell-n",
        type=int,
        default=2,
        help=(
            "Minimum transitions required for a subject x dt x bin covariance "
            "cell to be considered structurally valid. The default 2 is the "
            "mathematical minimum for ddof=1 covariance."
        ),
    )

    p.add_argument(
        "--core-bins",
        nargs="*",
        type=int,
        default=None,
        help=(
            "Optional explicit Step-11 bin IDs. When omitted, the largest "
            "contiguous qualifying core is selected automatically."
        ),
    )

    p.add_argument(
        "--min-complete-case-subjects",
        type=int,
        default=3,
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

    return p


def validate_args(
    args: argparse.Namespace,
) -> None:
    dt_values = parse_int_list(
        args.dt_values
    )

    if len(dt_values) < 3:
        raise ValueError(
            "Step 13 requires at least three temporal lags."
        )

    if any(
        dt < 1
        for dt in dt_values
    ):
        raise ValueError(
            "--dt-values must all be >=1"
        )

    if not (
        0 < float(
            args.core_min_subject_fraction
        ) <= 1
    ):
        raise ValueError(
            "--core-min-subject-fraction must be in (0,1]"
        )

    if int(
        args.min_subject_cell_n
    ) < 2:
        raise ValueError(
            "--min-subject-cell-n must be >=2"
        )

    if int(
        args.min_complete_case_subjects
    ) < 2:
        raise ValueError(
            "--min-complete-case-subjects must be >=2"
        )

    if int(
        args.n_bootstrap
    ) < 1:
        raise ValueError(
            "--n-bootstrap must be >=1"
        )


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(
        argv
    )
    validate_args(args)

    dt_values = parse_int_list(
        args.dt_values
    )

    metrics = list(
        dict.fromkeys(
            args.metrics
        )
    )

    if PRIMARY_METRIC not in metrics:
        raise ValueError(
            "cross_cov must be included in --metrics because it defines "
            "the primary temporal core."
        )

    step11_dir = (
        args.step11_dir
        .expanduser()
        .resolve(strict=True)
    )

    outdir = ensure_dir(
        args.outdir
        .expanduser()
        .resolve()
    )

    config_path = (
        step11_dir
        / "00_run_config.json"
    )
    signature_path = (
        step11_dir
        / "00_run_signature.json"
    )
    subject_path = (
        step11_dir
        / "04_subject_fluctuation_by_bin_dt.csv"
    )
    support_path = (
        step11_dir
        / "06_support_by_dt.csv"
    )

    step11_config = read_json(
        config_path
    )
    step11_signature = read_json(
        signature_path
    )
    validate_step11_contract(
        step11_config,
        step11_signature,
    )

    subject_table = normalize_subject_table(
        read_csv_required(
            subject_path,
            "Step-11 subject fluctuation table",
        ),
        metrics,
    )

    support = read_csv_required(
        support_path,
        "Step-11 support table",
    )

    available_dts = sorted(
        int(x)
        for x in subject_table[
            "dt"
        ].unique()
    )

    missing_dts = sorted(
        set(dt_values)
        - set(available_dts)
    )

    if missing_dts:
        raise ValueError(
            f"Requested temporal lags missing from Step 11: {missing_dts}"
        )

    core_summary, core_bins, requirements = (
        build_core_selection(
            subject_table,
            support,
            dt_values,
            core_min_subject_fraction=float(
                args.core_min_subject_fraction
            ),
            min_subject_cell_n=int(
                args.min_subject_cell_n
            ),
            explicit_bins=(
                args.core_bins
                if args.core_bins
                else []
            ),
        )
    )

    complete_subjects, subject_completeness = (
        complete_case_subjects(
            subject_table,
            dt_values,
            core_bins,
            [PRIMARY_METRIC],
            min_subject_cell_n=int(
                args.min_subject_cell_n
            ),
        )
    )

    if len(
        complete_subjects
    ) < int(
        args.min_complete_case_subjects
    ):
        raise ValueError(
            "Too few complete-case subjects: "
            f"{len(complete_subjects)} < {args.min_complete_case_subjects}"
        )

    validate_metrics_on_fixed_primary_cohort(
        subject_table,
        complete_subjects,
        dt_values,
        core_bins,
        metrics,
        min_subject_cell_n=int(
            args.min_subject_cell_n
        ),
    )

    core_support_audit = core_cell_support_audit(
        subject_table,
        complete_subjects,
        dt_values,
        core_bins,
    )

    subject_completeness[
        "selected_for_primary_complete_case"
    ] = subject_completeness[
        "subject"
    ].isin(
        complete_subjects
    )

    subject_core = build_subject_core_table(
        subject_table,
        complete_subjects,
        dt_values,
        core_bins,
        metrics,
    )

    cohort = cohort_table_from_subject_core(
        subject_core
    )

    arrays = prepare_subject_arrays(
        subject_core,
        complete_subjects,
        dt_values,
        metrics,
    )

    bootstrap, bootstrap_draws = (
        joint_subject_bootstrap(
            arrays,
            complete_subjects,
            dt_values,
            n_bootstrap=int(
                args.n_bootstrap
            ),
            seed=int(
                args.bootstrap_seed
            ),
        )
    )

    cohort = add_bootstrap_ci_to_cohort(
        cohort,
        bootstrap,
        dt_values,
    )

    subject_slopes = build_subject_slopes(
        subject_core,
        dt_values,
    )

    slope_summary, signflip, model_comparison = (
        build_slope_summary(
            cohort,
            bootstrap,
            subject_slopes,
            dt_values,
        )
    )

    slope_summary[
        "n_core_bins"
    ] = int(
        len(core_bins)
    )

    slope_summary[
        "core_x_min"
    ] = float(
        core_summary.loc[
            core_summary[
                "selected_for_primary_core"
            ],
            "x_left",
        ].min()
    )

    slope_summary[
        "core_x_max"
    ] = float(
        core_summary.loc[
            core_summary[
                "selected_for_primary_core"
            ],
            "x_right",
        ].max()
    )

    binwise = binwise_primary_slopes(
        subject_table,
        complete_subjects,
        core_bins,
        dt_values,
        bootstrap_draws,
    )

    # Compact analysis summary.
    summary = slope_summary.merge(
        signflip[
            [
                "estimator",
                "metric",
                "mean_subject_slope",
                "two_sided_p",
                "one_sided_negative_p",
                "one_sided_positive_p",
                "n_negative_subject_slopes",
                "n_positive_subject_slopes",
            ]
        ],
        on=[
            "estimator",
            "metric",
        ],
        how="left",
        validate="one_to_one",
    )

    summary[
        "primary_analysis"
    ] = (
        summary["estimator"].eq(
            PRIMARY_ESTIMATOR
        )
        & summary["metric"].eq(
            PRIMARY_METRIC
        )
    )

    summary[
        "interpretation"
    ] = np.where(
        summary[
            "bootstrap_q975"
        ] < 0,
        "negative temporal slope; no positive accumulation",
        np.where(
            summary[
                "bootstrap_q025"
            ] > 0,
            "positive temporal accumulation",
            "no resolved signed temporal slope",
        ),
    )

    # Step-14 contract: one row per selected core bin.
    selected_core = core_summary[
        core_summary[
            "selected_for_primary_core"
        ]
    ].copy()

    complete_subject_string = ";".join(
        str(x)
        for x in complete_subjects
    )

    selected_core[
        "complete_case_subjects"
    ] = complete_subject_string
    selected_core[
        "n_complete_case_subjects"
    ] = int(
        len(complete_subjects)
    )
    selected_core[
        "dt_values"
    ] = ",".join(
        str(x)
        for x in dt_values
    )
    selected_core[
        "primary_metric"
    ] = PRIMARY_METRIC
    selected_core[
        "primary_estimator"
    ] = PRIMARY_ESTIMATOR

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "step11_analysis_signature": step11_config.get(
            "analysis_signature"
        ),
        "step11_subject_table": path_identity(subject_path),
        "step11_support_table": path_identity(support_path),
        "parameters": {
            "dt_values": [int(x) for x in dt_values],
            "metrics": list(metrics),
            "primary_metric": PRIMARY_METRIC,
            "primary_estimator": PRIMARY_ESTIMATOR,
            "sensitivity_estimator": SENSITIVITY_ESTIMATOR,
            "core_min_subject_fraction": float(
                args.core_min_subject_fraction
            ),
            "min_subject_cell_n": int(
                args.min_subject_cell_n
            ),
            "min_complete_case_subjects": int(
                args.min_complete_case_subjects
            ),
            "n_bootstrap": int(args.n_bootstrap),
            "bootstrap_seed": int(args.bootstrap_seed),
        },
        "selected_core_bins": [
            int(x) for x in core_bins
        ],
        "complete_case_subjects": [
            int(x) for x in complete_subjects
        ],
        "primary_complete_case_defined_by": PRIMARY_METRIC,
    }
    analysis_signature = stable_sha256(
        signature_payload
    )

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

    input_manifest = pd.DataFrame(
        [
            {
                "input": "step11_run_config",
                **path_identity(config_path),
            },
            {
                "input": "step11_run_signature",
                **path_identity(signature_path),
            },
            {
                "input": "step11_subject_fluctuation_by_bin_dt",
                **path_identity(subject_path),
            },
            {
                "input": "step11_support_by_dt",
                **path_identity(support_path),
            },
        ]
    )
    input_manifest.to_csv(
        outdir / "00_input_manifest.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Write outputs.
    # -----------------------------------------------------------------
    core_summary.to_csv(
        outdir
        / "01_core_bin_selection.csv",
        index=False,
    )

    subject_completeness.to_csv(
        outdir
        / "02_complete_case_subjects.csv",
        index=False,
    )

    subject_core.to_csv(
        outdir
        / "03_subject_core_metric_by_dt.csv",
        index=False,
    )

    cohort.to_csv(
        outdir
        / "04_cohort_core_metric_by_dt.csv",
        index=False,
    )

    slope_summary.to_csv(
        outdir
        / "05_temporal_slope_summary.csv",
        index=False,
    )

    subject_slopes.to_csv(
        outdir
        / "06_subject_slopes.csv",
        index=False,
    )

    signflip.to_csv(
        outdir
        / "07_exact_signflip_tests.csv",
        index=False,
    )

    binwise.to_csv(
        outdir
        / "08_binwise_temporal_slopes.csv",
        index=False,
    )

    model_comparison.to_csv(
        outdir
        / "09_model_comparison.csv",
        index=False,
    )

    bootstrap.to_csv(
        outdir
        / "10_joint_subject_bootstrap.csv",
        index=False,
    )

    selected_core[
        "step13_analysis_signature"
    ] = analysis_signature

    selected_core.to_csv(
        outdir
        / "11_step14_contract.csv",
        index=False,
    )

    core_support_audit.to_csv(
        outdir
        / "12_core_cell_support_audit.csv",
        index=False,
    )

    summary.to_csv(
        outdir
        / "00_analysis_summary.csv",
        index=False,
    )

    run_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_signature": analysis_signature,
        "step11_dir": str(
            step11_dir
        ),
        "outdir": str(
            outdir
        ),
        "upstream_step11_script_version":
            step11_config.get(
                "script_version"
            ),
        "upstream_step11_analysis_signature":
            step11_config.get(
                "analysis_signature"
            ),
        "primary_estimand": (
            "cross_cov = Cov(dx_observed_rep1, dx_observed_rep2 | "
            "xmid_latent, dt, common4)"
        ),
        "step12_used_computationally": False,
        "dt_values": dt_values,
        "metrics": metrics,
        "primary_metric": PRIMARY_METRIC,
        "primary_estimator": PRIMARY_ESTIMATOR,
        "sensitivity_estimator":
            SENSITIVITY_ESTIMATOR,
        "core_selection": {
            "rule": (
                "largest contiguous Step-11 bin run meeting the lag-specific "
                "subject-coverage requirement"
            ),
            "core_min_subject_fraction": float(
                args.core_min_subject_fraction
            ),
            "min_subject_cell_n": int(
                args.min_subject_cell_n
            ),
            "requirements_by_dt": {
                str(k): int(v)
                for k, v in requirements.items()
            },
            "selected_bins": [
                int(x)
                for x in core_bins
            ],
            "x_min": float(
                selected_core[
                    "x_left"
                ].min()
            ),
            "x_max": float(
                selected_core[
                    "x_right"
                ].max()
            ),
        },
        "complete_case": {
            "defined_by_metric": PRIMARY_METRIC,
            "technical_comparators_cannot_redefine_cohort": True,
            "n_subjects": int(
                len(complete_subjects)
            ),
            "subjects": [
                int(x)
                for x in complete_subjects
            ],
            "requires_all_selected_bins_at_all_selected_dt":
                True,
        },
        "bootstrap": {
            "unit": "biological subject",
            "n_bootstrap": int(
                args.n_bootstrap
            ),
            "seed": int(
                args.bootstrap_seed
            ),
            "same_subject_draw_across_dt_metrics_estimators":
                True,
        },
        "signed_slope": {
            "model": (
                "M(dt)=intercept_at_dt1+slope_per_week*(dt-1)"
            ),
            "primary_temporal_descriptor":
                True,
        },
        "model_comparison": {
            "models": list(
                MODELS
            ),
            "positive_incremental_constraint": "slope D >= 0 only; intercept K unconstrained",
            "aicc_includes_residual_variance":
                True,
            "interpretive_role":
                "descriptive finite-lag comparison; does not override signed-slope evidence",
        },
        "exact_signflip": {
            "unit": "subject-specific slope",
            "two_sided": True,
            "one_sided_negative": True,
            "one_sided_positive": True,
        },
    }

    (
        outdir
        / "00_run_config.json"
    ).write_text(
        json.dumps(
            json_safe(
                run_config
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_readme(
        outdir
    )
    write_manifest(
        outdir
    )

    # -----------------------------------------------------------------
    # Console summary.
    # -----------------------------------------------------------------
    print(
        "\n[DONE] Step 13 temporal fluctuation scaling"
    )
    print(
        "Selected core bins       :",
        core_bins,
    )
    print(
        "Core xmid range          :",
        f"{selected_core['x_left'].min():.6f}",
        "to",
        f"{selected_core['x_right'].max():.6f}",
    )
    print(
        "Complete-case subjects   :",
        complete_subjects,
    )

    primary = summary[
        summary[
            "primary_analysis"
        ]
    ]

    if not primary.empty:
        row = primary.iloc[0]

        print(
            "Primary cross_cov slope :",
            f"{row['observed_slope_per_week']:.6f}",
        )
        print(
            "Bootstrap 95% interval  :",
            f"[{row['bootstrap_q025']:.6f}, "
            f"{row['bootstrap_q975']:.6f}]",
        )
        print(
            "Bootstrap P(slope > 0)  :",
            f"{row['bootstrap_positive_slope_fraction']:.6f}",
        )
        print(
            "Exact sign-flip P (2s)  :",
            f"{row['two_sided_p']:.6f}",
        )

    print(
        "Outputs:"
    )
    for filename in [
        "00_analysis_summary.csv",
        "01_core_bin_selection.csv",
        "02_complete_case_subjects.csv",
        "03_subject_core_metric_by_dt.csv",
        "04_cohort_core_metric_by_dt.csv",
        "05_temporal_slope_summary.csv",
        "06_subject_slopes.csv",
        "07_exact_signflip_tests.csv",
        "08_binwise_temporal_slopes.csv",
        "09_model_comparison.csv",
        "10_joint_subject_bootstrap.csv",
        "11_step14_contract.csv",
        "12_core_cell_support_audit.csv",
        "00_run_signature.json",
        "00_input_manifest.csv",
    ]:
        print(
            " -",
            outdir / filename,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
