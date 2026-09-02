#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-pseudo_longitudinal_representation_assessment.py
===================================================

Pseudo-longitudinal displacement-representation assessment for ClonoDynamics.

Overview
--------
This script implements Step 8 of the ClonoDynamics pipeline.

It evaluates how alternative finite-time displacement representations behave
under the pseudo-longitudinal technical null, using the uncertainty-weighted
transition coordinate `xstar_latent` as the default conditioning coordinate.

The analysis is intentionally an ASSESSMENT rather than a ranking procedure.

It addresses three distinct questions:

    1. Directional-null adequacy

       Do alternative displacement representations preserve the near-null
       directional structure expected when no genuine longitudinal biological
       change is present?

    2. Residual sequencing-depth dependence

       Among relative-frequency representations that are directionally
       compatible with the pseudo-longitudinal null, how strongly does signed
       displacement remain coupled to sequencing-depth change?

    3. Closure / count decomposition

       How different are changes in observed relative frequency from changes in
       pooled log counts, and how strongly is that difference associated with
       sequencing-depth change?

Dispersion descriptors are also reported to document the scale and fluctuation
magnitude of each representation. Smaller dispersion is NOT interpreted as
superior performance because the candidate displacement representations are not
necessarily on equivalent stochastic scales.

No automatic representation winner is assigned.


Position in the ClonoDynamics pipeline
--------------------------------------
The relevant public workflow is:

    Step 1
    1-repertoire_characterization.py
                |
                v

    Step 2
    2-noise_aware_latent_inference.py
                |
                v

    Step 3
    3-latent_posterior_quality_characterization.py
                |
                v

    Step 4
    4-latent_trajectory_construction.py
                |
                v

    Step 5
    5-latent_transition_construction.py
                |
                v

    Step 6
    6-latent_transition_dataset_characterization.py

    Step 7
    7-pseudo_longitudinal_conditioning_benchmark.py

        technical-null evaluation of:
            conditioning geometry,
            replicate decoupling,
            pseudo-time ordering
                |
                v

    Step 8
    8-pseudo_longitudinal_representation_assessment.py

        xstar-conditioned comparison of:
            latent-frequency displacement,
            observed-frequency displacement,
            log-count displacement,
            closure/depth effects
                |
                v

        representation choice/justification for subsequent
        longitudinal fluctuation and dispersion analyses.

Step 7 determines the appropriate estimand-specific conditioning strategy.
Step 8 then keeps the conditioning coordinate fixed and evaluates the
DISPLACEMENT REPRESENTATION itself.


Scientific null
---------------
For a pseudo transition t0 -> t1, define a displacement representation r as

    Delta x_i^(r)
        = x_i^(r)(t1) - x_i^(r)(t0).

Because pseudo-time is constructed from technical measurements without genuine
biological longitudinal evolution, the primary directional-null expectations
are approximately:

    median(Delta x^(r)) = 0

and

    P(Delta x^(r) > 0) = 0.5.

The analysis therefore evaluates residual directional structure without
treating absolute displacement amplitude alone as a measure of inferential
quality.


Primary input
-------------
The required input is a pseudo-longitudinal transition table following the
Step-5 transition schema.

It is supplied with:

    --input

Typical input:

    latent_transitions.parquet

Supported input formats are:

    .parquet
    .pq
    .csv
    .tsv
    .txt.


Required transition identifiers
-------------------------------
The input must contain:

    subject
        Subject identifier.

    aaSeqCDR3
        Amino-acid CDR3 clonotype identifier.

    dt
        Temporal separation of the pseudo transition.

    obs_class
        Endpoint observability class (TT, TF, FT or FF).


Conditioning coordinate
-----------------------
The default conditioning coordinate is:

    xstar_latent

selected with:

    --conditioning-col xstar_latent.

For a transition with endpoint latent log-frequencies x0 and x1 and endpoint
posterior standard deviations sigma0 and sigma1, Step 5 defines:

    w0 = 1 / sigma0^2
    w1 = 1 / sigma1^2

and

                         w0*x0 + w1*x1
    xstar_latent = -----------------------------.
                              w0 + w1

When valid endpoint uncertainties are unavailable, Step 5 may use the midpoint
fallback.

In this Step, xstar is treated as a transition-centered conditioning coordinate
for displacement-magnitude/fluctuation benchmarking.

The conditioning coordinate is held fixed across all displacement
representations.


Default displacement representations
------------------------------------
The default mapping is:

    latent_frequency = dx_latent
    observed_frequency = dx_freq_obs
    log_count = dx_count_sum.

It can be changed with:

    --representations label=column,...

At least two representations are required.


1. Latent-frequency displacement
--------------------------------
The primary latent displacement is:

    Delta x_latent
        = x1_latent - x0_latent

and is represented by:

    dx_latent.

This quantity is based on latent clonotype-frequency estimates inferred from
the technical-replicate noise model.


2. Observed-frequency displacement
----------------------------------
The directly observed relative-frequency displacement is:

    Delta x_freq_obs
        = x1_obs - x0_obs

and is represented by:

    dx_freq_obs.

It retains the compositional normalization by sequencing depth and therefore
may remain coupled to changes in the observation denominator.


3. Log-count displacement
-------------------------
When pooled replicate counts are available, Step 5 defines:

    x_count_sum(t)
        = ln[count_sum(t) + 1]

and:

    Delta x_count
        = x_count_sum(t1) - x_count_sum(t0).

This is represented by:

    dx_count_sum.

Log-count displacement is retained as a NON-COMPOSITIONAL CONTROL.

It is not treated as a candidate relative-frequency representation.


Relative-frequency candidate policy
-----------------------------------
By default, the representations treated as relative-frequency candidates are:

    latent_frequency
    observed_frequency.

This is controlled by:

    --relative-frequency-representations.

Representations not included in this list remain available as controls.

The default `log_count` representation is therefore a control rather than a
candidate for replacing latent relative frequency.


Closure component
-----------------
When available, the transition table contains:

    closure_component
        = dx_freq_obs - dx_count_sum.

This quantity is supplied with:

    --closure-col

and defaults to:

    closure_component.

The closure component is analyzed separately.

It is NEVER treated as a candidate displacement representation.

It quantifies the difference between the observed relative-frequency change and
the pooled log-count change, thereby exposing the contribution associated with
relative-frequency normalization and sequencing-depth variation.


Sequencing-depth change
-----------------------
The default depth columns are:

    depth_sum_t0
    depth_sum_t1.

For valid positive depths D0 and D1, the script calculates:

    Delta log D
        = ln(D1) - ln(D0)

and:

    mean_log_depth
        = [ln(D0) + ln(D1)] / 2.

The primary descriptor of residual observation dependence is:

    |Spearman(Delta x, Delta log D)|.

Smaller absolute correlation indicates weaker residual coupling between signed
displacement and sequencing-depth change, but this criterion is interpreted
only together with directional-null compatibility.


Transition subset
-----------------
The default primary analysis uses:

    --dt-values 1
    --classes TT.

Thus, the primary pseudo benchmark is restricted to one-step transitions that
are observable at both endpoints.

Alternative dt values or observability classes can be supplied explicitly.


Common-row policy
-----------------
All displacement representations are evaluated on exactly the SAME transition
rows.

A row is retained only when:

    - the conditioning coordinate is finite;
    - every selected displacement representation is finite.

This complete-case policy prevents representation-specific missingness from
creating artificial differences between candidates.

Depth and closure variables are optional and do not determine inclusion in the
representation-comparison row set.


Fixed conditioning bins
-----------------------
The conditioning coordinate is divided into fixed abundance bins shared by all
representations.

Bin edges can be supplied explicitly using:

    --edges-csv.

Accepted edge-table formats contain either:

    x_lo, x_hi

or:

    bin_left, bin_right

or:

    edge.

If no edge table is provided, an abundance interval is derived from the
conditioning coordinate using:

    --range-q-low
    --range-q-high

with defaults:

    0.01
    0.99.

The interval is divided into:

    --n-bins

equal-width bins, with default:

    30.

The same rows and bin edges are then used for every displacement
representation.


Per-bin displacement statistics
-------------------------------
For representation r and abundance bin b, let:

    Delta x_{i,b}^{(r)}

denote the finite displacement values.

The script calculates:


1. Median displacement

    median_dx_b
        = median_i(Delta x_{i,b}).


2. Mean displacement

    mean_dx_b
        = mean_i(Delta x_{i,b}).


3. Positive-displacement probability

    ppos_b
        = P(Delta x_{i,b} > 0).


4. Mean squared displacement

    MSD_b
        = mean_i[(Delta x_{i,b})^2].


5. Displacement variance

    Var_b
        = Var_i(Delta x_{i,b}).

The implementation uses population variance (`ddof=0`).


6. Median absolute deviation

    MAD_b
        = median_i(
            |Delta x_{i,b} - median_dx_b|
          ).


7. Empirical displacement quantiles

    q025_dx_b
    q975_dx_b.

The script additionally stores:

    abs_median_dx_b
        = |median_dx_b|

and:

    ppos_bias_b
        = |ppos_b - 0.5|.


Common-bin requirement
----------------------
A pooled bin is considered sufficiently populated when:

    n_b >= --min-n

with default:

    200.

Only bin IDs meeting this requirement for EVERY selected representation are
included in representation-level summary metrics.

At least:

    --min-common-bins

common valid bins are required.

The default is:

    10.


Directional-null summary metrics
--------------------------------
Across the B common valid bins, the script calculates:


1. Mean absolute mean displacement

                                 1
    mean_abs_mean_dx = ---------------- sum_b |mean_dx_b|.
                                 B


2. Mean absolute median displacement

                                   1
    mean_abs_median_dx = ---------------- sum_b |median_dx_b|.
                                   B


3. RMS median displacement

                                  1
    rms_median_dx = sqrt( ---------------- sum_b median_dx_b^2 ).
                                  B


4. Maximum absolute median displacement

    max_abs_median_dx
        = max_b |median_dx_b|.


5. Median-displacement range

    median_dx_range
        = max_b(median_dx_b) - min_b(median_dx_b).


6. Mean absolute probability bias

                                  1
    mean_abs_ppos_bias = ---------------- sum_b |ppos_b - 0.5|.
                                  B


7. RMS probability bias

                                1
    rms_ppos_bias = sqrt( ---------------- sum_b (ppos_b - 0.5)^2 ).
                                B


8. Maximum probability bias

    max_abs_ppos_bias
        = max_b |ppos_b - 0.5|.

These metrics quantify directional-null adequacy.

They are descriptive and are not combined into a global representation score.


Residual depth-dependence metrics
---------------------------------
For each representation, the script calculates:


1. Spearman correlation

    rho_S
        = Spearman(Delta x, Delta log D).


2. Absolute Spearman correlation

    |rho_S|.

This is the primary empirical descriptor of residual observation dependence.


3. Pearson correlation

    r
        = Pearson(Delta x, Delta log D).


4. OLS slope

The script fits the univariate slope:

    Delta x = beta * Delta log D + intercept

and reports:

    beta
        = Cov(Delta log D, Delta x)
          / Var(Delta log D).


5. Displacement-magnitude depth dependence

The script also calculates:

    Spearman[(Delta x)^2, mean_log_depth].

This assesses whether displacement magnitude varies systematically with
sequencing depth even when signed displacement does not.


Dispersion descriptors
----------------------
The following representation-level quantities are retained:

    mean_msd
    median_msd

    mean_mad_dx
    median_mad_dx

    msd_log_roughness.

For the binwise MSD sequence M_b, log-roughness is defined from the absolute
second difference of log(M_b):

    roughness
        = mean |Delta^2 log(M_b)|.

These quantities characterize fluctuation magnitude and smoothness.

They are NOT representation-selection criteria.

In particular, a representation is not considered superior merely because it
has a smaller MSD or MAD.


Representation roles
--------------------
The assessment table explicitly labels each representation as either:

    relative_frequency_candidate

or:

    non_compositional_control.

It also records:

    directional_null_metrics_role
        = descriptive_adequacy_not_global_ranking

    depth_dependence_role
        = residual_observation_dependence_descriptor

    dispersion_metrics_role
        = descriptive_not_selection_criterion

and:

    automatic_winner_assigned
        = False.


Closure/depth diagnostic
------------------------
For the closure component C:

    C = dx_freq_obs - dx_count_sum,

the script evaluates:

    Pearson(C, Delta log D)

    Spearman(C, Delta log D)

and the OLS slope:

    C = beta_C * Delta log D + intercept.

It also reports closure behavior across the same xstar conditioning bins:

    median_closure
    mean_closure
    MAD_closure
    MSD_closure
    P(closure > 0).

These analyses are diagnostic and are kept separate from the representation
candidate assessment.


Pairwise displacement agreement
-------------------------------
For every pair of displacement representations A and B, evaluated on common
finite rows, the script reports:

    Pearson(A, B)
    Spearman(A, B)

    median(A - B)
    median|A - B|

    RMS(A - B).

This quantifies how similarly different representations describe individual
pseudo transitions.


Subject-level representation assessment
---------------------------------------
The same fixed pooled bin edges are reused separately within each subject.

Subject-level bins require:

    --min-n-subject

observations, with default:

    50.

At least:

    --min-common-bins-subject

common valid bins are required, with default:

    6.

The same directional, dispersion and depth-dependence descriptors are then
calculated at subject level.


Subject bootstrap
-----------------
Pairwise representation contrasts are evaluated through SUBJECT-level
bootstrap resampling.

The default number of bootstrap replicates is:

    --n-subject-bootstrap 2000.

For representations A and B and metric M, the subject-level difference is:

    Delta M_s
        = M_s(A) - M_s(B).

Subjects are resampled with replacement and the mean subject-level difference
is calculated for every bootstrap draw.

The output includes:

    subject_mean_delta_a_minus_b
    subject_median_delta_a_minus_b

    bootstrap_mean_delta_a_minus_b
    bootstrap_median_delta_a_minus_b

    bootstrap_q025
    bootstrap_q975

    P(delta < 0)
    P(delta > 0)
    P(delta = 0).

For the artifact descriptors:

    mean_abs_median_dx
    mean_abs_ppos_bias
    abs_spearman_dx_vs_depth_change,

a negative A-minus-B difference means that A has the smaller artifact
descriptor.

For dispersion descriptors such as MSD or MAD, the signed bootstrap difference
is descriptive only and has no better/worse interpretation.

No bootstrap winner fraction is computed.


Sensitivity analyses
--------------------
In addition to the primary abundance range, the script evaluates:

    low_boundary_trim

and:

    central_range.

These are defined from:

    --sensitivity-low-q
    --sensitivity-high-q

with defaults:

    0.05
    0.95.

An optional:

    reliable_range

can also be added using:

    --reliable-x-min
    --reliable-x-max.

All scenarios reuse the SAME primary fixed bin edges.

Thus, sensitivity analyses modify the included abundance range rather than
re-estimating the binning geometry.


Output
------
The Step-8 output directory contains:


00_run_config.json
------------------
Complete run configuration including:

    script and version;
    input path;
    conditioning coordinate;
    representation mapping;
    relative-frequency candidate labels;
    control representation labels;
    closure/depth columns;
    dt and observability filters;
    binning parameters;
    subject-bootstrap settings;
    sensitivity-range settings;
    random seed;
    metric definitions;
    automatic_winner_assigned = False.


01_analysis_row_audit.csv
-------------------------
Records:

    input row count;
    rows after dt filtering;
    rows after observability-class filtering;
    complete-case representation rows;
    rows retained within conditioning range;
    rows with valid depth pairs;
    subject count;
    conditioning and representation definitions;
    bin-edge source;
    analysis abundance range;
    number of fixed bins;
    winner policy.


02_conditioning_bin_edges.csv
-----------------------------
Fixed xstar conditioning bins:

    bin_id
    x_lo
    x_hi
    x_center.


03_pooled_representation_curves.csv
-----------------------------------
Per-representation, per-bin displacement statistics including:

    n
    n_subjects
    x_median

    median_dx
    mean_dx
    ppos

    msd
    var_dx
    mad_dx

    q025_dx
    q975_dx

    abs_median_dx
    ppos_bias
    meets_min_n.


04_pooled_representation_assessment.csv
---------------------------------------
Primary pooled representation summary containing directional-null, dispersion
and depth-dependence descriptors together with explicit representation roles.

No automatic ranking is assigned.


05_subject_representation_assessment.csv
----------------------------------------
Subject-level representation summaries using the same fixed pooled bin edges.


06_subject_bootstrap_pairwise_effects.csv
-----------------------------------------
Subject-bootstrap pairwise representation effect-size distributions.


07_pairwise_displacement_agreement.csv
--------------------------------------
Transition-level pairwise agreement between displacement representations.


08_depth_dependence_summary.csv
-------------------------------
Pooled residual sequencing-depth dependence for every displacement
representation.


09_closure_diagnostic_curves.csv
--------------------------------
Closure-component summaries across conditioning-abundance bins.


10_closure_depth_summary.csv
----------------------------
Whole-dataset association between the closure component and sequencing-depth
change.


11_representation_assessment.csv
--------------------------------
A manuscript-oriented copy of the primary pooled representation assessment.


12_sensitivity_representation_assessment.csv
--------------------------------------------
Representation assessment across the primary and alternative abundance-range
scenarios.


13_sensitivity_bootstrap_pairwise_effects.csv
---------------------------------------------
Subject-bootstrap pairwise representation contrasts across sensitivity
scenarios.


representation_assessment_report.md
-----------------------------------
Human-readable summary of:

    data audit;
    primary representation assessment;
    closure/depth diagnostic;
    sensitivity analyses;
    non-ranking interpretation policy.


Primary run
-----------
A standard Step-8 run is:

    python3 8-pseudo_longitudinal_representation_assessment.py \
        --input \
        ./results_pseudo/5-latent_transition_construction/latent_transitions.parquet \
        --outdir \
        ./results_pseudo/8-pseudo_longitudinal_representation_assessment \
        --conditioning-col xstar_latent \
        --representations \
        latent_frequency=dx_latent,observed_frequency=dx_freq_obs,log_count=dx_count_sum \
        --relative-frequency-representations \
        latent_frequency,observed_frequency \
        --closure-col closure_component \
        --dt-values 1 \
        --classes TT \
        --min-n 200 \
        --min-common-bins 10 \
        --min-n-subject 50 \
        --min-common-bins-subject 6 \
        --n-subject-bootstrap 2000 \
        --seed 123


Interpretation
--------------
Step 8 asks:

    "Once transition conditioning is fixed at xstar, which properties of the
     inferred pseudo-longitudinal displacement depend on whether abundance
     change is represented by latent frequency, observed frequency or counts?"

The answer is deliberately decomposed into:

    directional-null adequacy;
    residual depth dependence;
    closure/count behavior;
    descriptive fluctuation magnitude.

No single scalar score is constructed and no representation is automatically
declared the winner.

For the current ClonoDynamics framework, this separation is essential:
directional bias, observation-depth dependence and fluctuation amplitude are
different statistical properties and should not be collapsed into one
optimization criterion.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import polars as pl
except Exception:  # optional fallback
    pl = None


SCRIPT_VERSION = "v3-xstar-pseudo-longitudinal-representation-assessment-2026-08-24"

DEFAULT_REPRESENTATIONS = (
    "latent_frequency=dx_latent,"
    "observed_frequency=dx_freq_obs,"
    "log_count=dx_count_sum"
)
DEFAULT_RELATIVE_FREQUENCY_REPRESENTATIONS = (
    "latent_frequency,observed_frequency"
)

DIRECTIONAL_METRICS = [
    "mean_abs_mean_dx",
    "mean_abs_median_dx",
    "mean_abs_ppos_bias",
    "rms_median_dx",
    "max_abs_median_dx",
    "rms_ppos_bias",
    "max_abs_ppos_bias",
]
DEPTH_METRICS = [
    "abs_spearman_dx_vs_depth_change",
    "spearman_dx_vs_depth_change",
    "pearson_dx_vs_depth_change",
    "ols_slope_dx_vs_depth_change",
    "spearman_sqdx_vs_mean_log_depth",
]
DISPERSION_METRICS = [
    "mean_msd",
    "median_msd",
    "mean_mad_dx",
    "median_mad_dx",
    "msd_log_roughness",
]
BOOTSTRAP_EFFECT_METRICS = [
    "mean_abs_mean_dx",
    "mean_abs_median_dx",
    "mean_abs_ppos_bias",
    "abs_spearman_dx_vs_depth_change",
    "mean_msd",
    "mean_mad_dx",
]
LOWER_IS_SMALLER_ARTIFACT = {
    "mean_abs_median_dx",
    "mean_abs_ppos_bias",
    "abs_spearman_dx_vs_depth_change",
}


# =============================================================================
# Generic utilities
# =============================================================================


def parse_csv_values(text: str) -> List[str]:
    return [value.strip() for value in str(text).split(",") if value.strip()]


def parse_int_values(text: str) -> List[int]:
    return [int(value) for value in parse_csv_values(text)]


def parse_mapping(text: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for item in parse_csv_values(text):
        if "=" not in item:
            raise ValueError(
                "Representation mappings must use label=column; got {!r}".format(item)
            )
        label, column = item.split("=", 1)
        label = label.strip()
        column = column.strip()
        if not label or not column:
            raise ValueError("Invalid representation mapping: {!r}".format(item))
        if label in mapping:
            raise ValueError("Duplicate representation label: {}".format(label))
        mapping[label] = column
    if len(mapping) < 2:
        raise ValueError("At least two representations are required")
    return mapping


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if np.sum(valid) < 3:
        return np.nan
    xr = pd.Series(x[valid]).rank(method="average").to_numpy(dtype=float)
    yr = pd.Series(y[valid]).rank(method="average").to_numpy(dtype=float)
    if np.nanstd(xr) == 0 or np.nanstd(yr) == 0:
        return np.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if np.sum(valid) < 3:
        return np.nan
    xv = x[valid]
    yv = y[valid]
    if np.nanstd(xv) == 0 or np.nanstd(yv) == 0:
        return np.nan
    return float(np.corrcoef(xv, yv)[0, 1])


def safe_ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if np.sum(valid) < 3:
        return np.nan
    xv = x[valid]
    yv = y[valid]
    var = float(np.var(xv, ddof=0))
    if not np.isfinite(var) or var <= 0:
        return np.nan
    cov = float(np.mean((xv - np.mean(xv)) * (yv - np.mean(yv))))
    return cov / var


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        if np.isnan(value):
            return None
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


# =============================================================================
# Input loading
# =============================================================================


def _polars_schema_names(lf) -> List[str]:
    try:
        return list(lf.collect_schema().names())
    except Exception:
        try:
            return list(lf.schema.keys())
        except Exception:
            return []


def _collect_polars(lf):
    try:
        return lf.collect(engine="streaming")
    except Exception:
        try:
            return lf.collect(streaming=True)
        except Exception:
            return lf.collect()


def load_transition_subset(
    path: Path,
    wanted_columns: Sequence[str],
    dt_values: Sequence[int],
    classes: Sequence[str],
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    suffix = path.suffix.lower()
    audit = {
        "input_rows": -1,
        "rows_after_dt_filter": -1,
        "rows_after_class_filter": -1,
    }

    if pl is not None and suffix in {".parquet", ".pq", ".csv", ".tsv", ".txt"}:
        if suffix in {".parquet", ".pq"}:
            lf = pl.scan_parquet(str(path))
        else:
            separator = "\t" if suffix in {".tsv", ".txt"} else ","
            lf = pl.scan_csv(str(path), separator=separator, infer_schema_length=1000)

        schema_names = _polars_schema_names(lf)
        missing = [c for c in wanted_columns if c not in schema_names]
        if missing:
            raise ValueError("Missing required projected columns: {}".format(missing))

        try:
            audit["input_rows"] = int(_collect_polars(lf.select(pl.len())).item())
        except Exception:
            pass

        lf_dt = lf.filter(pl.col("dt").is_in([int(v) for v in dt_values]))
        try:
            audit["rows_after_dt_filter"] = int(
                _collect_polars(lf_dt.select(pl.len())).item()
            )
        except Exception:
            pass

        lf_class = lf_dt
        if classes:
            lf_class = lf_class.filter(pl.col("obs_class").is_in(list(classes)))
        try:
            audit["rows_after_class_filter"] = int(
                _collect_polars(lf_class.select(pl.len())).item()
            )
        except Exception:
            pass

        df = _collect_polars(lf_class.select(list(wanted_columns))).to_pandas()
        return df, audit

    if suffix in {".parquet", ".pq"}:
        df = pd.read_parquet(path, columns=list(wanted_columns))
        audit["input_rows"] = int(len(pd.read_parquet(path, columns=["dt"])))
    elif suffix == ".csv":
        df = pd.read_csv(path, usecols=list(wanted_columns))
        audit["input_rows"] = int(len(df))
    elif suffix in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t", usecols=list(wanted_columns))
        audit["input_rows"] = int(len(df))
    else:
        raise ValueError("Unsupported input format: {}".format(path))

    dt_mask = pd.to_numeric(df["dt"], errors="coerce").isin(dt_values)
    audit["rows_after_dt_filter"] = int(dt_mask.sum())
    df = df.loc[dt_mask].copy()
    if classes:
        class_mask = df["obs_class"].astype(str).isin(classes)
        audit["rows_after_class_filter"] = int(class_mask.sum())
        df = df.loc[class_mask].copy()
    else:
        audit["rows_after_class_filter"] = int(len(df))
    return df, audit


def available_columns(path: Path) -> List[str]:
    suffix = path.suffix.lower()
    if pl is not None and suffix in {".parquet", ".pq", ".csv", ".tsv", ".txt"}:
        if suffix in {".parquet", ".pq"}:
            lf = pl.scan_parquet(str(path))
        else:
            separator = "\t" if suffix in {".tsv", ".txt"} else ","
            lf = pl.scan_csv(str(path), separator=separator, infer_schema_length=1000)
        return _polars_schema_names(lf)
    if suffix == ".csv":
        return list(pd.read_csv(path, nrows=0).columns)
    if suffix in {".tsv", ".txt"}:
        return list(pd.read_csv(path, sep="\t", nrows=0).columns)
    if suffix in {".parquet", ".pq"}:
        try:
            import pyarrow.parquet as pq
            return list(pq.ParquetFile(path).schema.names)
        except Exception:
            return list(pd.read_parquet(path).columns)
    raise ValueError("Unsupported input format: {}".format(path))


# =============================================================================
# Fixed bins and common rows
# =============================================================================


def read_edges_csv(path: Path) -> np.ndarray:
    d = pd.read_csv(path)
    if {"x_lo", "x_hi"}.issubset(d.columns):
        lo = pd.to_numeric(d["x_lo"], errors="coerce").to_numpy(dtype=float)
        hi = pd.to_numeric(d["x_hi"], errors="coerce").to_numpy(dtype=float)
        edges = np.concatenate([lo[:1], hi])
    elif {"bin_left", "bin_right"}.issubset(d.columns):
        lo = pd.to_numeric(d["bin_left"], errors="coerce").to_numpy(dtype=float)
        hi = pd.to_numeric(d["bin_right"], errors="coerce").to_numpy(dtype=float)
        edges = np.concatenate([lo[:1], hi])
    elif "edge" in d.columns:
        edges = pd.to_numeric(d["edge"], errors="coerce").to_numpy(dtype=float)
    else:
        raise ValueError(
            "Edges CSV must contain x_lo/x_hi, bin_left/bin_right, or edge"
        )
    edges = np.asarray(edges, dtype=float)
    edges = np.unique(edges[np.isfinite(edges)])
    if len(edges) < 3 or not np.all(np.diff(edges) > 0):
        raise ValueError("Invalid fixed bin edges in {}".format(path))
    return edges


def derive_edges(x: np.ndarray, n_bins: int, q_low: float, q_high: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        raise ValueError("No finite conditioning values available for binning")
    x_min = float(np.quantile(x, q_low))
    x_max = float(np.quantile(x, q_high))
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min >= x_max:
        raise ValueError("Could not derive a valid conditioning range")
    return np.linspace(x_min, x_max, int(n_bins) + 1)


def assign_bin_ids(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    out = np.zeros(len(x), dtype=np.int32)
    finite = np.isfinite(x)
    if not np.any(finite):
        return out
    idx = np.searchsorted(edges, x[finite], side="right") - 1
    idx[x[finite] == edges[-1]] = len(edges) - 2
    valid = (idx >= 0) & (idx < len(edges) - 1)
    positions = np.where(finite)[0]
    out[positions[valid]] = idx[valid] + 1
    return out


def build_edges_table(edges: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "bin_id": np.arange(1, len(edges), dtype=int),
            "x_lo": edges[:-1],
            "x_hi": edges[1:],
            "x_center": 0.5 * (edges[:-1] + edges[1:]),
        }
    )


def prepare_common_rows(
    df: pd.DataFrame,
    conditioning_col: str,
    rep_map: Dict[str, str],
    closure_col: str,
    depth0_col: str,
    depth1_col: str,
    edges: np.ndarray,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    d = df.copy()
    numeric_columns = [conditioning_col] + list(rep_map.values())
    for column in [closure_col, depth0_col, depth1_col]:
        if column in d.columns:
            numeric_columns.append(column)
    for column in dict.fromkeys(numeric_columns):
        d[column] = pd.to_numeric(d[column], errors="coerce")

    common_mask = np.isfinite(d[conditioning_col].to_numpy(dtype=float))
    for column in rep_map.values():
        common_mask &= np.isfinite(d[column].to_numpy(dtype=float))
    d = d.loc[common_mask].copy()
    rows_complete = int(len(d))

    x = d[conditioning_col].to_numpy(dtype=float)
    range_mask = (x >= float(edges[0])) & (x <= float(edges[-1]))
    d = d.loc[range_mask].copy()
    d["bin_id"] = assign_bin_ids(d[conditioning_col].to_numpy(dtype=float), edges)
    d = d[d["bin_id"] > 0].copy()

    if depth0_col in d.columns and depth1_col in d.columns:
        depth0 = d[depth0_col].to_numpy(dtype=float)
        depth1 = d[depth1_col].to_numpy(dtype=float)
        valid = (
            np.isfinite(depth0)
            & np.isfinite(depth1)
            & (depth0 > 0)
            & (depth1 > 0)
        )
        delta = np.full(len(d), np.nan, dtype=float)
        mean_depth = np.full(len(d), np.nan, dtype=float)
        delta[valid] = np.log(depth1[valid]) - np.log(depth0[valid])
        mean_depth[valid] = 0.5 * (np.log(depth1[valid]) + np.log(depth0[valid]))
        d["delta_log_depth"] = delta
        d["mean_log_depth"] = mean_depth
        n_depth_valid = int(np.sum(valid))
    else:
        d["delta_log_depth"] = np.nan
        d["mean_log_depth"] = np.nan
        n_depth_valid = 0

    audit = {
        "rows_after_complete_representation_case_filter": rows_complete,
        "rows_after_conditioning_range_filter": int(len(d)),
        "rows_with_valid_depth_pair": n_depth_valid,
        "n_subjects": int(d["subject"].nunique()),
    }
    return d, audit


# =============================================================================
# Binned displacement summaries
# =============================================================================


def displacement_metrics(values: np.ndarray) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {
            "median_dx": np.nan,
            "mean_dx": np.nan,
            "ppos": np.nan,
            "msd": np.nan,
            "var_dx": np.nan,
            "mad_dx": np.nan,
            "q025_dx": np.nan,
            "q975_dx": np.nan,
        }
    median = float(np.median(x))
    return {
        "median_dx": median,
        "mean_dx": float(np.mean(x)),
        "ppos": float(np.mean(x > 0)),
        "msd": float(np.mean(x * x)),
        "var_dx": float(np.var(x, ddof=0)),
        "mad_dx": float(np.median(np.abs(x - median))),
        "q025_dx": float(np.quantile(x, 0.025)),
        "q975_dx": float(np.quantile(x, 0.975)),
    }


def pooled_curves(
    d: pd.DataFrame,
    conditioning_col: str,
    rep_map: Dict[str, str],
    edges: np.ndarray,
    min_n: int,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    edge_table = build_edges_table(edges).set_index("bin_id")
    for representation, dx_col in rep_map.items():
        for bin_id, group in d.groupby("bin_id", sort=True):
            bin_id = int(bin_id)
            values = group[dx_col].to_numpy(dtype=float)
            metrics = displacement_metrics(values)
            row: Dict[str, object] = {
                "representation": representation,
                "dx_col": dx_col,
                "conditioning_col": conditioning_col,
                "bin_id": bin_id,
                "bin_lo": float(edge_table.loc[bin_id, "x_lo"]),
                "bin_hi": float(edge_table.loc[bin_id, "x_hi"]),
                "x_center": float(edge_table.loc[bin_id, "x_center"]),
                "x_median": float(np.median(group[conditioning_col].to_numpy(dtype=float))),
                "n": int(len(group)),
                "n_subjects": int(group["subject"].nunique()),
            }
            row.update(metrics)
            row["abs_median_dx"] = abs(float(row["median_dx"]))
            row["ppos_bias"] = abs(float(row["ppos"]) - 0.5)
            row["meets_min_n"] = bool(len(group) >= int(min_n))
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["representation", "bin_id"]).reset_index(drop=True)


def curve_roughness_log(values: np.ndarray) -> float:
    y = np.asarray(values, dtype=float)
    y = y[np.isfinite(y) & (y > 0)]
    if len(y) < 3:
        return np.nan
    second_diff = np.diff(np.log(y), n=2)
    if len(second_diff) == 0:
        return np.nan
    return float(np.mean(np.abs(second_diff)))


def summarize_curves(curves: pd.DataFrame, min_common_bins: int) -> pd.DataFrame:
    if curves.empty:
        return pd.DataFrame()
    valid_sets: List[set] = []
    for _, group in curves.groupby("representation", sort=False):
        valid_sets.append(set(group.loc[group["meets_min_n"], "bin_id"].astype(int)))
    common_bins = set.intersection(*valid_sets) if valid_sets else set()
    if len(common_bins) < int(min_common_bins):
        raise ValueError(
            "Only {} common bins meet min_n; require at least {}".format(
                len(common_bins), min_common_bins
            )
        )

    rows = []
    for representation, group in curves.groupby("representation", sort=False):
        g = group[group["bin_id"].isin(sorted(common_bins))].sort_values("bin_id")
        mean_dx = g["mean_dx"].to_numpy(dtype=float)
        median_dx = g["median_dx"].to_numpy(dtype=float)
        ppos_bias = g["ppos_bias"].to_numpy(dtype=float)
        msd = g["msd"].to_numpy(dtype=float)
        mad = g["mad_dx"].to_numpy(dtype=float)
        rows.append(
            {
                "representation": representation,
                "dx_col": str(g["dx_col"].iloc[0]),
                "n_common_bins": int(len(g)),
                "common_bin_ids": ",".join(str(int(v)) for v in g["bin_id"]),
                "mean_abs_mean_dx": float(np.mean(np.abs(mean_dx))),
                "mean_abs_median_dx": float(np.mean(np.abs(median_dx))),
                "rms_median_dx": float(np.sqrt(np.mean(median_dx * median_dx))),
                "max_abs_median_dx": float(np.max(np.abs(median_dx))),
                "median_dx_range": float(np.max(median_dx) - np.min(median_dx)),
                "mean_abs_ppos_bias": float(np.mean(ppos_bias)),
                "rms_ppos_bias": float(np.sqrt(np.mean(ppos_bias * ppos_bias))),
                "max_abs_ppos_bias": float(np.max(ppos_bias)),
                "mean_msd": float(np.mean(msd)),
                "median_msd": float(np.median(msd)),
                "mean_mad_dx": float(np.mean(mad)),
                "median_mad_dx": float(np.median(mad)),
                "msd_log_roughness": curve_roughness_log(msd),
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# Depth dependence, representation roles, and closure
# =============================================================================


def depth_dependence_for_values(
    dx: np.ndarray,
    delta_log_depth: np.ndarray,
    mean_log_depth: np.ndarray,
) -> Dict[str, float]:
    dx = np.asarray(dx, dtype=float)
    delta = np.asarray(delta_log_depth, dtype=float)
    mean_depth = np.asarray(mean_log_depth, dtype=float)
    valid_delta = np.isfinite(dx) & np.isfinite(delta)
    valid_mean = np.isfinite(dx) & np.isfinite(mean_depth)
    spearman = safe_spearman(dx[valid_delta], delta[valid_delta])
    pearson = safe_pearson(dx[valid_delta], delta[valid_delta])
    slope = safe_ols_slope(delta[valid_delta], dx[valid_delta])
    sqdx = dx * dx
    magnitude_spearman = safe_spearman(sqdx[valid_mean], mean_depth[valid_mean])
    return {
        "n_depth_valid": int(np.sum(valid_delta)),
        "pearson_dx_vs_depth_change": pearson,
        "spearman_dx_vs_depth_change": spearman,
        "abs_spearman_dx_vs_depth_change": (
            abs(spearman) if np.isfinite(spearman) else np.nan
        ),
        "ols_slope_dx_vs_depth_change": slope,
        "spearman_sqdx_vs_mean_log_depth": magnitude_spearman,
    }


def pooled_depth_summary(d: pd.DataFrame, rep_map: Dict[str, str]) -> pd.DataFrame:
    rows = []
    delta = d["delta_log_depth"].to_numpy(dtype=float)
    mean_depth = d["mean_log_depth"].to_numpy(dtype=float)
    for representation, dx_col in rep_map.items():
        metrics = depth_dependence_for_values(
            d[dx_col].to_numpy(dtype=float), delta, mean_depth
        )
        row = {"representation": representation, "dx_col": dx_col}
        row.update(metrics)
        rows.append(row)
    return pd.DataFrame(rows)


def build_assessment_table(
    curve_summary: pd.DataFrame,
    depth_summary: pd.DataFrame,
    relative_frequency_representations: Sequence[str],
) -> pd.DataFrame:
    out = curve_summary.merge(depth_summary, on=["representation", "dx_col"], how="left")
    relative_set = set(str(v) for v in relative_frequency_representations)
    out["representation_role"] = np.where(
        out["representation"].astype(str).isin(relative_set),
        "relative_frequency_candidate",
        "non_compositional_control",
    )
    out["directional_null_metrics_role"] = (
        "descriptive_adequacy_not_global_ranking"
    )
    out["depth_dependence_role"] = (
        "residual_observation_dependence_descriptor"
    )
    out["dispersion_metrics_role"] = (
        "descriptive_not_selection_criterion"
    )
    out["automatic_winner_assigned"] = False
    preferred_order = list(relative_frequency_representations) + [
        r for r in out["representation"].astype(str).tolist()
        if r not in relative_set
    ]
    order_map = {name: i for i, name in enumerate(preferred_order)}
    out["__order"] = out["representation"].map(order_map).fillna(999)
    out = out.sort_values(["__order", "representation"]).drop(columns="__order")
    return out.reset_index(drop=True)


def closure_curves(
    d: pd.DataFrame,
    conditioning_col: str,
    closure_col: str,
    edges: np.ndarray,
    min_n: int,
) -> pd.DataFrame:
    if closure_col not in d.columns:
        return pd.DataFrame()
    edge_table = build_edges_table(edges).set_index("bin_id")
    rows = []
    for bin_id, group in d.groupby("bin_id", sort=True):
        raw = pd.to_numeric(group[closure_col], errors="coerce")
        valid = np.isfinite(raw.to_numpy(dtype=float))
        values = raw.to_numpy(dtype=float)[valid]
        if len(values) == 0:
            continue
        median = float(np.median(values))
        bin_id = int(bin_id)
        rows.append(
            {
                "bin_id": bin_id,
                "bin_lo": float(edge_table.loc[bin_id, "x_lo"]),
                "bin_hi": float(edge_table.loc[bin_id, "x_hi"]),
                "x_center": float(edge_table.loc[bin_id, "x_center"]),
                "x_median": float(np.median(group[conditioning_col].to_numpy(dtype=float))),
                "n": int(len(values)),
                "n_subjects": int(group.loc[valid, "subject"].nunique()),
                "median_closure": median,
                "mean_closure": float(np.mean(values)),
                "mad_closure": float(np.median(np.abs(values - median))),
                "msd_closure": float(np.mean(values * values)),
                "ppos_closure": float(np.mean(values > 0)),
                "meets_min_n": bool(len(values) >= int(min_n)),
            }
        )
    return pd.DataFrame(rows)


def closure_depth_summary(d: pd.DataFrame, closure_col: str) -> pd.DataFrame:
    if closure_col not in d.columns:
        return pd.DataFrame()
    closure = pd.to_numeric(d[closure_col], errors="coerce").to_numpy(dtype=float)
    delta = d["delta_log_depth"].to_numpy(dtype=float)
    valid = np.isfinite(closure) & np.isfinite(delta)
    return pd.DataFrame(
        [
            {
                "closure_col": closure_col,
                "n_depth_valid": int(np.sum(valid)),
                "pearson_closure_vs_depth_change": safe_pearson(
                    closure[valid], delta[valid]
                ),
                "spearman_closure_vs_depth_change": safe_spearman(
                    closure[valid], delta[valid]
                ),
                "ols_slope_closure_vs_depth_change": safe_ols_slope(
                    delta[valid], closure[valid]
                ),
            }
        ]
    )


def pairwise_displacement_agreement(
    d: pd.DataFrame,
    rep_map: Dict[str, str],
) -> pd.DataFrame:
    rows = []
    representations = list(rep_map.keys())
    for a, b in combinations(representations, 2):
        xa = pd.to_numeric(d[rep_map[a]], errors="coerce").to_numpy(dtype=float)
        xb = pd.to_numeric(d[rep_map[b]], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(xa) & np.isfinite(xb)
        if np.sum(valid) == 0:
            continue
        delta = xa[valid] - xb[valid]
        rows.append(
            {
                "candidate_a": a,
                "candidate_b": b,
                "pair": "{}_vs_{}".format(a, b),
                "n": int(np.sum(valid)),
                "pearson": safe_pearson(xa[valid], xb[valid]),
                "spearman": safe_spearman(xa[valid], xb[valid]),
                "median_difference_a_minus_b": float(np.median(delta)),
                "median_abs_difference": float(np.median(np.abs(delta))),
                "rms_difference": float(np.sqrt(np.mean(delta * delta))),
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# Subject summaries and pairwise bootstrap effect sizes
# =============================================================================


def subject_representation_assessment(
    d: pd.DataFrame,
    conditioning_col: str,
    rep_map: Dict[str, str],
    relative_frequency_representations: Sequence[str],
    edges: np.ndarray,
    min_n_subject: int,
    min_common_bins_subject: int,
) -> pd.DataFrame:
    rows = []
    for subject, sd in d.groupby("subject", sort=True):
        curves = pooled_curves(
            sd,
            conditioning_col=conditioning_col,
            rep_map=rep_map,
            edges=edges,
            min_n=min_n_subject,
        )
        try:
            curve_summary = summarize_curves(curves, min_common_bins_subject)
        except ValueError:
            continue
        depth_summary = pooled_depth_summary(sd, rep_map)
        assessment = build_assessment_table(
            curve_summary, depth_summary, relative_frequency_representations
        )
        assessment.insert(0, "subject", subject)
        rows.append(assessment)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def bootstrap_pairwise_effects(
    subject_assessment: pd.DataFrame,
    representations: Sequence[str],
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    if subject_assessment.empty:
        return pd.DataFrame()

    rng = np.random.default_rng(int(seed))
    rows: List[Dict[str, object]] = []
    for a, b in combinations(representations, 2):
        pair = "{}_vs_{}".format(a, b)
        for metric in BOOTSTRAP_EFFECT_METRICS:
            if metric not in subject_assessment.columns:
                continue
            pivot = subject_assessment.pivot(
                index="subject", columns="representation", values=metric
            )
            if a not in pivot.columns or b not in pivot.columns:
                continue
            pair_data = pivot[[a, b]].dropna()
            if pair_data.empty:
                continue
            delta_subject = (
                pair_data[a].to_numpy(dtype=float)
                - pair_data[b].to_numpy(dtype=float)
            )
            n_subjects = len(delta_subject)
            boot = np.empty(int(n_boot), dtype=float)
            for i in range(int(n_boot)):
                idx = rng.integers(0, n_subjects, size=n_subjects)
                boot[i] = float(np.mean(delta_subject[idx]))

            metric_role = (
                "artifact_descriptor_lower_is_smaller"
                if metric in LOWER_IS_SMALLER_ARTIFACT
                else "descriptive_effect_size_no_preferred_direction"
            )
            interpretation = (
                "negative delta means candidate_a has smaller artifact descriptor"
                if metric in LOWER_IS_SMALLER_ARTIFACT
                else "signed descriptive difference; no better/worse interpretation"
            )
            rows.append(
                {
                    "candidate_a": a,
                    "candidate_b": b,
                    "pair": pair,
                    "metric": metric,
                    "metric_role": metric_role,
                    "n_complete_subjects": int(n_subjects),
                    "subject_mean_delta_a_minus_b": float(np.mean(delta_subject)),
                    "subject_median_delta_a_minus_b": float(np.median(delta_subject)),
                    "bootstrap_mean_delta_a_minus_b": float(np.mean(boot)),
                    "bootstrap_median_delta_a_minus_b": float(np.median(boot)),
                    "bootstrap_q025": float(np.quantile(boot, 0.025)),
                    "bootstrap_q975": float(np.quantile(boot, 0.975)),
                    "p_delta_lt_zero": float(np.mean(boot < 0)),
                    "p_delta_gt_zero": float(np.mean(boot > 0)),
                    "p_delta_eq_zero": float(np.mean(boot == 0)),
                    "n_boot": int(n_boot),
                    "interpretation": interpretation,
                }
            )
    return pd.DataFrame(rows)


# =============================================================================
# Scenario analysis
# =============================================================================


def analyze_scenario(
    scenario: str,
    d: pd.DataFrame,
    conditioning_col: str,
    rep_map: Dict[str, str],
    relative_frequency_representations: Sequence[str],
    edges: np.ndarray,
    x_min: float,
    x_max: float,
    min_n: int,
    min_common_bins: int,
    min_n_subject: int,
    min_common_bins_subject: int,
    n_boot: int,
    seed: int,
) -> Dict[str, pd.DataFrame]:
    mask = (
        (d[conditioning_col].to_numpy(dtype=float) >= float(x_min))
        & (d[conditioning_col].to_numpy(dtype=float) <= float(x_max))
    )
    sd = d.loc[mask].copy()
    if sd.empty:
        raise ValueError("Scenario {} contains no rows".format(scenario))

    curves = pooled_curves(sd, conditioning_col, rep_map, edges, min_n)
    curve_summary = summarize_curves(curves, min_common_bins)
    depth_summary = pooled_depth_summary(sd, rep_map)
    assessment = build_assessment_table(
        curve_summary, depth_summary, relative_frequency_representations
    )
    subject_assessment = subject_representation_assessment(
        sd,
        conditioning_col,
        rep_map,
        relative_frequency_representations,
        edges,
        min_n_subject,
        min_common_bins_subject,
    )
    pairwise = bootstrap_pairwise_effects(
        subject_assessment,
        list(rep_map.keys()),
        n_boot=n_boot,
        seed=seed,
    )

    metadata = {
        "scenario": scenario,
        "scenario_x_min": float(x_min),
        "scenario_x_max": float(x_max),
        "scenario_n_rows": int(len(sd)),
        "scenario_n_subjects": int(sd["subject"].nunique()),
    }
    for table in [curves, curve_summary, depth_summary, assessment, subject_assessment, pairwise]:
        if table.empty:
            continue
        for pos, (key, value) in enumerate(metadata.items()):
            table.insert(pos, key, value)

    return {
        "data": sd,
        "curves": curves,
        "curve_summary": curve_summary,
        "depth_summary": depth_summary,
        "assessment": assessment,
        "subject_assessment": subject_assessment,
        "pairwise": pairwise,
    }


def drop_scenario_columns(table: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "scenario",
        "scenario_x_min",
        "scenario_x_max",
        "scenario_n_rows",
        "scenario_n_subjects",
    ]
    return table.drop(columns=[c for c in cols if c in table.columns])


# =============================================================================
# Report
# =============================================================================


def _format_md_value(value) -> str:
    if isinstance(value, (float, np.floating)):
        return "{:.6g}".format(value) if np.isfinite(value) else "NA"
    return str(value)


def write_report(
    path: Path,
    audit: pd.DataFrame,
    primary_assessment: pd.DataFrame,
    closure_depth: pd.DataFrame,
    sensitivity_assessment: pd.DataFrame,
) -> None:
    lines = [
        "# Step 8 — pseudo-longitudinal displacement-representation assessment",
        "",
        "## Interpretation framework",
        "",
        "This analysis does not assign a global winner. The uncertainty-weighted transition coordinate xstar is used as the conditioning coordinate for finite-time transition magnitude. Directional-null descriptors (mean |mean Δx|, mean |median Δx| and mean |P(Δx>0) − 0.5|) assess residual pseudo-longitudinal directional structure. Because displacement representations have different scales, absolute mean/median displacement amplitudes are descriptive within representation and are not interpreted as cross-scale optimization criteria. Residual observation dependence is quantified separately using |Spearman(Δx, Δlog depth)|. MSD and MAD are descriptive fluctuation-amplitude measures and are not treated as optimization criteria.",
        "",
        "## Data audit",
        "",
    ]
    for _, row in audit.iterrows():
        lines.append("- **{}:** {}".format(row["item"], row["value"]))

    lines.extend(["", "## Primary representation assessment", ""])
    if not primary_assessment.empty:
        cols = [
            "representation",
            "representation_role",
            "mean_abs_mean_dx",
            "mean_abs_median_dx",
            "mean_abs_ppos_bias",
            "abs_spearman_dx_vs_depth_change",
            "mean_msd",
            "mean_mad_dx",
        ]
        cols = [c for c in cols if c in primary_assessment.columns]
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
        for _, row in primary_assessment.iterrows():
            lines.append(
                "| " + " | ".join(_format_md_value(row[c]) for c in cols) + " |"
            )

    lines.extend(["", "## Closure diagnostic", ""])
    if closure_depth.empty:
        lines.append("Closure/depth diagnostic unavailable.")
    else:
        row = closure_depth.iloc[0]
        lines.append(
            "- Spearman(closure component, Δlog depth) = {:.6g}".format(
                row["spearman_closure_vs_depth_change"]
            )
        )
        lines.append(
            "- OLS slope(closure component ~ Δlog depth) = {:.6g}".format(
                row["ols_slope_closure_vs_depth_change"]
            )
        )

    lines.extend(["", "## Sensitivity ranges", ""])
    if not sensitivity_assessment.empty:
        compact_cols = [
            "scenario",
            "representation",
            "scenario_x_min",
            "scenario_x_max",
            "mean_abs_mean_dx",
            "mean_abs_median_dx",
            "mean_abs_ppos_bias",
            "abs_spearman_dx_vs_depth_change",
        ]
        compact_cols = [c for c in compact_cols if c in sensitivity_assessment.columns]
        lines.append("| " + " | ".join(compact_cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(compact_cols)) + " |")
        for _, row in sensitivity_assessment.iterrows():
            lines.append(
                "| "
                + " | ".join(_format_md_value(row[c]) for c in compact_cols)
                + " |"
            )

    lines.extend(
        [
            "",
            "## Decision note",
            "",
            "No automatic representation selection is encoded in this script. The scientific interpretation should distinguish residual directional-null structure from residual depth dependence and from descriptive fluctuation magnitude. Cross-representation comparisons of absolute displacement amplitude are scale-dependent and must not be used as a winner criterion. Log-count is retained as a non-compositional control, not as a relative-frequency candidate.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# =============================================================================
# CLI and main
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Step 8: assess pseudo-longitudinal displacement representations using "
            "xstar transition-centered conditioning without assigning a global winner."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--conditioning-col", default="xstar_latent")
    parser.add_argument("--representations", default=DEFAULT_REPRESENTATIONS)
    parser.add_argument(
        "--relative-frequency-representations",
        default=DEFAULT_RELATIVE_FREQUENCY_REPRESENTATIONS,
        help=(
            "Comma-separated labels treated as relative-frequency candidates. "
            "Other representations are retained as controls; no winner is assigned."
        ),
    )
    parser.add_argument("--closure-col", default="closure_component")
    parser.add_argument("--depth0-col", default="depth_sum_t0")
    parser.add_argument("--depth1-col", default="depth_sum_t1")
    parser.add_argument("--dt-values", default="1")
    parser.add_argument("--classes", default="TT")
    parser.add_argument("--edges-csv", type=Path, default=None)
    parser.add_argument("--n-bins", type=int, default=30)
    parser.add_argument("--range-q-low", type=float, default=0.01)
    parser.add_argument("--range-q-high", type=float, default=0.99)
    parser.add_argument("--min-n", type=int, default=200)
    parser.add_argument("--min-common-bins", type=int, default=10)
    parser.add_argument("--min-n-subject", type=int, default=50)
    parser.add_argument("--min-common-bins-subject", type=int, default=6)
    parser.add_argument("--n-subject-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--sensitivity-low-q", type=float, default=0.05)
    parser.add_argument("--sensitivity-high-q", type=float, default=0.95)
    parser.add_argument("--reliable-x-min", type=float, default=None)
    parser.add_argument("--reliable-x-max", type=float, default=None)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.input.exists():
        raise FileNotFoundError("Input transition table not found: {}".format(args.input))
    if args.edges_csv is not None and not args.edges_csv.exists():
        raise FileNotFoundError("Edges CSV not found: {}".format(args.edges_csv))
    if args.n_bins < 3:
        raise ValueError("--n-bins must be >= 3")
    if not (0 <= args.range_q_low < args.range_q_high <= 1):
        raise ValueError("Require 0 <= range-q-low < range-q-high <= 1")
    if not (0 <= args.sensitivity_low_q < args.sensitivity_high_q <= 1):
        raise ValueError(
            "Require 0 <= sensitivity-low-q < sensitivity-high-q <= 1"
        )
    if args.min_n < 1 or args.min_n_subject < 1:
        raise ValueError("Minimum bin counts must be positive")
    if args.n_subject_bootstrap < 1:
        raise ValueError("--n-subject-bootstrap must be positive")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)

    rep_map = parse_mapping(args.representations)
    relative_reps = parse_csv_values(args.relative_frequency_representations)
    missing_relative = [v for v in relative_reps if v not in rep_map]
    if missing_relative:
        raise ValueError(
            "Relative-frequency representations not present in --representations: {}".format(
                missing_relative
            )
        )
    if len(relative_reps) < 2:
        raise ValueError(
            "At least two --relative-frequency-representations are required"
        )

    dt_values = parse_int_values(args.dt_values)
    classes = parse_csv_values(args.classes)
    args.outdir.mkdir(parents=True, exist_ok=True)

    available = available_columns(args.input)
    wanted = [
        "subject",
        "aaSeqCDR3",
        "dt",
        "obs_class",
        args.conditioning_col,
    ] + list(rep_map.values())
    required_missing = [column for column in wanted if column not in available]
    if required_missing:
        raise ValueError("Missing required columns: {}".format(required_missing))
    for optional in [args.closure_col, args.depth0_col, args.depth1_col]:
        if optional in available:
            wanted.append(optional)
    wanted = list(dict.fromkeys(wanted))

    df, load_audit = load_transition_subset(
        args.input,
        wanted_columns=wanted,
        dt_values=dt_values,
        classes=classes,
    )

    complete_mask = np.isfinite(
        pd.to_numeric(df[args.conditioning_col], errors="coerce")
    )
    for column in rep_map.values():
        complete_mask &= np.isfinite(pd.to_numeric(df[column], errors="coerce"))
    complete_for_edges = df.loc[complete_mask].copy()
    if complete_for_edges.empty:
        raise ValueError("No common complete rows across representations")

    if args.edges_csv is not None:
        edges = read_edges_csv(args.edges_csv)
        edge_source = "fixed_edges_csv"
    else:
        edges = derive_edges(
            pd.to_numeric(
                complete_for_edges[args.conditioning_col], errors="coerce"
            ).to_numpy(dtype=float),
            args.n_bins,
            args.range_q_low,
            args.range_q_high,
        )
        edge_source = "conditioning_quantile_range"

    d, prep_audit = prepare_common_rows(
        df,
        conditioning_col=args.conditioning_col,
        rep_map=rep_map,
        closure_col=args.closure_col,
        depth0_col=args.depth0_col,
        depth1_col=args.depth1_col,
        edges=edges,
    )
    if d.empty:
        raise ValueError("No rows remain after common-row and conditioning-range filters")

    primary = analyze_scenario(
        scenario="primary",
        d=d,
        conditioning_col=args.conditioning_col,
        rep_map=rep_map,
        relative_frequency_representations=relative_reps,
        edges=edges,
        x_min=float(edges[0]),
        x_max=float(edges[-1]),
        min_n=args.min_n,
        min_common_bins=args.min_common_bins,
        min_n_subject=args.min_n_subject,
        min_common_bins_subject=args.min_common_bins_subject,
        n_boot=args.n_subject_bootstrap,
        seed=args.seed,
    )

    closure_curve = closure_curves(
        primary["data"],
        args.conditioning_col,
        args.closure_col,
        edges,
        args.min_n,
    )
    closure_depth = closure_depth_summary(primary["data"], args.closure_col)
    agreement = pairwise_displacement_agreement(primary["data"], rep_map)

    x_all = d[args.conditioning_col].to_numpy(dtype=float)
    x_finite = x_all[np.isfinite(x_all)]
    q_low = float(np.quantile(x_finite, args.sensitivity_low_q))
    q_high = float(np.quantile(x_finite, args.sensitivity_high_q))
    scenarios = [
        ("low_boundary_trim", max(float(edges[0]), q_low), float(edges[-1])),
        (
            "central_range",
            max(float(edges[0]), q_low),
            min(float(edges[-1]), q_high),
        ),
    ]
    if args.reliable_x_min is not None or args.reliable_x_max is not None:
        reliable_min = (
            float(edges[0])
            if args.reliable_x_min is None
            else max(float(edges[0]), float(args.reliable_x_min))
        )
        reliable_max = (
            float(edges[-1])
            if args.reliable_x_max is None
            else min(float(edges[-1]), float(args.reliable_x_max))
        )
        if reliable_min >= reliable_max:
            raise ValueError("Reliable-range limits do not overlap the primary range")
        scenarios.append(("reliable_range", reliable_min, reliable_max))

    scenario_results = [primary]
    for index, (name, x_min, x_max) in enumerate(scenarios, start=1):
        scenario_results.append(
            analyze_scenario(
                scenario=name,
                d=d,
                conditioning_col=args.conditioning_col,
                rep_map=rep_map,
                relative_frequency_representations=relative_reps,
                edges=edges,
                x_min=x_min,
                x_max=x_max,
                min_n=args.min_n,
                min_common_bins=args.min_common_bins,
                min_n_subject=args.min_n_subject,
                min_common_bins_subject=args.min_common_bins_subject,
                n_boot=args.n_subject_bootstrap,
                seed=args.seed + index,
            )
        )

    audit_rows = []
    for key, value in load_audit.items():
        audit_rows.append({"item": key, "value": value})
    for key, value in prep_audit.items():
        audit_rows.append({"item": key, "value": value})
    audit_rows.extend(
        [
            {"item": "dt_values", "value": ",".join(str(v) for v in dt_values)},
            {"item": "obs_classes", "value": ",".join(classes)},
            {"item": "conditioning_col", "value": args.conditioning_col},
            {"item": "conditioning_role", "value": "transition_centered_fluctuation_benchmark"},
            {"item": "representation_map", "value": json.dumps(rep_map, sort_keys=True)},
            {
                "item": "relative_frequency_representations",
                "value": ",".join(relative_reps),
            },
            {
                "item": "closure_col_available",
                "value": bool(args.closure_col in d.columns),
            },
            {
                "item": "depth_columns_available",
                "value": bool(
                    args.depth0_col in d.columns and args.depth1_col in d.columns
                ),
            },
            {"item": "bin_edge_source", "value": edge_source},
            {"item": "x_min", "value": float(edges[0])},
            {"item": "x_max", "value": float(edges[-1])},
            {"item": "n_fixed_bins", "value": int(len(edges) - 1)},
            {"item": "n_subjects", "value": int(d["subject"].nunique())},
            {"item": "automatic_winner_assigned", "value": False},
        ]
    )
    audit_df = pd.DataFrame(audit_rows)

    primary_assessment = drop_scenario_columns(primary["assessment"])
    primary_subject = drop_scenario_columns(primary["subject_assessment"])
    primary_pairwise = drop_scenario_columns(primary["pairwise"])
    primary_depth = drop_scenario_columns(primary["depth_summary"])
    primary_curves = drop_scenario_columns(primary["curves"])

    sensitivity_assessment = pd.concat(
        [r["assessment"] for r in scenario_results], ignore_index=True
    )
    sensitivity_pairwise = pd.concat(
        [r["pairwise"] for r in scenario_results if not r["pairwise"].empty],
        ignore_index=True,
    ) if any(not r["pairwise"].empty for r in scenario_results) else pd.DataFrame()

    write_csv(audit_df, args.outdir / "01_analysis_row_audit.csv")
    write_csv(build_edges_table(edges), args.outdir / "02_conditioning_bin_edges.csv")
    write_csv(primary_curves, args.outdir / "03_pooled_representation_curves.csv")
    write_csv(primary_assessment, args.outdir / "04_pooled_representation_assessment.csv")
    write_csv(primary_subject, args.outdir / "05_subject_representation_assessment.csv")
    write_csv(primary_pairwise, args.outdir / "06_subject_bootstrap_pairwise_effects.csv")
    write_csv(agreement, args.outdir / "07_pairwise_displacement_agreement.csv")
    write_csv(primary_depth, args.outdir / "08_depth_dependence_summary.csv")
    write_csv(closure_curve, args.outdir / "09_closure_diagnostic_curves.csv")
    write_csv(closure_depth, args.outdir / "10_closure_depth_summary.csv")
    write_csv(primary_assessment, args.outdir / "11_representation_assessment.csv")
    write_csv(
        sensitivity_assessment,
        args.outdir / "12_sensitivity_representation_assessment.csv",
    )
    write_csv(
        sensitivity_pairwise,
        args.outdir / "13_sensitivity_bootstrap_pairwise_effects.csv",
    )

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "input": str(args.input),
        "outdir": str(args.outdir),
        "conditioning_col": args.conditioning_col,
        "conditioning_role": "transition_centered_fluctuation_benchmark",
        "representations": rep_map,
        "relative_frequency_representations": relative_reps,
        "control_representations": [r for r in rep_map if r not in set(relative_reps)],
        "closure_col": args.closure_col,
        "depth0_col": args.depth0_col,
        "depth1_col": args.depth1_col,
        "dt_values": dt_values,
        "classes": classes,
        "edges_csv": str(args.edges_csv) if args.edges_csv is not None else "",
        "edge_source": edge_source,
        "n_bins": int(len(edges) - 1),
        "range_q_low": args.range_q_low,
        "range_q_high": args.range_q_high,
        "min_n": args.min_n,
        "min_common_bins": args.min_common_bins,
        "min_n_subject": args.min_n_subject,
        "min_common_bins_subject": args.min_common_bins_subject,
        "n_subject_bootstrap": args.n_subject_bootstrap,
        "seed": args.seed,
        "sensitivity_low_q": args.sensitivity_low_q,
        "sensitivity_high_q": args.sensitivity_high_q,
        "reliable_x_min": args.reliable_x_min,
        "reliable_x_max": args.reliable_x_max,
        "automatic_winner_assigned": False,
        "directional_null_metrics": DIRECTIONAL_METRICS,
        "depth_dependence_metrics": DEPTH_METRICS,
        "dispersion_metrics": DISPERSION_METRICS,
        "bootstrap_effect_metrics": BOOTSTRAP_EFFECT_METRICS,
    }
    (args.outdir / "00_run_config.json").write_text(
        json.dumps(json_safe(config), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    write_report(
        args.outdir / "representation_assessment_report.md",
        audit_df,
        primary_assessment,
        closure_depth,
        sensitivity_assessment,
    )

    print("[INFO] Representation assessment complete")
    print("[INFO] Conditioning coordinate: {}".format(args.conditioning_col))
    print("[INFO] No automatic winner was assigned")
    print("[INFO] Output directory: {}".format(args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
