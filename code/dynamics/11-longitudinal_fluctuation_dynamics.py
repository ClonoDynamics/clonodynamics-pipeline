#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
11-longitudinal_fluctuation_dynamics.py
========================================

Abundance-resolved finite-time longitudinal fluctuation analysis for
ClonoDynamics.

Overview
--------
This script implements Step 11 of the ClonoDynamics pipeline.

It quantifies the magnitude and abundance dependence of genuine longitudinal
clonotype fluctuations across multiple temporal lags after the displacement
representation and conditioning geometry have been fixed by the preceding
technical-null benchmarks.

The manuscript-primary configuration is:

    displacement
        latent relative-frequency displacement:

            dx_latent

    conditioning coordinate
        uncertainty-weighted transition-centred latent abundance:

            xstar_latent

    transition class
        TT

    abundance coordinate
        one shared ABSOLUTE xstar grid used across all selected dt values.

For each temporal lag dt and shared abundance bin, the primary fluctuation
quantity is:

    Var(Delta x | xstar, dt).

The script additionally reports:

    MSD(xstar, dt)
        = E[(Delta x)^2 | xstar, dt]

    MAD(Delta x | xstar, dt)

and directional diagnostics:

    E[Delta x | xstar, dt]
    median(Delta x | xstar, dt)
    P(Delta x > 0 | xstar, dt).

Because genuine longitudinal forward drift is not assumed to be zero, the
conditional variance and the mean-squared displacement are kept separate.

For a conditional displacement distribution:

    MSD
        = Var(Delta x) + [E(Delta x)]^2

up to the distinction between finite-sample estimators.

Thus, MSD contains both fluctuation around the conditional mean and the squared
directional component, whereas Var(Delta x) is the primary dispersion measure
around the conditional mean.

This Step estimates abundance-resolved finite-time fluctuation statistics.
It does NOT fit a temporal scaling law or a predefined stochastic process.


Position in the ClonoDynamics pipeline
--------------------------------------
Relevant public workflow:

    Step 1
    1-repertoire_characterization.py

    Step 2
    2-noise_aware_latent_inference.py

    Step 3
    3-latent_posterior_quality_characterization.py

    Step 4
    4-latent_trajectory_construction.py

    Step 5
    5-latent_transition_construction.py
                |
                +---------------------------------------------+
                |                                             |
                v                                             v

    Steps 7-8                                    Steps 9-10
    pseudo-longitudinal technical-null           genuine forward drift
    conditioning/representation assessment       and longitudinal-vs-pseudo test
                |                                             |
                +--------------------+------------------------+
                                     |
                                     v

    Step 11
    11-longitudinal_fluctuation_dynamics.py

        latent dx
        +
        xstar conditioning
        +
        selected temporal lags
                |
                v
        shared absolute abundance bins
                |
                v
        Var(dx), MSD, MAD and directional diagnostics
                |
                v
        transition-weighted primary curves
        + equal-subject and LOSO sensitivities
        + temporal-comparability diagnostics


Step 11 reads the generic Step-5 transition table directly.

Steps 7-8 justify the primary latent/xstar architecture, whereas Steps 9-10
address forward directional structure. They are not direct file inputs to this
script.


Input
-----
The required input is a longitudinal transition table following the Step-5
schema:

    latent_transitions.parquet

supplied with:

    --transitions.

Supported formats are:

    .parquet
    .csv
    .tsv
    .txt.


Required core identifiers
-------------------------
The input must contain:

    subject
        Subject identifier.

    aaSeqCDR3
        Amino-acid CDR3 clonotype identifier.

    t0
        Initial sampling time.

    t1
        Final sampling time.

The temporal lag is read from:

    dt

when present, otherwise reconstructed as:

    dt = t1 - t0.


Effective temporal lag
----------------------
If available, the script uses the first matching field among:

    effective_dt
    dt_effective
    delta_t_effective

as a descriptive effective temporal interval.

If none is present:

    effective_dt = dt.

`effective_dt` is retained for temporal-comparability diagnostics. The primary
curves remain stratified by the nominal integer `dt`.


Observability class
-------------------
The transition class is read from:

    obs_class

when available.

If `obs_class` is absent but:

    obs0
    obs1

are present, the four classes are reconstructed as:

    TT
        observable at both endpoints;

    TF
        observable only at t0;

    FT
        observable only at t1;

    FF
        observable at neither endpoint.

The primary configuration uses:

    --mode TT.


Supported transition modes
--------------------------
The script accepts:

    TT
    TF
    FT
    FF
    ALL
    NON_TT
    NON_FF.

Definitions:

    TT
        retain TT only;

    TF / FT / FF
        retain the selected class only;

    ALL
        retain all classes;

    NON_TT
        retain every class except TT;

    NON_FF
        retain TT + TF + FT while excluding FF.

`NON_FF` is intended as a detectability-inclusive sensitivity analysis.


Representation
--------------
One displacement representation is analyzed per run.

Available choices are:

    latent
    observed
    count.

The manuscript-primary representation is:

    latent.

Non-primary representations require:

    --allow-nonprimary.


Latent representation
---------------------
The latent representation resolves:

    x0
        from x0_latent

    x1
        from x1_latent

    dx
        from dx_latent when available,
        otherwise x1_latent - x0_latent

    midpoint
        from xmid_latent when available,
        otherwise (x0 + x1)/2

    xstar
        from xstar_latent.


Observed representation
-----------------------
The observed representation searches compatible observed-frequency fields for:

    x0
    x1
    dx
    midpoint
    xstar.

No latent-xstar fallback is applied automatically.

Therefore, if:

    --representation observed
    --conditioning xstar

is requested, an observed-representation-specific xstar field must exist.


Count representation
--------------------
The count representation analogously searches count/log-count transition fields
for:

    x0
    x1
    dx
    midpoint
    xstar.

Again, xstar must belong to the selected representation; latent xstar is not
silently substituted.


Conditioning coordinate
-----------------------
The supported conditioning coordinates are:

    xstar
    x0
    midpoint.

The manuscript-primary choice is:

    xstar.

Non-primary conditioning requires:

    --allow-nonprimary.


Primary xstar conditioning
--------------------------
For the latent representation, Step 5 defines an uncertainty-weighted
transition-centred coordinate of the form:

                     w0*x0 + w1*x1
    xstar = -----------------------------,
                         w0 + w1

where the endpoint weights are based on posterior uncertainty.

Step 11 does not recompute xstar.

It uses the Step-5 transition-level `xstar_latent` value as the conditioning
coordinate.


Core transition preparation
---------------------------
For each retained transition, Step 11 standardizes:

    subject
    aaSeqCDR3
    t0
    t1
    dt
    effective_dt
    obs_class

    x0
    x1
    xmid
    x_condition
    dx.

Rows are retained only when:

    subject is defined;
    t0 and t1 are defined;
    dt is defined;
    x_condition is finite;
    dx is finite.

The requested dt and observability filters are then applied.


Posterior and detectability diagnostics
---------------------------------------
When present in the Step-5 transition table, the following quantities are
propagated:

    p_detect_state_t0
    p_detect_state_t1

    p_dropout_state_t0
    p_dropout_state_t1

    dx_post_sd
    dx_post_q025
    dx_post_q975

    p_dx_gt0
    p_dx_lt0.

Endpoint detectability is summarized per transition as:

    p_detect_mean
        = mean(p_detect_t0, p_detect_t1)

and:

    p_dropout_mean
        = mean(p_dropout_t0, p_dropout_t1).

These posterior quantities are DIAGNOSTICS.

They are not treated as alternative displacement representations and they do
not redefine the primary fluctuation estimator.


Temporal-lag stratification
---------------------------
Temporal intervals are always analyzed separately.

Selected lags can be supplied with:

    --dt 1 2 3 4 5

for example.

If `--dt` is omitted, every available dt is analyzed.

Curves never mix distinct dt values.


Shared absolute abundance range
-------------------------------
A common abundance range is constructed across the selected dt values.

For each lag d, let:

    xmin_d
        = minimum x_condition at dt=d

    xmax_d
        = maximum x_condition at dt=d.


Intersection mode
-----------------
With the default:

    --shared-range intersection

the common range is:

    x_min
        = max_d xmin_d

    x_max
        = min_d xmax_d.

This retains only abundance values represented at every selected temporal lag.


Union mode
----------
With:

    --shared-range union

the range is:

    x_min
        = min_d xmin_d

    x_max
        = max_d xmax_d.

This spans the pooled abundance support but does not guarantee that every bin is
represented at every dt.


Explicit range
--------------
The lower and/or upper limits can be overridden using:

    --x-min
    --x-max.


Shared absolute bins
--------------------
The common interval:

    [x_min, x_max]

is divided into:

    K = --n-bins

equal-width absolute-abundance bins.

The default is:

    K = 30.

Bin width is:

                 x_max - x_min
    h = -----------------------------.
                         K

For bin k:

    bin_left_k
        = x_min + k*h

    bin_right_k
        = x_min + (k+1)*h

    bin_mid_k
        = x_min + (k+0.5)*h.

The same absolute bins are used for every selected dt.

This is essential for abundance-matched comparison of fluctuation statistics
across temporal lags.


Bin-retention criteria
----------------------
After point statistics are calculated, a dt x abundance bin is retained only
when:

    n >= --min-n

and:

    n_subjects >= --min-subjects.

Defaults are:

    min_n = 30
    min_subjects = 2.


Primary point estimates
-----------------------
For temporal lag d and abundance bin b, let:

    Delta x_i

denote the retained transition displacements.


1. Mean displacement
--------------------
    mean_dx
        = (1/n) * sum_i Delta x_i.


2. Median displacement
----------------------
    median_dx
        = median_i(Delta x_i).


3. Positive-displacement probability
------------------------------------
    ppos
        = (1/n) * sum_i I(Delta x_i > 0).


4. Mean-squared displacement
----------------------------
    msd
        = (1/n) * sum_i (Delta x_i)^2.

This is the finite-time second moment:

    E[(Delta x)^2 | x_condition, dt].


5. Conditional displacement variance
------------------------------------
The primary fluctuation metric is the sample variance:

                           1
    var_dx = ----------------------------- *
                         n - 1

             sum_i (Delta x_i - mean_dx)^2.

The implementation uses:

    ddof = 1.

This measures dispersion around the conditional mean and therefore separates
fluctuation magnitude from the directional mean component.


6. Median absolute deviation
----------------------------
Let:

    m = median_dx.

Then:

    mad_dx
        = median_i |Delta x_i - m|.

This is retained as a robust descriptive dispersion measure.


Directional diagnostics
-----------------------
The primary fluctuation analysis is accompanied by:

    mean_dx
    median_dx
    ppos.

These are not alternative fluctuation definitions.

They document directional asymmetry within each abundance/lag bin and make the
distinction between variance and MSD interpretable.


Per-bin posterior/detectability summaries
-----------------------------------------
For each retained dt x abundance bin, the script additionally reports:

    median_p_detect
    median_p_dropout

    median_dx_post_sd

    median_dx_post_95_width
        = median(dx_post_q975 - dx_post_q025)

    median_p_dx_gt0
    median_p_dx_lt0.

These characterize inferential support within the same bins used for the
fluctuation analysis.


Primary estimator: transition weighted
--------------------------------------
The manuscript-primary pooled point estimates are transition weighted.

Each transition contributes equally within a dt x abundance bin.

Consequently, subjects contributing more retained transitions contribute more
weight to the pooled estimator.

This is intentional for the primary analysis and is explicitly evaluated
against equal-subject weighting as a sensitivity analysis.


Subject-bin sufficient summaries
---------------------------------
For each:

    dt x bin x subject

the script calculates:

    n_subject_bin
    sum_dx
    sum_dx2
    n_pos

    subject_mean_dx
    subject_median_dx
    subject_msd
    subject_var_dx
    subject_ppos
    subject_mad_dx.

These summaries support the subject-cluster bootstrap, equal-subject
sensitivity and leave-one-subject-out analyses.


Subject-cluster bootstrap
-------------------------
The primary uncertainty analysis resamples SUBJECTS with replacement within
each dt x abundance bin.

For a bootstrap sample, subject-specific sufficient statistics are pooled using
their bootstrap multiplicities.

For the transition-weighted estimator:

    N*
        = sum_s n_s*

    S1*
        = sum_s sum_dx_s*

    S2*
        = sum_s sum_dx2_s*

    P*
        = sum_s n_pos_s*.


Bootstrap mean
--------------
    mean_dx*
        = S1* / N*.


Bootstrap positive probability
------------------------------
    ppos*
        = P* / N*.


Bootstrap MSD
-------------
    msd*
        = S2* / N*.


Bootstrap variance
------------------
For N* > 1:

                          S2* - N*(mean_dx*)^2
    var_dx* = -----------------------------------------.
                                   N* - 1


Confidence intervals
--------------------
The requested confidence level is controlled by:

    --ci

with default:

    95.

For the default 95% interval, the 2.5th and 97.5th percentiles of the
subject-cluster bootstrap distribution are reported for:

    mean_dx
    ppos
    msd
    var_dx.

The default number of bootstrap replicates is:

    --n-bootstrap 2000

with seed:

    --seed 123.

Transition-weighted bootstrap intervals are not calculated for:

    median_dx
    mad_dx.


Equal-subject-weighted sensitivity
----------------------------------
When sensitivity analysis is enabled, Step 11 additionally constructs curves
in which each contributing subject receives equal weight within each
dt x abundance bin.

For example:

    mean_dx_equal
        = mean_s(subject_mean_dx_s)

and:

    msd_equal
        = mean_s(subject_msd_s).

The equal-subject variance is calculated as:

    var_dx_equal
        = msd_equal - mean_dx_equal^2.

This is a marginal second-moment decomposition of the equal-subject estimator;
it is not the arithmetic mean of the subject-specific variances.

Equal-subject weighting is a sensitivity analysis, not the manuscript-primary
estimator.


Weighting comparison
--------------------
For each point metric, the script compares:

    equal_subject_weighted_value

with:

    transition_weighted_value.

The signed difference is:

    difference_equal_minus_transition
        = equal_subject_value - transition_weighted_value.

When the transition-weighted value is sufficiently far from zero, the script
also reports:

    relative_difference
        = difference_equal_minus_transition
          / |transition_weighted_value|.


Leave-one-subject-out sensitivity
---------------------------------
For every subject, Step 11 recomputes the principal inferential metrics after
omitting that subject.

LOSO is generated for:

    transition_weighted
    equal_subject_weighted

and for the four primary bootstrap metrics:

    mean_dx
    ppos
    msd
    var_dx.

For each omission the script records:

    full_value
    loo_value

    signed_change
        = loo_value - full_value

    abs_change

    relative_change
        = signed_change / |full_value|

when the denominator is non-negligible.


LOSO summary
------------
Across subject omissions, the summary includes:

    n_loo_runs
    n_contributing_subjects_tested

    loo_mean
    loo_sd
    loo_min
    loo_max

    median_abs_change
    max_abs_change
    max_abs_relative_change

    most_influential_subject.

These are robustness diagnostics and do not modify the primary curves.


Temporal-comparability diagnostics
----------------------------------
Because the number of available subjects, intervals and clonotypes can differ
with dt, Step 11 explicitly characterizes the information structure available
at each temporal lag.


By dt x abundance bin
---------------------
The script reports:

    n_transitions
    n_subjects
    n_unique_clonotypes
    n_subject_clonotype_pairs
    n_subject_interval_pairs
    n_interval_pairs

    mean/sd/median/min/max x_condition

    mean/sd/min/max effective_dt.


By dt
-----
The same structural quantities are summarized across each temporal lag.


Subject availability by dt
--------------------------
For each dt x subject:

    n_transitions
    n_interval_pairs
    n_unique_clonotypes
    mean_effective_dt

are reported.

These diagnostics are descriptive.

They do not redefine the shared bins and are not used to reweight the primary
transition-weighted estimator.


Summary by temporal lag
-----------------------
For each dt, the retained curves are summarized by:

    n_bins_retained
    total_n_in_bins
    median_bin_n

    min_subjects_per_bin
    median_subjects_per_bin

    amplitude_mean_dx
        = max_b(mean_dx_b) - min_b(mean_dx_b)

    mean_abs_mean_dx

    mean_msd_across_bins
    mean_var_dx_across_bins
    mean_mad_dx_across_bins

    median_p_detect_across_bins
    median_p_dropout_across_bins.


Shared-bin coverage across dt
-----------------------------
For each abundance bin, the script reports:

    n_dt_present
    dt_present

    min_n_across_dt
    min_subjects_across_dt

    bin_left
    bin_right
    bin_mid

    present_in_all_selected_dt.

This table identifies which retained abundance bins remain empirically
supported across the full set of selected temporal lags.


Filtering and interpretation policy
-----------------------------------
The manuscript-primary Step-11 architecture is frozen to:

    representation = latent
    conditioning   = xstar
    mode           = TT.

Changing representation or conditioning requires explicit:

    --allow-nonprimary.

This prevents accidental substitution of a sensitivity analysis for the primary
estimand.

Changing the observability mode does not require this flag but should be
reported explicitly because it changes the transition estimand.


Outputs
-------
Unless disabled, tabular outputs are written in BOTH CSV and Parquet form.


00_analysis_metadata.json
-------------------------
Complete analysis metadata including:

    input arguments;
    script version;
    resolved input columns;
    selected dt values;
    shared abundance range;
    primary architecture;
    primary estimator;
    bootstrap unit;
    temporal-comparability outputs;
    sensitivity estimators;
    row counts.


00_run_config.json
------------------
Pipeline-wide copy of the same immutable run metadata.


01_input_schema.csv / .parquet
------------------------------
Input transition-table schema:

    column
    dtype.


02_filter_counts_by_dt_obs_class.csv / .parquet
------------------------------------------------
Counts after representation, dt and observability filtering, stratified by:

    dt
    obs_class

with transition, subject and subject-clonotype counts.


03_x_condition_range_by_dt.csv / .parquet
-----------------------------------------
Native conditioning-abundance range for each selected dt:

    x_min_dt
    x_max_dt
    n_transitions
    n_subjects.


04_shared_absolute_bin_definitions.csv / .parquet
--------------------------------------------------
Shared equal-width absolute conditioning bins:

    bin_id
    bin_left
    bin_right
    bin_mid.


05_filtered_binned_transitions.parquet
--------------------------------------
Optional compact transition-level table after filtering and shared-range
assignment.

Written only with:

    --write-filtered-transitions.


06_temporal_comparability_by_dt_bin.csv / .parquet
---------------------------------------------------
Sampling/composition comparability by:

    dt x abundance bin.


07_temporal_comparability_by_dt.csv / .parquet
-----------------------------------------------
Sampling/composition comparability summarized by temporal lag.


08_subject_availability_by_dt.csv / .parquet
--------------------------------------------
Subject-specific structural availability at each temporal lag.


09_subject_bin_dynamics.csv / .parquet
--------------------------------------
Subject-by-dt-by-bin sufficient and descriptive statistics used for bootstrap
and sensitivity analyses.


10_binned_dynamics_long.csv / .parquet
--------------------------------------
PRIMARY Step-11 abundance-resolved fluctuation table.

One row represents one retained:

    dt x abundance bin

and includes:

    representation
    conditioning
    mode

    n
    n_subjects

    mean_dx
    bootstrap CI

    median_dx

    ppos
    bootstrap CI

    msd
    bootstrap CI

    var_dx
    bootstrap CI

    mad_dx

    detectability/dropout diagnostics

    posterior transition-uncertainty diagnostics.


11_summary_by_dt.csv / .parquet
--------------------------------
Compact summary of retained fluctuation curves by temporal lag.


12_shared_bin_coverage_across_dt.csv / .parquet
------------------------------------------------
Coverage of each absolute abundance bin across the selected dt values.


13_weighting_sensitivity_curves_long.csv / .parquet
----------------------------------------------------
Primary transition-weighted and equal-subject-weighted curves in one common
schema.

Generated when sensitivity analysis is enabled.


14_weighting_sensitivity_comparison_long.csv / .parquet
--------------------------------------------------------
Point-by-point comparison of equal-subject and transition-weighted metrics.

Generated when sensitivity analysis is enabled.


15_leave_one_subject_out_long.csv / .parquet
---------------------------------------------
All subject omissions for both weighting estimators and the four primary
inferential metrics.

Generated when sensitivity analysis is enabled.


16_leave_one_subject_out_summary.csv / .parquet
------------------------------------------------
LOSO robustness summary, including the most influential subject in each
dt x abundance-bin x estimator x metric combination.

Generated when sensitivity analysis is enabled.


DEBUG_REPORT.txt
----------------
Human-readable run summary recording:

    input;
    dataset label;
    script version;
    representation;
    conditioning;
    fluctuation metrics;
    effective-dt source;
    selected dt;
    shared abundance range;
    bin criteria;
    bootstrap settings;
    row counts;
    sensitivity status.

The filename is retained for backward compatibility even though the content is
a reproducibility/run report rather than a debugging analysis.


Output-format controls
----------------------
CSV output is enabled by default and can be disabled with:

    --no-write-csv.

Parquet output is enabled by default and can be disabled with:

    --no-write-parquet.

Parquet compression defaults to:

    zstd.

Large transition input can be collected using Polars streaming mode with:

    --streaming.


Primary run
-----------
A standard manuscript-primary Step-11 run is:

    python3 11-longitudinal_fluctuation_dynamics.py \
        --transitions \
        ./results/5-latent_transition_construction/latent_transitions.parquet \
        --dataset-label healthy \
        --outdir \
        ./results/11-longitudinal_fluctuation_dynamics \
        --representation latent \
        --conditioning xstar \
        --dt 1 2 3 4 5 \
        --mode TT \
        --shared-range intersection \
        --n-bins 30 \
        --min-n 30 \
        --min-subjects 2 \
        --n-bootstrap 2000 \
        --ci 95 \
        --seed 123 \
        --streaming


Interpretation
--------------
Step 11 asks:

    "At a fixed absolute latent abundance, how large are genuine longitudinal
     clonotype fluctuations, and how do their conditional variance, second
     moment and robust dispersion change across temporal lags?"

The primary fluctuation quantity is:

    Var(dx_latent | xstar_latent, dt).

The complementary MSD is reported separately because forward drift is non-zero.

The analysis preserves distinct dt values and a common absolute abundance grid,
allowing subsequent temporal-scaling analysis to ask whether fluctuation
magnitude systematically accumulates with elapsed time.

Step 11 itself does not fit that temporal law.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import polars as pl


SCRIPT_VERSION = "v1-step11-latent-xstar-longitudinal-fluctuation-dynamics-2026-08-24"
PRIMARY_REPRESENTATION = "latent"
PRIMARY_CONDITIONING = "xstar"
PRIMARY_MODE = "TT"

PRIMARY_BOOT_METRICS = ("mean_dx", "ppos", "msd", "var_dx")
ALL_POINT_METRICS = (
    "mean_dx",
    "median_dx",
    "ppos",
    "msd",
    "var_dx",
    "mad_dx",
)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 11: analyze Step-5 clonotype transitions using one explicit representation "
            "and shared absolute abundance bins across temporal intervals."
        ),
    )

    ap.add_argument("--transitions", required=True)
    ap.add_argument("--dataset-label", default="dataset")
    ap.add_argument("--outdir", "--out-dir", dest="out_dir", required=True)

    ap.add_argument(
        "--representation",
        choices=["latent", "observed", "count"],
        default="latent",
        help=(
            "One representation per run. Posterior is not a representation: "
            "posterior fields are summarized as uncertainty descriptors."
        ),
    )
    ap.add_argument(
        "--conditioning",
        choices=["xstar", "x0", "midpoint"],
        default=PRIMARY_CONDITIONING,
        help=(
            "Primary conditioning is xstar. x0/midpoint are sensitivity-only and require "
            "--allow-nonprimary."
        ),
    )
    ap.add_argument(
        "--allow-nonprimary",
        action="store_true",
        help=(
            "Explicitly permit non-primary representation or conditioning settings. "
            "The manuscript-primary Step-11 run leaves this flag unset."
        ),
    )
    ap.add_argument(
        "--mode",
        choices=["TT", "ALL", "TF", "FT", "FF", "NON_TT", "NON_FF"],
        default=PRIMARY_MODE,
    )
    ap.add_argument(
        "--dt",
        type=int,
        nargs="*",
        default=None,
        help="Selected temporal intervals. Omit to analyze every available dt separately.",
    )

    ap.add_argument("--n-bins", type=int, default=30)
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument(
        "--min-subjects",
        type=int,
        default=2,
        help="Minimum subjects required to retain a dt × abundance bin.",
    )
    ap.add_argument(
        "--shared-range",
        choices=["intersection", "union"],
        default="intersection",
        help=(
            "Range used for common absolute bins. Intersection retains abundance "
            "regions represented in every selected dt; union spans the full pooled range."
        ),
    )
    ap.add_argument(
        "--x-min",
        type=float,
        default=None,
        help="Optional explicit lower bound for shared absolute bins.",
    )
    ap.add_argument(
        "--x-max",
        type=float,
        default=None,
        help="Optional explicit upper bound for shared absolute bins.",
    )

    ap.add_argument("--n-bootstrap", "--n-boot", dest="n_boot", type=int, default=2000)
    ap.add_argument("--ci", type=float, default=95.0)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument(
        "--sensitivity-analysis",
        action="store_true",
        default=True,
        help=(
            "Generate transition-weighted versus equal-subject-weighted curves "
            "and leave-one-subject-out diagnostics."
        ),
    )
    ap.add_argument(
        "--no-sensitivity-analysis",
        action="store_false",
        dest="sensitivity_analysis",
        help="Disable weighting and leave-one-subject-out sensitivity outputs.",
    )

    ap.add_argument("--streaming", action="store_true")
    ap.add_argument("--write-parquet", action="store_true", default=True)
    ap.add_argument("--no-write-parquet", action="store_false", dest="write_parquet")
    ap.add_argument("--write-csv", action="store_true", default=True)
    ap.add_argument("--no-write-csv", action="store_false", dest="write_csv")
    ap.add_argument(
        "--write-filtered-transitions",
        action="store_true",
        help="Write the compact filtered transition table. Disabled by default.",
    )
    ap.add_argument("--compression", default="zstd")

    return ap


# -----------------------------------------------------------------------------
# I/O and schema
# -----------------------------------------------------------------------------

def collect(lf: pl.LazyFrame, streaming: bool) -> pl.DataFrame:
    if streaming:
        try:
            return lf.collect(engine="streaming")
        except TypeError:
            return lf.collect(streaming=True)
    return lf.collect()


def scan_table(path: Path) -> pl.LazyFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pl.scan_parquet(str(path))
    if suffix == ".csv":
        return pl.scan_csv(str(path), infer_schema_length=10000)
    if suffix in {".tsv", ".txt"}:
        return pl.scan_csv(str(path), separator="\t", infer_schema_length=10000)
    raise ValueError("Unsupported input format: {}".format(path))


def _csv_safe_dataframe(df: pl.DataFrame) -> pl.DataFrame:
    """
    Return a CSV-compatible DataFrame.

    Step 11 output tables are expected to be flat. Nested columns are resolved
    where they are created rather than cast generically here.
    """
    nested = []

    for name, dtype in df.schema.items():
        dtype_name = str(dtype)
        if (
            dtype_name.startswith("List")
            or dtype_name.startswith("Array")
            or dtype_name.startswith("Struct")
            or dtype_name == "Object"
        ):
            nested.append((name, dtype_name))

    if nested:
        raise ValueError(
            "CSV output still contains nested columns: {}. "
            "Resolve them before write_table().".format(nested)
        )

    return df


def write_table(
    df: pl.DataFrame,
    path_no_suffix: Path,
    write_csv: bool,
    write_parquet: bool,
    compression: str,
) -> None:
    path_no_suffix.parent.mkdir(parents=True, exist_ok=True)

    if write_csv:
        csv_df = _csv_safe_dataframe(df)
        csv_df.write_csv(path_no_suffix.with_suffix(".csv"))

    if write_parquet:
        df.write_parquet(
            path_no_suffix.with_suffix(".parquet"),
            compression=compression,
        )


def first_existing(names: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    for candidate in candidates:
        if candidate in names:
            return candidate
    return None


def normalize_bool_expr(column: str) -> pl.Expr:
    value = pl.col(column).cast(pl.String, strict=False).str.strip_chars().str.to_uppercase()
    return (
        pl.when(value.is_in(["T", "TRUE", "1", "YES", "Y", "PRESENT", "DETECTED"]))
        .then(pl.lit(True))
        .when(value.is_in(["F", "FALSE", "0", "NO", "N", "ABSENT"]))
        .then(pl.lit(False))
        .otherwise(pl.lit(None, dtype=pl.Boolean))
    )


def representation_columns(
    schema_names: Sequence[str],
    representation: str,
    conditioning: str,
) -> Dict[str, Optional[str]]:
    catalog = {
        "latent": {
            "x0": ["x0_latent"],
            "x1": ["x1_latent"],
            "dx": ["dx_latent"],
            "xmid": ["xmid_latent"],
            "xstar": ["xstar_latent"],
        },
        "observed": {
            "x0": ["x0_obs", "x0_freq_obs", "x0_observed"],
            "x1": ["x1_obs", "x1_freq_obs", "x1_observed"],
            "dx": ["dx_obs", "dx_freq_obs", "dx_observed"],
            "xmid": ["xmid_obs", "xmid_freq_obs", "xmid_observed"],
            "xstar": ["xstar_obs", "xstar_freq_obs", "xstar_observed"],
        },
        "count": {
            "x0": ["x0_count", "x_count_sum_t0", "x0_counts", "x0_log_count"],
            "x1": ["x1_count", "x_count_sum_t1", "x1_counts", "x1_log_count"],
            "dx": ["dx_count", "dx_count_sum", "dx_counts", "dx_log_count"],
            "xmid": ["xmid_count", "xmid_count_sum", "xmid_counts", "xmid_log_count"],
            "xstar": ["xstar_count", "xstar_count_sum", "xstar_counts", "xstar_log_count"],
        },
    }

    spec = catalog[representation]
    out = {
        key: first_existing(schema_names, candidates)
        for key, candidates in spec.items()
    }

    if out["x0"] is None or out["x1"] is None:
        raise ValueError(
            "Representation '{}' is unavailable: x0={} x1={}".format(
                representation, out["x0"], out["x1"]
            )
        )

    if conditioning == "xstar" and out["xstar"] is None:
        raise ValueError(
            "Representation '{}' has no representation-specific xstar column. "
            "No latent-xstar fallback is applied automatically.".format(representation)
        )

    return out


# -----------------------------------------------------------------------------
# Input preparation
# -----------------------------------------------------------------------------

def prepare_transitions(
    lf: pl.LazyFrame,
    schema_names: Sequence[str],
    args: argparse.Namespace,
) -> Tuple[pl.DataFrame, Dict[str, Optional[str]], pl.DataFrame]:
    required_core = ["subject", "aaSeqCDR3", "t0", "t1"]
    missing = [c for c in required_core if c not in schema_names]
    if missing:
        raise ValueError("Missing required core columns: {}".format(missing))

    colmap = representation_columns(
        schema_names=schema_names,
        representation=args.representation,
        conditioning=args.conditioning,
    )

    dt_expr = (
        pl.col("dt").cast(pl.Int32, strict=False)
        if "dt" in schema_names
        else (
            pl.col("t1").cast(pl.Int32, strict=False)
            - pl.col("t0").cast(pl.Int32, strict=False)
        )
    )

    effective_dt_source = first_existing(
        schema_names,
        ["effective_dt", "dt_effective", "delta_t_effective"],
    )
    effective_dt_expr = (
        pl.col(effective_dt_source).cast(pl.Float64, strict=False)
        if effective_dt_source is not None
        else dt_expr.cast(pl.Float64)
    )

    if "obs_class" in schema_names:
        obs_class_expr = pl.col("obs_class").cast(pl.String).str.to_uppercase()
    elif "obs0" in schema_names and "obs1" in schema_names:
        obs0 = normalize_bool_expr("obs0")
        obs1 = normalize_bool_expr("obs1")
        obs_class_expr = (
            pl.when(obs0 & obs1).then(pl.lit("TT"))
            .when(obs0 & ~obs1).then(pl.lit("TF"))
            .when(~obs0 & obs1).then(pl.lit("FT"))
            .otherwise(pl.lit("FF"))
        )
    else:
        raise ValueError("Input requires obs_class or both obs0 and obs1.")

    x0_expr = pl.col(colmap["x0"]).cast(pl.Float64, strict=False)
    x1_expr = pl.col(colmap["x1"]).cast(pl.Float64, strict=False)

    if colmap["dx"] is not None:
        dx_expr = pl.col(colmap["dx"]).cast(pl.Float64, strict=False)
        dx_source = colmap["dx"]
    else:
        dx_expr = x1_expr - x0_expr
        dx_source = "__computed_x1_minus_x0"

    if colmap["xmid"] is not None:
        xmid_expr = pl.col(colmap["xmid"]).cast(pl.Float64, strict=False)
        xmid_source = colmap["xmid"]
    else:
        xmid_expr = 0.5 * (x0_expr + x1_expr)
        xmid_source = "__computed_midpoint"

    if args.conditioning == "xstar":
        x_condition_expr = pl.col(colmap["xstar"]).cast(pl.Float64, strict=False)
        x_condition_source = colmap["xstar"]
    elif args.conditioning == "x0":
        x_condition_expr = x0_expr
        x_condition_source = colmap["x0"]
    else:
        x_condition_expr = xmid_expr
        x_condition_source = xmid_source

    optional_exprs = []

    # Current Step-5 detectability/dropout and posterior-transition fields.
    for alias, source in [
        ("p_detect_t0", "p_detect_state_t0"),
        ("p_detect_t1", "p_detect_state_t1"),
        ("p_dropout_t0", "p_dropout_state_t0"),
        ("p_dropout_t1", "p_dropout_state_t1"),
        ("dx_post_sd", "dx_post_sd"),
        ("dx_post_q025", "dx_post_q025"),
        ("dx_post_q975", "dx_post_q975"),
        ("p_dx_gt0", "p_dx_gt0"),
        ("p_dx_lt0", "p_dx_lt0"),
    ]:
        if source in schema_names:
            optional_exprs.append(
                pl.col(source).cast(pl.Float64, strict=False).alias(alias)
            )
        else:
            optional_exprs.append(pl.lit(None, dtype=pl.Float64).alias(alias))

    selected = lf.select(
        [
            pl.col("subject").cast(pl.String, strict=False).alias("subject"),
            pl.col("aaSeqCDR3").cast(pl.String).alias("aaSeqCDR3"),
            pl.col("t0").cast(pl.Int32, strict=False).alias("t0"),
            pl.col("t1").cast(pl.Int32, strict=False).alias("t1"),
            dt_expr.alias("dt"),
            effective_dt_expr.alias("effective_dt"),
            obs_class_expr.alias("obs_class"),
            x0_expr.alias("x0"),
            x1_expr.alias("x1"),
            xmid_expr.alias("xmid"),
            x_condition_expr.alias("x_condition"),
            dx_expr.alias("dx"),
            *optional_exprs,
        ]
    ).with_columns(
        [
            pl.mean_horizontal(["p_detect_t0", "p_detect_t1"]).alias(
                "p_detect_mean"
            ),
            pl.mean_horizontal(["p_dropout_t0", "p_dropout_t1"]).alias(
                "p_dropout_mean"
            ),
        ]
    )

    # Core validity filter before collection.
    selected = selected.filter(
        pl.col("subject").is_not_null()
        & pl.col("t0").is_not_null()
        & pl.col("t1").is_not_null()
        & pl.col("dt").is_not_null()
        & pl.col("x_condition").is_not_null()
        & pl.col("x_condition").is_finite()
        & pl.col("dx").is_not_null()
        & pl.col("dx").is_finite()
    )

    if args.dt:
        selected = selected.filter(pl.col("dt").is_in(args.dt))

    mode = args.mode.upper()
    if mode == "TT":
        selected = selected.filter(pl.col("obs_class") == "TT")
    elif mode in {"TF", "FT", "FF"}:
        selected = selected.filter(pl.col("obs_class") == mode)
    elif mode == "NON_TT":
        selected = selected.filter(pl.col("obs_class") != "TT")
    elif mode == "NON_FF":
        # Detectability-inclusive continuous-dynamics sensitivity:
        # retain TT + TF + FT while excluding FF transitions.
        selected = selected.filter(pl.col("obs_class") != "FF")
    elif mode != "ALL":
        raise ValueError("Unsupported mode: {}".format(args.mode))

    # Counts are computed from the already standardized/filterable lazy frame.
    filter_counts = collect(
        selected.group_by(["dt", "obs_class"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.struct(["subject", "aaSeqCDR3"])
                .n_unique()
                .alias("n_subject_clonotype_pairs"),
            ]
        )
        .sort(["dt", "obs_class"]),
        args.streaming,
    )

    data = collect(selected, args.streaming)

    resolved = {
        "representation": args.representation,
        "x0": colmap["x0"],
        "x1": colmap["x1"],
        "dx": dx_source,
        "xmid": xmid_source,
        "conditioning": args.conditioning,
        "x_condition": x_condition_source,
        "effective_dt": effective_dt_source or "__nominal_dt",
    }

    return data, resolved, filter_counts


# -----------------------------------------------------------------------------
# Shared bins
# -----------------------------------------------------------------------------

def determine_shared_range(
    data: pl.DataFrame,
    shared_range: str,
    explicit_min: Optional[float],
    explicit_max: Optional[float],
) -> Tuple[float, float, pl.DataFrame]:
    dt_ranges = (
        data.group_by("dt")
        .agg(
            [
                pl.col("x_condition").min().alias("x_min_dt"),
                pl.col("x_condition").max().alias("x_max_dt"),
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
            ]
        )
        .sort("dt")
    )

    if explicit_min is not None:
        x_min = float(explicit_min)
    elif shared_range == "intersection":
        x_min = float(dt_ranges["x_min_dt"].max())
    else:
        x_min = float(dt_ranges["x_min_dt"].min())

    if explicit_max is not None:
        x_max = float(explicit_max)
    elif shared_range == "intersection":
        x_max = float(dt_ranges["x_max_dt"].min())
    else:
        x_max = float(dt_ranges["x_max_dt"].max())

    if not math.isfinite(x_min) or not math.isfinite(x_max) or x_max <= x_min:
        raise ValueError(
            "Invalid shared abundance range: x_min={} x_max={}. "
            "Try --shared-range union or explicit --x-min/--x-max.".format(
                x_min, x_max
            )
        )

    return x_min, x_max, dt_ranges


def assign_shared_bins(
    data: pl.DataFrame,
    x_min: float,
    x_max: float,
    n_bins: int,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    if n_bins < 2:
        raise ValueError("--n-bins must be >= 2")

    width = (x_max - x_min) / float(n_bins)
    filtered = data.filter(
        (pl.col("x_condition") >= x_min)
        & (pl.col("x_condition") <= x_max)
    ).with_columns(
        (
            ((pl.col("x_condition") - x_min) / width)
            .floor()
            .clip(0, n_bins - 1)
            .cast(pl.Int32)
        ).alias("bin_id")
    ).with_columns(
        [
            (pl.lit(x_min) + pl.col("bin_id") * width).alias("bin_left"),
            (pl.lit(x_min) + (pl.col("bin_id") + 1) * width).alias("bin_right"),
            (pl.lit(x_min) + (pl.col("bin_id") + 0.5) * width).alias(
                "bin_mid"
            ),
        ]
    )

    bin_table = pl.DataFrame(
        {
            "bin_id": list(range(n_bins)),
            "bin_left": [x_min + i * width for i in range(n_bins)],
            "bin_right": [x_min + (i + 1) * width for i in range(n_bins)],
            "bin_mid": [x_min + (i + 0.5) * width for i in range(n_bins)],
        }
    )

    return filtered, bin_table


# -----------------------------------------------------------------------------
# Metrics and cluster bootstrap
# -----------------------------------------------------------------------------

def add_point_metrics(data: pl.DataFrame) -> Tuple[pl.DataFrame, pl.DataFrame]:
    keys = ["dt", "bin_id", "bin_left", "bin_right", "bin_mid"]

    base = (
        data.group_by(keys)
        .agg(
            [
                pl.len().alias("n"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.col("x_condition").min().alias("x_raw_min"),
                pl.col("x_condition").median().alias("x_raw_median"),
                pl.col("x_condition").max().alias("x_raw_max"),
                pl.col("dx").mean().alias("mean_dx"),
                pl.col("dx").median().alias("median_dx"),
                (pl.col("dx") > 0).mean().alias("ppos"),
                (pl.col("dx") ** 2).mean().alias("msd"),
                pl.col("dx").var(ddof=1).alias("var_dx"),
                pl.col("p_detect_mean").median().alias("median_p_detect"),
                pl.col("p_dropout_mean").median().alias("median_p_dropout"),
                pl.col("dx_post_sd").median().alias("median_dx_post_sd"),
                (
                    pl.col("dx_post_q975") - pl.col("dx_post_q025")
                ).median().alias("median_dx_post_95_width"),
                pl.col("p_dx_gt0").median().alias("median_p_dx_gt0"),
                pl.col("p_dx_lt0").median().alias("median_p_dx_lt0"),
            ]
        )
        .sort(["dt", "bin_id"])
    )

    medians = base.select(["dt", "bin_id", "median_dx"])
    mad = (
        data.join(medians, on=["dt", "bin_id"], how="left")
        .with_columns(
            (pl.col("dx") - pl.col("median_dx")).abs().alias("_abs_dev")
        )
        .group_by(["dt", "bin_id"])
        .agg(pl.col("_abs_dev").median().alias("mad_dx"))
    )

    base = base.join(mad, on=["dt", "bin_id"], how="left")

    # Subject-bin summaries are reusable for cluster bootstrap, equal-subject
    # weighting, and leave-one-subject-out diagnostics.
    subject_keys = ["dt", "bin_id", "subject"]
    subject = (
        data.group_by(subject_keys)
        .agg(
            [
                pl.len().alias("n_subject_bin"),
                pl.col("dx").sum().alias("sum_dx"),
                (pl.col("dx") ** 2).sum().alias("sum_dx2"),
                (pl.col("dx") > 0).sum().alias("n_pos"),
                pl.col("dx").mean().alias("subject_mean_dx"),
                pl.col("dx").median().alias("subject_median_dx"),
                (pl.col("dx") ** 2).mean().alias("subject_msd"),
                pl.col("dx").var(ddof=1).alias("subject_var_dx"),
                (pl.col("dx") > 0).mean().alias("subject_ppos"),
            ]
        )
        .sort(subject_keys)
    )

    subject_medians = subject.select(
        subject_keys + ["subject_median_dx"]
    )
    subject_mad = (
        data.join(subject_medians, on=subject_keys, how="left")
        .with_columns(
            (pl.col("dx") - pl.col("subject_median_dx"))
            .abs()
            .alias("_subject_abs_dev")
        )
        .group_by(subject_keys)
        .agg(
            pl.col("_subject_abs_dev")
            .median()
            .alias("subject_mad_dx")
        )
    )
    subject = (
        subject.join(subject_mad, on=subject_keys, how="left")
        .sort(subject_keys)
    )

    return base, subject


def bootstrap_subject_clusters(
    subject_stats: pl.DataFrame,
    n_boot: int,
    ci: float,
    seed: int,
) -> pl.DataFrame:
    columns = [
        "dt",
        "bin_id",
        "mean_dx_ci_low",
        "mean_dx_ci_high",
        "ppos_ci_low",
        "ppos_ci_high",
        "msd_ci_low",
        "msd_ci_high",
        "var_dx_ci_low",
        "var_dx_ci_high",
        "bootstrap_n_subjects",
        "bootstrap_replicates",
    ]

    if n_boot <= 0 or subject_stats.height == 0:
        return pl.DataFrame(schema={c: pl.Float64 for c in columns})

    rng = np.random.default_rng(seed)
    rows = []

    for group in subject_stats.partition_by(["dt", "bin_id"], maintain_order=True):
        dt = int(group["dt"][0])
        bin_id = int(group["bin_id"][0])
        n_subjects = group.height

        if n_subjects < 2:
            rows.append(
                {
                    "dt": dt,
                    "bin_id": bin_id,
                    "mean_dx_ci_low": None,
                    "mean_dx_ci_high": None,
                    "ppos_ci_low": None,
                    "ppos_ci_high": None,
                    "msd_ci_low": None,
                    "msd_ci_high": None,
                    "var_dx_ci_low": None,
                    "var_dx_ci_high": None,
                    "bootstrap_n_subjects": n_subjects,
                    "bootstrap_replicates": n_boot,
                }
            )
            continue

        n = group["n_subject_bin"].to_numpy().astype(float)
        sum_dx = group["sum_dx"].to_numpy().astype(float)
        sum_dx2 = group["sum_dx2"].to_numpy().astype(float)
        n_pos = group["n_pos"].to_numpy().astype(float)

        boot_mean = np.empty(n_boot, dtype=float)
        boot_ppos = np.empty(n_boot, dtype=float)
        boot_msd = np.empty(n_boot, dtype=float)
        boot_var = np.empty(n_boot, dtype=float)

        for b in range(n_boot):
            idx = rng.integers(0, n_subjects, size=n_subjects)

            total_n = float(n[idx].sum())
            total_sum = float(sum_dx[idx].sum())
            total_sum2 = float(sum_dx2[idx].sum())
            total_pos = float(n_pos[idx].sum())

            mean_dx = total_sum / total_n
            msd = total_sum2 / total_n
            variance = (
                (total_sum2 - total_n * mean_dx * mean_dx) / (total_n - 1.0)
                if total_n > 1
                else np.nan
            )

            boot_mean[b] = mean_dx
            boot_ppos[b] = total_pos / total_n
            boot_msd[b] = msd
            boot_var[b] = variance

        alpha = (100.0 - ci) / 2.0

        def interval(values: np.ndarray) -> Tuple[float, float]:
            return (
                float(np.nanpercentile(values, alpha)),
                float(np.nanpercentile(values, 100.0 - alpha)),
            )

        mean_lo, mean_hi = interval(boot_mean)
        ppos_lo, ppos_hi = interval(boot_ppos)
        msd_lo, msd_hi = interval(boot_msd)
        var_lo, var_hi = interval(boot_var)

        rows.append(
            {
                "dt": dt,
                "bin_id": bin_id,
                "mean_dx_ci_low": mean_lo,
                "mean_dx_ci_high": mean_hi,
                "ppos_ci_low": ppos_lo,
                "ppos_ci_high": ppos_hi,
                "msd_ci_low": msd_lo,
                "msd_ci_high": msd_hi,
                "var_dx_ci_low": var_lo,
                "var_dx_ci_high": var_hi,
                "bootstrap_n_subjects": n_subjects,
                "bootstrap_replicates": n_boot,
            }
        )

    return pl.DataFrame(rows).sort(["dt", "bin_id"])


def bootstrap_equal_subject_curves(
    subject_stats: pl.DataFrame,
    n_boot: int,
    ci: float,
    seed: int,
) -> pl.DataFrame:
    """Bootstrap equal-subject-weighted metrics within each dt × abundance bin."""
    metric_columns = {
        "mean_dx": "subject_mean_dx",
        "median_dx": "subject_median_dx",
        "ppos": "subject_ppos",
        "msd": "subject_msd",
        "mad_dx": "subject_mad_dx",
    }
    output_columns = ["dt", "bin_id"]
    for metric in ALL_POINT_METRICS:
        output_columns.extend(
            ["{}_ci_low".format(metric), "{}_ci_high".format(metric)]
        )
    output_columns.extend(["bootstrap_n_subjects", "bootstrap_replicates"])

    if n_boot <= 0 or subject_stats.height == 0:
        return pl.DataFrame(schema={c: pl.Float64 for c in output_columns})

    rng = np.random.default_rng(seed)
    alpha = (100.0 - ci) / 2.0
    rows = []

    for group in subject_stats.partition_by(["dt", "bin_id"], maintain_order=True):
        dt = int(group["dt"][0])
        bin_id = int(group["bin_id"][0])
        n_subjects = group.height
        row = {
            "dt": dt,
            "bin_id": bin_id,
            "bootstrap_n_subjects": n_subjects,
            "bootstrap_replicates": n_boot,
        }

        if n_subjects < 2:
            for metric in ALL_POINT_METRICS:
                row["{}_ci_low".format(metric)] = None
                row["{}_ci_high".format(metric)] = None
            rows.append(row)
            continue

        bootstrap_indices = rng.integers(
            0, n_subjects, size=(n_boot, n_subjects)
        )
        for metric, source in metric_columns.items():
            values = group[source].to_numpy().astype(float)
            sampled = values[bootstrap_indices]
            if metric in {"median_dx", "mad_dx"}:
                estimates = np.nanmedian(sampled, axis=1)
            else:
                valid_counts = np.isfinite(sampled).sum(axis=1)
                sampled_sums = np.nansum(sampled, axis=1)
                estimates = np.divide(
                    sampled_sums,
                    valid_counts,
                    out=np.full(n_boot, np.nan, dtype=float),
                    where=valid_counts > 0,
                )
            finite = estimates[np.isfinite(estimates)]
            if finite.size == 0:
                low = high = None
            else:
                low = float(np.percentile(finite, alpha))
                high = float(np.percentile(finite, 100.0 - alpha))
            row["{}_ci_low".format(metric)] = low
            row["{}_ci_high".format(metric)] = high

        sampled_mean = group["subject_mean_dx"].to_numpy().astype(float)[
            bootstrap_indices
        ]
        sampled_msd = group["subject_msd"].to_numpy().astype(float)[
            bootstrap_indices
        ]
        mean_counts = np.isfinite(sampled_mean).sum(axis=1)
        msd_counts = np.isfinite(sampled_msd).sum(axis=1)
        boot_equal_mean = np.divide(
            np.nansum(sampled_mean, axis=1),
            mean_counts,
            out=np.full(n_boot, np.nan, dtype=float),
            where=mean_counts > 0,
        )
        boot_equal_msd = np.divide(
            np.nansum(sampled_msd, axis=1),
            msd_counts,
            out=np.full(n_boot, np.nan, dtype=float),
            where=msd_counts > 0,
        )
        boot_equal_var = boot_equal_msd - boot_equal_mean ** 2
        finite_var = boot_equal_var[np.isfinite(boot_equal_var)]
        if finite_var.size == 0:
            row["var_dx_ci_low"] = None
            row["var_dx_ci_high"] = None
        else:
            row["var_dx_ci_low"] = float(np.percentile(finite_var, alpha))
            row["var_dx_ci_high"] = float(
                np.percentile(finite_var, 100.0 - alpha)
            )

        rows.append(row)

    return pl.DataFrame(rows).sort(["dt", "bin_id"])


def build_equal_subject_curves(
    subject_stats: pl.DataFrame,
    point_reference: pl.DataFrame,
    n_boot: int,
    ci: float,
    seed: int,
) -> pl.DataFrame:
    """Construct curves in which each contributing subject has equal weight."""
    equal = (
        subject_stats.group_by(["dt", "bin_id"])
        .agg(
            [
                pl.col("n_subject_bin").sum().alias("n"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.col("subject_mean_dx").mean().alias("mean_dx"),
                pl.col("subject_median_dx").median().alias("median_dx"),
                pl.col("subject_ppos").mean().alias("ppos"),
                pl.col("subject_msd").mean().alias("msd"),
                pl.col("subject_mad_dx").median().alias("mad_dx"),
            ]
        )
        .with_columns(
            (pl.col("msd") - pl.col("mean_dx") ** 2).alias("var_dx")
        )
        .sort(["dt", "bin_id"])
    )

    reference_columns = [
        "dt",
        "bin_id",
        "bin_left",
        "bin_mid",
        "bin_right",
        "x_raw_min",
        "x_raw_median",
        "x_raw_max",
    ]
    equal = equal.join(
        point_reference.select(reference_columns),
        on=["dt", "bin_id"],
        how="left",
    )

    boot = bootstrap_equal_subject_curves(
        subject_stats=subject_stats,
        n_boot=n_boot,
        ci=ci,
        seed=seed,
    )
    return equal.join(boot, on=["dt", "bin_id"], how="left").sort(
        ["dt", "bin_id"]
    )


def assemble_weighting_sensitivity_curves(
    transition_curves: pl.DataFrame,
    equal_curves: pl.DataFrame,
) -> pl.DataFrame:
    """Return both estimators in one schema-compatible sensitivity table."""
    identity_columns = [
        "dataset_label",
        "representation",
        "conditioning",
        "abundance_coordinate",
        "mode",
    ]
    base_columns = [
        "dt",
        "bin_id",
        "bin_left",
        "bin_mid",
        "bin_right",
        "x_raw_min",
        "x_raw_median",
        "x_raw_max",
        "n",
        "n_subjects",
    ]
    ci_columns = []
    for metric in ALL_POINT_METRICS:
        ci_columns.extend(
            ["{}_ci_low".format(metric), "{}_ci_high".format(metric)]
        )
    output_columns = (
        identity_columns
        + ["estimator"]
        + base_columns
        + list(ALL_POINT_METRICS)
        + ci_columns
        + ["bootstrap_n_subjects", "bootstrap_replicates"]
    )

    transition = transition_curves
    for metric in ("median_dx", "mad_dx"):
        for suffix in ("ci_low", "ci_high"):
            column = "{}_{}".format(metric, suffix)
            if column not in transition.columns:
                transition = transition.with_columns(
                    pl.lit(None, dtype=pl.Float64).alias(column)
                )
    transition = transition.with_columns(
        pl.lit("transition_weighted").alias("estimator")
    ).select(output_columns)

    identity_values = {
        column: transition_curves[column][0] for column in identity_columns
    }
    equal = equal_curves.with_columns(
        [
            pl.lit(identity_values[column]).alias(column)
            for column in identity_columns
        ]
        + [pl.lit("equal_subject_weighted").alias("estimator")]
    ).select(output_columns)

    return pl.concat([transition, equal], how="vertical").sort(
        ["estimator", "dt", "bin_id"]
    )


def build_weighting_comparison(
    weighting_curves: pl.DataFrame,
) -> pl.DataFrame:
    """Compare equal-subject and transition-weighted point estimates."""
    keys = [
        "dataset_label",
        "representation",
        "conditioning",
        "abundance_coordinate",
        "mode",
        "dt",
        "bin_id",
        "bin_left",
        "bin_mid",
        "bin_right",
        "n",
        "n_subjects",
    ]
    transition = weighting_curves.filter(
        pl.col("estimator") == "transition_weighted"
    )
    equal = weighting_curves.filter(
        pl.col("estimator") == "equal_subject_weighted"
    )

    rows = []
    for metric in ALL_POINT_METRICS:
        joined = transition.select(
            keys + [pl.col(metric).alias("transition_weighted_value")]
        ).join(
            equal.select(
                ["dt", "bin_id", pl.col(metric).alias("equal_subject_value")]
            ),
            on=["dt", "bin_id"],
            how="inner",
        )
        joined = joined.with_columns(
            [
                pl.lit(metric).alias("metric"),
                (
                    pl.col("equal_subject_value")
                    - pl.col("transition_weighted_value")
                ).alias("difference_equal_minus_transition"),
                pl.when(pl.col("transition_weighted_value").abs() > 1e-12)
                .then(
                    (
                        pl.col("equal_subject_value")
                        - pl.col("transition_weighted_value")
                    )
                    / pl.col("transition_weighted_value").abs()
                )
                .otherwise(None)
                .alias("relative_difference"),
            ]
        )
        rows.append(joined)

    if not rows:
        return pl.DataFrame()
    return pl.concat(rows, how="vertical").sort(["metric", "dt", "bin_id"])


def _transition_weighted_metrics_from_subject_arrays(
    n: np.ndarray,
    sum_dx: np.ndarray,
    sum_dx2: np.ndarray,
    n_pos: np.ndarray,
) -> Dict[str, float]:
    total_n = float(np.nansum(n))
    if total_n <= 0:
        return {metric: np.nan for metric in PRIMARY_BOOT_METRICS}
    total_sum = float(np.nansum(sum_dx))
    total_sum2 = float(np.nansum(sum_dx2))
    mean_dx = total_sum / total_n
    variance = (
        (total_sum2 - total_n * mean_dx * mean_dx) / (total_n - 1.0)
        if total_n > 1.0
        else np.nan
    )
    return {
        "mean_dx": mean_dx,
        "ppos": float(np.nansum(n_pos)) / total_n,
        "msd": total_sum2 / total_n,
        "var_dx": variance,
    }


def _equal_subject_metrics_from_subject_arrays(
    subject_mean: np.ndarray,
    subject_ppos: np.ndarray,
    subject_msd: np.ndarray,
) -> Dict[str, float]:
    arrays = {
        "mean_dx": subject_mean,
        "ppos": subject_ppos,
        "msd": subject_msd,
    }
    output = {}
    for metric, values in arrays.items():
        finite = values[np.isfinite(values)]
        output[metric] = float(np.mean(finite)) if finite.size else np.nan
    if np.isfinite(output["mean_dx"]) and np.isfinite(output["msd"]):
        output["var_dx"] = output["msd"] - output["mean_dx"] ** 2
    else:
        output["var_dx"] = np.nan
    return output


def build_leave_one_subject_out(
    subject_stats: pl.DataFrame,
    point_reference: pl.DataFrame,
    all_subjects: Sequence[str],
    dataset_label: str,
) -> pl.DataFrame:
    """Compute leave-one-subject-out curves for the four inferential metrics."""
    reference = {
        (int(row["dt"]), int(row["bin_id"])): row
        for row in point_reference.select(
            ["dt", "bin_id", "bin_left", "bin_mid", "bin_right"]
        ).iter_rows(named=True)
    }
    rows = []

    for group in subject_stats.partition_by(["dt", "bin_id"], maintain_order=True):
        dt = int(group["dt"][0])
        bin_id = int(group["bin_id"][0])
        bin_info = reference[(dt, bin_id)]
        subjects = group["subject"].to_numpy().astype(str)
        n = group["n_subject_bin"].to_numpy().astype(float)
        sum_dx = group["sum_dx"].to_numpy().astype(float)
        sum_dx2 = group["sum_dx2"].to_numpy().astype(float)
        n_pos = group["n_pos"].to_numpy().astype(float)
        subject_mean = group["subject_mean_dx"].to_numpy().astype(float)
        subject_ppos = group["subject_ppos"].to_numpy().astype(float)
        subject_msd = group["subject_msd"].to_numpy().astype(float)

        full_by_estimator = {
            "transition_weighted": _transition_weighted_metrics_from_subject_arrays(
                n, sum_dx, sum_dx2, n_pos
            ),
            "equal_subject_weighted": _equal_subject_metrics_from_subject_arrays(
                subject_mean, subject_ppos, subject_msd
            ),
        }

        for omitted_subject in all_subjects:
            keep = subjects != str(omitted_subject)
            contributed = bool(np.any(~keep))
            n_subjects_remaining = int(np.sum(keep))
            n_transitions_remaining = int(np.nansum(n[keep]))
            loo_by_estimator = {
                "transition_weighted": _transition_weighted_metrics_from_subject_arrays(
                    n[keep], sum_dx[keep], sum_dx2[keep], n_pos[keep]
                ),
                "equal_subject_weighted": _equal_subject_metrics_from_subject_arrays(
                    subject_mean[keep],
                    subject_ppos[keep],
                    subject_msd[keep],
                ),
            }

            for estimator in (
                "transition_weighted",
                "equal_subject_weighted",
            ):
                for metric in PRIMARY_BOOT_METRICS:
                    full_value = full_by_estimator[estimator][metric]
                    loo_value = loo_by_estimator[estimator][metric]
                    signed_change = loo_value - full_value
                    relative_change = (
                        signed_change / abs(full_value)
                        if np.isfinite(full_value) and abs(full_value) > 1e-12
                        else np.nan
                    )
                    rows.append(
                        {
                            "dataset_label": dataset_label,
                            "estimator": estimator,
                            "metric": metric,
                            "dt": dt,
                            "bin_id": bin_id,
                            "bin_left": float(bin_info["bin_left"]),
                            "bin_mid": float(bin_info["bin_mid"]),
                            "bin_right": float(bin_info["bin_right"]),
                            "omitted_subject": str(omitted_subject),
                            "omitted_subject_contributed": contributed,
                            "n_subjects_full": int(group.height),
                            "n_subjects_remaining": n_subjects_remaining,
                            "n_transitions_remaining": n_transitions_remaining,
                            "full_value": full_value,
                            "loo_value": loo_value,
                            "signed_change": signed_change,
                            "abs_change": abs(signed_change),
                            "relative_change": relative_change,
                        }
                    )

    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows).sort(
        ["estimator", "metric", "dt", "bin_id", "omitted_subject"]
    )


def summarize_leave_one_subject_out(loo: pl.DataFrame) -> pl.DataFrame:
    if loo.height == 0:
        return pl.DataFrame()
    numeric_columns = [
        "full_value",
        "loo_value",
        "signed_change",
        "abs_change",
        "relative_change",
    ]
    loo = loo.with_columns(
        [pl.col(column).fill_nan(None).alias(column) for column in numeric_columns]
    )
    keys = [
        "dataset_label",
        "estimator",
        "metric",
        "dt",
        "bin_id",
        "bin_left",
        "bin_mid",
        "bin_right",
        "n_subjects_full",
        "full_value",
    ]
    return (
        loo.group_by(keys)
        .agg(
            [
                pl.len().alias("n_loo_runs"),
                pl.col("omitted_subject_contributed")
                .sum()
                .alias("n_contributing_subjects_tested"),
                pl.col("loo_value").mean().alias("loo_mean"),
                pl.col("loo_value").std(ddof=1).alias("loo_sd"),
                pl.col("loo_value").min().alias("loo_min"),
                pl.col("loo_value").max().alias("loo_max"),
                pl.col("abs_change").median().alias("median_abs_change"),
                pl.col("abs_change").max().alias("max_abs_change"),
                pl.col("relative_change")
                .abs()
                .max()
                .alias("max_abs_relative_change"),
                pl.col("omitted_subject")
                .sort_by("abs_change", descending=True)
                .first()
                .alias("most_influential_subject"),
            ]
        )
        .sort(["estimator", "metric", "dt", "bin_id"])
    )



# -----------------------------------------------------------------------------
# Temporal comparability diagnostics
# -----------------------------------------------------------------------------


def build_temporal_comparability_by_dt_bin(data: pl.DataFrame) -> pl.DataFrame:
    """
    Summarize sampling/composition comparability for each dt x shared abundance bin.

    The table is deliberately descriptive and is not used to redefine bins or
    reweight the primary estimator.
    """
    if data.height == 0:
        return pl.DataFrame()

    return (
        data.group_by(["dt", "bin_id", "bin_left", "bin_right", "bin_mid"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes"),
                pl.struct(["subject", "aaSeqCDR3"])
                .n_unique()
                .alias("n_subject_clonotype_pairs"),
                pl.struct(["subject", "t0", "t1"])
                .n_unique()
                .alias("n_subject_interval_pairs"),
                pl.struct(["t0", "t1"]).n_unique().alias("n_interval_pairs"),
                pl.col("x_condition").mean().alias("mean_x_condition"),
                pl.col("x_condition").std(ddof=1).alias("sd_x_condition"),
                pl.col("x_condition").median().alias("median_x_condition"),
                pl.col("x_condition").min().alias("min_x_condition"),
                pl.col("x_condition").max().alias("max_x_condition"),
                pl.col("effective_dt").mean().alias("mean_effective_dt"),
                pl.col("effective_dt").std(ddof=1).alias("sd_effective_dt"),
                pl.col("effective_dt").min().alias("min_effective_dt"),
                pl.col("effective_dt").max().alias("max_effective_dt"),
            ]
        )
        .sort(["dt", "bin_id"])
    )


def build_temporal_comparability_by_dt(data: pl.DataFrame) -> pl.DataFrame:
    """Summarize the structurally available longitudinal information at each dt."""
    if data.height == 0:
        return pl.DataFrame()

    return (
        data.group_by("dt")
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects_available"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes"),
                pl.struct(["subject", "aaSeqCDR3"])
                .n_unique()
                .alias("n_subject_clonotype_pairs"),
                pl.struct(["subject", "t0", "t1"])
                .n_unique()
                .alias("n_subject_interval_pairs"),
                pl.struct(["t0", "t1"]).n_unique().alias("n_interval_pairs"),
                pl.col("x_condition").mean().alias("mean_x_condition"),
                pl.col("x_condition").std(ddof=1).alias("sd_x_condition"),
                pl.col("x_condition").median().alias("median_x_condition"),
                pl.col("x_condition").min().alias("min_x_condition"),
                pl.col("x_condition").max().alias("max_x_condition"),
                pl.col("effective_dt").mean().alias("mean_effective_dt"),
                pl.col("effective_dt").std(ddof=1).alias("sd_effective_dt"),
                pl.col("effective_dt").min().alias("min_effective_dt"),
                pl.col("effective_dt").max().alias("max_effective_dt"),
            ]
        )
        .sort("dt")
    )


def build_subject_availability_by_dt(data: pl.DataFrame) -> pl.DataFrame:
    """List the subjects structurally contributing at each temporal lag."""
    if data.height == 0:
        return pl.DataFrame()

    return (
        data.group_by(["dt", "subject"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.struct(["t0", "t1"]).n_unique().alias("n_interval_pairs"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes"),
                pl.col("effective_dt").mean().alias("mean_effective_dt"),
            ]
        )
        .sort(["dt", "subject"])
    )


# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

def build_curve_summary(curves: pl.DataFrame) -> pl.DataFrame:
    if curves.height == 0:
        return pl.DataFrame()

    return (
        curves.group_by("dt")
        .agg(
            [
                pl.len().alias("n_bins_retained"),
                pl.col("n").sum().alias("total_n_in_bins"),
                pl.col("n").median().alias("median_bin_n"),
                pl.col("n_subjects").min().alias("min_subjects_per_bin"),
                pl.col("n_subjects").median().alias("median_subjects_per_bin"),
                (
                    pl.col("mean_dx").max() - pl.col("mean_dx").min()
                ).alias("amplitude_mean_dx"),
                pl.col("mean_dx").abs().mean().alias("mean_abs_mean_dx"),
                pl.col("msd").mean().alias("mean_msd_across_bins"),
                pl.col("var_dx").mean().alias("mean_var_dx_across_bins"),
                pl.col("mad_dx").mean().alias("mean_mad_dx_across_bins"),
                pl.col("median_p_detect").median().alias(
                    "median_p_detect_across_bins"
                ),
                pl.col("median_p_dropout").median().alias(
                    "median_p_dropout_across_bins"
                ),
            ]
        )
        .sort("dt")
    )


def build_bin_coverage(curves: pl.DataFrame, selected_dts: Sequence[int]) -> pl.DataFrame:
    if curves.height == 0:
        return pl.DataFrame()

    n_dt = len(selected_dts)

    coverage = (
        curves.group_by("bin_id")
        .agg(
            [
                pl.col("dt").n_unique().alias("n_dt_present"),
                pl.col("dt").sort().alias("_dt_present_list"),
                pl.col("n").min().alias("min_n_across_dt"),
                pl.col("n_subjects").min().alias("min_subjects_across_dt"),
                pl.col("bin_left").first().alias("bin_left"),
                pl.col("bin_right").first().alias("bin_right"),
                pl.col("bin_mid").first().alias("bin_mid"),
            ]
        )
        .with_columns(
            [
                pl.col("_dt_present_list")
                .list.eval(pl.element().cast(pl.String))
                .list.join(",")
                .alias("dt_present"),
                (pl.col("n_dt_present") == n_dt).alias(
                    "present_in_all_selected_dt"
                ),
            ]
        )
        .drop("_dt_present_list")
        .sort("bin_id")
    )

    return coverage


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    args = build_argparser().parse_args()

    nonprimary = []
    if args.representation != PRIMARY_REPRESENTATION:
        nonprimary.append("representation={}".format(args.representation))
    if args.conditioning != PRIMARY_CONDITIONING:
        nonprimary.append("conditioning={}".format(args.conditioning))
    if nonprimary and not args.allow_nonprimary:
        raise ValueError(
            "Step 11 primary architecture is frozen to representation='{}' and "
            "conditioning='{}'. Non-primary setting(s): {}. Use --allow-nonprimary "
            "only for an explicitly labelled sensitivity analysis.".format(
                PRIMARY_REPRESENTATION,
                PRIMARY_CONDITIONING,
                ", ".join(nonprimary),
            )
        )

    input_path = Path(args.transitions)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[INFO] 11-longitudinal_fluctuation_dynamics.py | {}".format(SCRIPT_VERSION))
    print("[INFO] Primary fluctuation estimand: Var(dx_latent | xstar_latent, dt)")
    print("[INFO] Complementary second moment: MSD = E[dx_latent^2 | xstar_latent, dt]")
    print("[INFO] Representation={} | conditioning={} | mode={}".format(
        args.representation, args.conditioning, args.mode
    ))

    if not input_path.exists():
        raise FileNotFoundError("Input not found: {}".format(input_path))

    lf = scan_table(input_path)
    schema = lf.collect_schema()
    schema_names = schema.names()

    schema_table = pl.DataFrame(
        {
            "column": schema_names,
            "dtype": [str(schema[name]) for name in schema_names],
        }
    )
    write_table(
        schema_table,
        out_dir / "01_input_schema",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    print("[READ] {}".format(input_path))
    data, resolved_columns, filter_counts = prepare_transitions(
        lf, schema_names, args
    )

    if data.height == 0:
        raise ValueError("No transitions remain after filtering.")

    selected_dts = sorted(int(v) for v in data["dt"].unique().to_list())
    print("[FILTER] rows={:,}; dt={}".format(data.height, selected_dts))

    write_table(
        filter_counts,
        out_dir / "02_filter_counts_by_dt_obs_class",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    x_min, x_max, dt_ranges = determine_shared_range(
        data=data,
        shared_range=args.shared_range,
        explicit_min=args.x_min,
        explicit_max=args.x_max,
    )
    write_table(
        dt_ranges,
        out_dir / "03_x_condition_range_by_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    binned, bin_table = assign_shared_bins(
        data=data,
        x_min=x_min,
        x_max=x_max,
        n_bins=args.n_bins,
    )
    write_table(
        bin_table,
        out_dir / "04_shared_absolute_bin_definitions",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    if args.write_filtered_transitions:
        write_table(
            binned,
            out_dir / "05_filtered_binned_transitions",
            write_csv=False,
            write_parquet=True,
            compression=args.compression,
        )

    comparability_dt_bin = build_temporal_comparability_by_dt_bin(binned)
    write_table(
        comparability_dt_bin,
        out_dir / "06_temporal_comparability_by_dt_bin",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    comparability_dt = build_temporal_comparability_by_dt(binned)
    write_table(
        comparability_dt,
        out_dir / "07_temporal_comparability_by_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    subject_availability = build_subject_availability_by_dt(binned)
    write_table(
        subject_availability,
        out_dir / "08_subject_availability_by_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    point, subject = add_point_metrics(binned)

    point = point.filter(
        (pl.col("n") >= args.min_n)
        & (pl.col("n_subjects") >= args.min_subjects)
    )

    # Keep subject rows only for retained bins.
    retained_keys = point.select(["dt", "bin_id"])
    subject = subject.join(retained_keys, on=["dt", "bin_id"], how="inner")

    boot = bootstrap_subject_clusters(
        subject_stats=subject,
        n_boot=args.n_boot,
        ci=args.ci,
        seed=args.seed,
    )

    curves = (
        point.join(boot, on=["dt", "bin_id"], how="left")
        .with_columns(
            [
                pl.lit(args.dataset_label).alias("dataset_label"),
                pl.lit(args.representation).alias("representation"),
                pl.lit(args.conditioning).alias("conditioning"),
                pl.lit("absolute_shared").alias("abundance_coordinate"),
                pl.lit(args.mode).alias("mode"),
            ]
        )
        .select(
            [
                "dataset_label",
                "representation",
                "conditioning",
                "abundance_coordinate",
                "mode",
                "dt",
                "bin_id",
                "bin_left",
                "bin_mid",
                "bin_right",
                "x_raw_min",
                "x_raw_median",
                "x_raw_max",
                "n",
                "n_subjects",
                "mean_dx",
                "mean_dx_ci_low",
                "mean_dx_ci_high",
                "median_dx",
                "ppos",
                "ppos_ci_low",
                "ppos_ci_high",
                "msd",
                "msd_ci_low",
                "msd_ci_high",
                "var_dx",
                "var_dx_ci_low",
                "var_dx_ci_high",
                "mad_dx",
                "median_p_detect",
                "median_p_dropout",
                "median_dx_post_sd",
                "median_dx_post_95_width",
                "median_p_dx_gt0",
                "median_p_dx_lt0",
                "bootstrap_n_subjects",
                "bootstrap_replicates",
            ]
        )
        .sort(["dt", "bin_id"])
    )

    write_table(
        subject,
        out_dir / "09_subject_bin_dynamics",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )
    write_table(
        curves,
        out_dir / "10_binned_dynamics_long",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    summary = build_curve_summary(curves)
    write_table(
        summary,
        out_dir / "11_summary_by_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    coverage = build_bin_coverage(curves, selected_dts)
    write_table(
        coverage,
        out_dir / "12_shared_bin_coverage_across_dt",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    sensitivity_curves = pl.DataFrame()
    weighting_comparison = pl.DataFrame()
    loo_long = pl.DataFrame()
    loo_summary = pl.DataFrame()

    if args.sensitivity_analysis:
        print("[SENSITIVITY] Building equal-subject-weighted curves")
        equal_curves = build_equal_subject_curves(
            subject_stats=subject,
            point_reference=point,
            n_boot=args.n_boot,
            ci=args.ci,
            seed=args.seed,
        )
        sensitivity_curves = assemble_weighting_sensitivity_curves(
            transition_curves=curves,
            equal_curves=equal_curves,
        )
        weighting_comparison = build_weighting_comparison(sensitivity_curves)

        all_subjects = sorted(
            str(value) for value in binned["subject"].unique().to_list()
        )
        print(
            "[SENSITIVITY] Leave-one-subject-out across {} subjects".format(
                len(all_subjects)
            )
        )
        loo_long = build_leave_one_subject_out(
            subject_stats=subject,
            point_reference=point,
            all_subjects=all_subjects,
            dataset_label=args.dataset_label,
        )
        loo_summary = summarize_leave_one_subject_out(loo_long)

        write_table(
            sensitivity_curves,
            out_dir / "13_weighting_sensitivity_curves_long",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )
        write_table(
            weighting_comparison,
            out_dir / "14_weighting_sensitivity_comparison_long",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )
        write_table(
            loo_long,
            out_dir / "15_leave_one_subject_out_long",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )
        write_table(
            loo_summary,
            out_dir / "16_leave_one_subject_out_summary",
            args.write_csv,
            args.write_parquet,
            args.compression,
        )

    metadata = vars(args).copy()
    metadata.update(
        {
            "script_version": SCRIPT_VERSION,
            "analysis_layer": "step11_longitudinal_fluctuation_dynamics",
            "primary_architecture": {
                "representation": PRIMARY_REPRESENTATION,
                "conditioning": PRIMARY_CONDITIONING,
                "mode": PRIMARY_MODE,
                "primary_fluctuation_metric": "var_dx",
                "complementary_second_moment": "msd",
                "robust_descriptive_dispersion": "mad_dx",
                "directional_diagnostics": ["mean_dx", "median_dx", "ppos"],
            },
            "is_primary_configuration": (
                args.representation == PRIMARY_REPRESENTATION
                and args.conditioning == PRIMARY_CONDITIONING
                and args.mode == PRIMARY_MODE
            ),
            "generated": datetime.now().isoformat(timespec="seconds"),
            "resolved_columns": resolved_columns,
            "selected_dt": selected_dts,
            "shared_x_min": x_min,
            "shared_x_max": x_max,
            "n_rows_after_filter": data.height,
            "n_rows_in_shared_range": binned.height,
            "n_curve_rows": curves.height,
            "temporal_comparability_outputs": [
                "06_temporal_comparability_by_dt_bin",
                "07_temporal_comparability_by_dt",
                "08_subject_availability_by_dt",
            ],
            "posterior_is_representation": False,
            "dt_mixed_in_curves": False,
            "bootstrap_unit": "subject_cluster",
            "primary_estimator": "transition_weighted",
            "sensitivity_estimators": [
                "transition_weighted",
                "equal_subject_weighted",
            ] if args.sensitivity_analysis else [],
            "leave_one_subject_out_metrics": list(PRIMARY_BOOT_METRICS)
            if args.sensitivity_analysis else [],
            "n_sensitivity_curve_rows": sensitivity_curves.height,
            "n_loo_rows": loo_long.height,
        }
    )
    with open(out_dir / "00_analysis_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)
    # Duplicate the same immutable run metadata under the pipeline-wide config name.
    with open(out_dir / "00_run_config.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)

    with open(out_dir / "DEBUG_REPORT.txt", "w", encoding="utf-8") as fh:
        fh.write("STEP 11: LONGITUDINAL FLUCTUATION DYNAMICS\n")
        fh.write("=========================================\n\n")
        fh.write("input: {}\n".format(input_path))
        fh.write("dataset_label: {}\n".format(args.dataset_label))
        fh.write("script_version: {}\n".format(SCRIPT_VERSION))
        fh.write("representation: {}\n".format(args.representation))
        fh.write("conditioning: {}\n".format(args.conditioning))
        fh.write("primary_fluctuation_metric: var_dx\n")
        fh.write("complementary_second_moment: msd\n")
        fh.write("robust_descriptive_dispersion: mad_dx\n")
        fh.write("effective_dt_source: {}\n".format(resolved_columns.get("effective_dt")))
        fh.write("mode: {}\n".format(args.mode))
        fh.write("selected_dt: {}\n".format(selected_dts))
        fh.write("shared_range_method: {}\n".format(args.shared_range))
        fh.write("shared_x_range: [{:.8g}, {:.8g}]\n".format(x_min, x_max))
        fh.write("n_bins_defined: {}\n".format(args.n_bins))
        fh.write("min_n: {}\n".format(args.min_n))
        fh.write("min_subjects: {}\n".format(args.min_subjects))
        fh.write("bootstrap_unit: subject_cluster\n")
        fh.write("n_boot: {}\n".format(args.n_boot))
        fh.write("rows_after_filter: {:,}\n".format(data.height))
        fh.write("rows_in_shared_range: {:,}\n".format(binned.height))
        fh.write("curve_rows: {:,}\n".format(curves.height))
        fh.write("sensitivity_analysis: {}\n".format(args.sensitivity_analysis))
        fh.write("sensitivity_curve_rows: {:,}\n".format(sensitivity_curves.height))
        fh.write("loo_rows: {:,}\n".format(loo_long.height))
        fh.write("Posterior quantities are summarized as latent-transition uncertainty.\n")
        fh.write("No figures and no cross-dataset comparisons are generated.\n")

    print("[DONE] Step 11 completed")
    print("[OUT] {}".format(out_dir))
    print("[MAIN] {}".format(out_dir / "10_binned_dynamics_long.parquet"))
    print("[SUBJECT] {}".format(out_dir / "09_subject_bin_dynamics.parquet"))
    if args.sensitivity_analysis:
        print("[SENSITIVITY] {}".format(
            out_dir / "13_weighting_sensitivity_curves_long.parquet"
        ))
        print("[LOO] {}".format(
            out_dir / "16_leave_one_subject_out_summary.parquet"
        ))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise