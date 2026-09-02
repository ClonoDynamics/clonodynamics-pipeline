#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
14-interval_position_structure.py
==================================

Interval-position and anchored-lag sensitivity analysis for ClonoDynamics.

Overview
--------
This script implements Step 14 of the ClonoDynamics pipeline.

It tests whether apparent temporal-lag structure in genuine longitudinal
clonotype fluctuations can be explained, in whole or in part, by the calendar
position of the contributing intervals rather than by elapsed lag alone.

The script is downstream of the finalized longitudinal fluctuation and
temporal-scaling architecture:

    Step 11
    11-longitudinal_fluctuation_dynamics.py

        -> latent displacement: dx_latent
        -> conditioning: xstar_latent
        -> transition class: TT
        -> shared absolute xstar bins

    Step 13
    13-temporal_fluctuation_scaling.py

        -> internally comparable longitudinal core bins
        -> temporal scaling across selected dt values

Step 14 then asks two complementary sensitivity questions:

    1. Calendar-position structure

       For a fixed nominal lag dt and abundance bin, do subject-level
       fluctuation metrics align systematically with the actual interval
       position (for example 1->2 versus 4->5)?

    2. Anchored lag structure

       Does a signed lag trend persist when interval composition is constrained
       by a common starting or ending time point?

The analysis is deliberately model-agnostic with respect to periodicity.

It does NOT fit a sinusoid or claim a periodic biological mechanism.

An exploratory alternating-position score is reported only as a descriptive
diagnostic when enough interval positions are available.


Position in the ClonoDynamics pipeline
--------------------------------------
Relevant public workflow:

    Step 5
    5-latent_transition_construction.py
                |
                v

    Step 11
    11-longitudinal_fluctuation_dynamics.py

        shared absolute xstar bins
        + abundance-resolved Var/MSD
                |
                v

    Step 13
    13-temporal_fluctuation_scaling.py

        selected internally comparable
        longitudinal core bins
                |
                v

    Step 14
    14-interval_position_structure.py

        Step-5 transition table
        + Step-11 bin geometry
        + Step-13 selected core bins
                |
                +--> interval-position permutation tests
                |
                +--> common-start anchored profiles
                |
                +--> common-end anchored profiles
                |
                +--> available-subject analysis
                |
                +--> matched-all-dt sensitivity


Step 14 does not replace Step 13 temporal-model analysis.

It is a composition-control and interval-position sensitivity analysis designed
to determine whether lag-dependent structure is robust to the calendar
placement and subject composition of contributing intervals.


Frozen primary architecture
---------------------------
The public primary configuration is:

    representation
        latent

    displacement
        dx_latent

    conditioning
        xstar_latent

    transition class
        TT

    primary fluctuation metric
        var_dx

    complementary fluctuation metric
        msd.

Other supported metrics may be requested explicitly as sensitivity analyses.


Inputs
------
Three inputs are required.


1. --transitions
----------------
The genuine longitudinal transition table generated in Step 5:

    latent_transitions.parquet

Required fields are:

    subject
    aaSeqCDR3
    t0
    t1
    xstar_latent
    dx_latent.

The script also uses:

    dt

when available, otherwise:

    dt = t1 - t0.

The transition class must be available either as:

    obs_class

or as the pair:

    obs0
    obs1.

Only TT transitions are retained.

If a `pattern` column is present, the script resolves one transition pattern
and avoids mixing duplicated all-pair/adjacent representations.


2. --step11-dir
----------------
The output directory from:

    11-longitudinal_fluctuation_dynamics.py.

Step 14 reads:

    00_analysis_metadata.json

and:

    04_shared_absolute_bin_definitions.csv/parquet.

The Step-11 metadata are used to verify the frozen latent/xstar/TT
configuration.

The shared absolute bin table supplies:

    bin_id
    bin_left
    bin_mid
    bin_right.


3. --step13-dir
----------------
The output directory from:

    13-temporal_fluctuation_scaling.py.

Step 14 reads:

    00_run_config.json

and:

    02_tested_abundance_bins.csv/parquet.

The selected Step-13 bins define the internally comparable longitudinal
abundance core used here.


Step-11 validation
------------------
Step 14 requires the Step-11 run to use:

    representation = latent
    conditioning   = xstar
    mode           = TT.

The resolved upstream displacement must be:

    dx_latent

or the equivalent explicitly reconstructed:

    x1_latent - x0_latent.

The resolved conditioning coordinate must be:

    xstar_latent.

The Step-11 public analysis layer is expected to be:

    step11_longitudinal_fluctuation_dynamics.

A non-primary Step-11 configuration causes Step 14 to stop.


Step-13 validation
------------------
The Step-13 configuration must confirm:

    upstream_representation_validated = True

    upstream_required_representation = latent

    upstream_required_conditioning = xstar

    upstream_required_x_condition_source = xstar_latent

    upstream_required_mode = TT

and:

    primary_metric = var_dx.

The Step-13 public analysis layer is expected to be:

    step13_temporal_fluctuation_scaling.

Step 14 additionally requires that Step 13 did NOT restrict temporal scaling to
the narrower pseudo-overlap abundance domain.


Core abundance bins
-------------------
By default, Step 14 reuses the bins marked:

    selected_for_analysis = True

in the Step-13:

    02_tested_abundance_bins

table.

An explicit subset can be supplied with:

    --core-bin-ids.

Explicit bin IDs must exist in the Step-11 shared-bin table.

The same Step-11 bin geometry is then used to assign transition-level xstar
values back to abundance bins.


Uniform-bin requirement
-----------------------
The current transition reassignment implementation requires the Step-11
absolute abundance bins to be:

    contiguous

and:

    uniformly spaced.

If the bin IDs are not contiguous or widths differ, execution stops rather
than silently assigning transitions incorrectly.


Transition preparation
----------------------
For every input transition, Step 14 standardizes:

    subject
    aaSeqCDR3
    t0
    t1
    dt
    effective_dt
    xstar_latent
    dx_latent
    obs_class

and optionally:

    pattern.

The effective temporal lag is read from the first available field among:

    effective_dt
    dt_effective
    delta_t_effective.

If none is present:

    effective_dt = dt.


Primary transition filtering
----------------------------
Rows are retained only when:

    dt is among the requested dt values;

    obs_class = TT;

    xstar_latent is finite;

    dx_latent is finite;

    xstar_latent lies within the Step-11 absolute-bin range;

    the assigned bin belongs to the Step-13 longitudinal core.

The default requested lags are:

    dt = 1,2,3,4,5.


Transition-pattern policy
-------------------------
When a `pattern` field exists:

    --pattern auto

prefers:

    all

when available.

If only one pattern is present, that pattern is used.

If multiple patterns are present and none is `all`, the user must specify the
desired pattern explicitly.

This prevents accidental duplication of the same biological transition through
multiple transition-construction modes.


Duplicate-transition guard
--------------------------
After filtering, Step 14 checks uniqueness of:

    subject
    aaSeqCDR3
    t0
    t1
    bin_id.

If duplicated keys remain, execution stops.

This protects the interval-position analysis from accidental repeated
all/adjacent transition representations or duplicated input rows.


Supported metrics
-----------------
The default metrics are:

    var_dx
    msd.

Additional supported sensitivity metrics are:

    mad_dx
    mean_dx
    median_dx
    ppos.

If `var_dx` is omitted, the script emits a warning because the run is no longer
the manuscript-primary Step-14 configuration.


Subject x interval x abundance-bin summaries
--------------------------------------------
For every:

    subject x dt x t0 x t1 x abundance bin

the script calculates:

    n_transitions

    n_unique_clonotypes

    mean_xstar_latent
    median_xstar_latent

    mean_effective_dt

    sum_dx
    sum_dx2

    mean_dx
    median_dx

    ppos
        = P(dx_latent > 0)

    msd
        = mean(dx_latent^2)

    var_dx
        = sample variance of dx_latent, ddof=1

    mad_dx
        = median |dx_latent - median(dx_latent)|.

The actual calendar interval is retained as:

    interval_label
        = t0 -> t1

and its calendar position is represented by:

    interval_position
        = t0.


Pooled interval summaries
-------------------------
For each:

    dt x t0 x t1 x abundance bin

Step 14 also reports transition-weighted pooled quantities including:

    n_transitions
    n_subjects

    n_unique_clonotypes_global
    n_subject_clonotype_pairs

    mean_xstar_latent
    median_xstar_latent

    mean_dx_transition_weighted
    median_dx_transition_weighted
    ppos_transition_weighted

    msd_transition_weighted
    var_dx_transition_weighted
    mad_dx_transition_weighted.


Equal-subject interval profiles
-------------------------------
Subject-level interval metrics are additionally aggregated with equal subject
weight.

For most metrics:

    cohort_value
        = mean_s(metric_s).

For:

    median_dx

the cohort summary uses the median across subject-specific median values.

The resulting profile is indexed by:

    dt
    t0
    t1
    abundance bin.

These equal-subject profiles form the basis of the interval-position
permutation test.


Primary interval-position question
----------------------------------
For fixed:

    dt
    abundance bin
    metric,

subject-level interval values are aligned by calendar position.

Let:

    M_s(j)

denote the metric for subject s at calendar interval position j.

Only interval positions with at least:

    --min-subjects-per-interval

contributing subjects are retained.

The default is:

    3.

At least:

    --min-interval-positions

positions are required for the primary permutation test.

The default is:

    3.


Observed cohort interval profile
--------------------------------
At each retained interval position j, the equal-subject cohort profile is:

                 1
    M_bar(j) = -------- sum_s M_s(j)
               n_j

over finite subject values.


Primary interval-position statistic
-----------------------------------
Let:

    M_bar_all
        = mean_j M_bar(j).

The primary synchronized-position statistic is:

    RMS_position
        = sqrt[
            mean_j (M_bar(j) - M_bar_all)^2
          ].

This measures the magnitude of calendar-position variation in the cohort
profile for a fixed lag and abundance bin.


Calendar slope
--------------
A linear calendar-position slope is also calculated:

    M_bar(j)
        = a + b * position_j.

The inferential statistic is:

    |b|.

This is descriptive of monotonic calendar-position structure and is separate
from the primary RMS statistic.


Exploratory alternating-position score
--------------------------------------
When at least four interval positions are available, the centered cohort
profile is compared with the alternating contrast:

    +1, -1, +1, -1, ...

after centering that contrast.

The score is:

    |corr(centered cohort profile, alternating contrast)|.

This is explicitly exploratory.

It is NOT a sinusoidal fit, a periodicity test or evidence of oscillatory
biology.


Centered sign changes
---------------------
For descriptive purposes, the script also counts sign changes in the centered
cohort profile after exact zero values are removed.


Permutation null
----------------
The null is generated by independently permuting interval labels WITHIN each
subject.

For one subject, only that subject's finite observed metric values are shuffled
among that subject's available calendar positions.

Therefore the permutation preserves:

    subject identity;

    nominal dt;

    abundance bin;

    the subject's observed metric-value set;

    the subject-specific missingness pattern;

    the number of available calendar positions.

It destroys only:

    cross-subject alignment of calendar interval position.


Scientific interpretation of the permutation
--------------------------------------------
The empirical null tests:

    synchronized calendar-position structure.

It does NOT test:

    generic interval heterogeneity

because the set of interval values within each subject is preserved.


Number of permutations
----------------------
The default number of interval-label permutations is:

    --n-perm 2000.

The global random seed is:

    --seed 123.


Permutation P-values
--------------------
For an observed non-negative statistic T_obs and permutation-null values T_j:

    p
        = [1 + sum_j I(T_j >= T_obs)] / (N + 1).

This upper-tail empirical P-value is used for:

    RMS position structure;

    absolute calendar slope;

    exploratory alternation score.


Permutation confidence descriptors
----------------------------------
For each null distribution the script reports:

    median
    lower CI quantile
    upper CI quantile

using the confidence level supplied by:

    --ci

with default:

    95%.


Multiple-testing correction
---------------------------
Benjamini-Hochberg correction is applied separately within each metric across
all valid:

    dt x abundance-bin

tests.

Adjusted q-values are reported for:

    RMS position structure;

    absolute calendar slope;

    alternation score.


Permutation-summary output
--------------------------
For each metric, Step 14 reports:

    n_tests
    n_dt_tested
    n_bins_tested

and, where available:

    number of BH-adjusted q values < 0.05

    minimum BH-adjusted q value

for each of the three interval-position statistics.


Anchored lag analysis
---------------------
Step 14 additionally constructs two interval-composition-constrained temporal
profiles.


1. Common-start profile
-----------------------
All transitions beginning at the earliest observed t0 are retained:

    t0 = min(t0).


2. Common-end profile
---------------------
All transitions ending at the latest observed t1 are retained:

    t1 = max(t1).

These analyses constrain calendar composition while allowing lag dt to vary.


Anchored signed temporal model
------------------------------
For each:

    anchor type
    metric
    abundance bin,

the equal-subject anchored profile is fit with:

    M(dt)
        = K + D * (dt - 1).

The fit is ordinary least squares.

The output reports:

    intercept_at_dt1 = K

and:

    signed_slope = D.

At least:

    --min-anchor-dt-points

finite lag values are required.

The default is:

    4.


Joint subject bootstrap for anchored slopes
-------------------------------------------
Subjects are resampled with replacement.

For each bootstrap draw, the same sampled-subject multiplicities are applied
jointly across all available dt values within the anchored profile.

This preserves subject-level covariance across lag.

The default number of anchor bootstrap replicates is:

    --n-boot-anchor 2000.


Available-subject anchored mode
-------------------------------
The first anchored analysis uses:

    subject_matching = available.

At each dt, all finite subject values available at that lag contribute.

Therefore the contributing subject set may differ across lag.

This reproduces the original availability-based anchor analysis.


Strict matched-all-dt mode
--------------------------
The second anchored analysis uses:

    subject_matching = matched_all_dt.

Within each:

    anchor type
    x metric
    x abundance bin,

a subject is retained only when that subject has a finite metric value at EVERY
requested dt.

The subject set is therefore identical across all temporal lags.


Minimum matched subjects
------------------------
The matched-all-dt inferential bootstrap is generated only when at least:

    --min-matched-anchor-subjects

complete-case subjects are available.

The default is:

    3.

Otherwise:

    test_status = insufficient_matched_subjects.


Anchored bootstrap summaries
----------------------------
For every anchored profile the script reports:

    observed_intercept_at_dt1
    observed_signed_slope

    n_bootstrap_requested
    n_valid_bootstrap

    signed_slope_median
    signed_slope_q_low
    signed_slope_q_high

    signed_slope_negative_fraction
    signed_slope_positive_fraction

    signed_slope_ci_below_zero
    signed_slope_ci_above_zero

    subject counts by dt

    contributing subject IDs

    observed profile values.


Anchor-subject matching diagnostic
----------------------------------
A separate structural diagnostic records, for every subject and anchor type:

    available dt values
    missing dt values
    complete_all_requested_dt.

This table is based on structural anchor availability.

The strict inferential matched-all-dt analysis applies the stronger criterion
of finite metric values within each metric x abundance bin.


Interpretive boundary
---------------------
Step 14 is a sensitivity analysis.

It should be interpreted as follows:

    - interval-position permutation tests ask whether calendar-position effects
      are synchronized across subjects;

    - anchored analyses ask whether signed lag trends persist when interval
      composition is constrained;

    - matched-all-dt anchored analyses additionally remove lag-dependent
      subject-composition differences;

    - the exploratory alternation score does not establish periodicity;

    - no sinusoidal model is fitted;

    - Step 14 does not replace the six-model temporal-scaling analysis of
      Step 13.


Outputs
-------
Unless disabled, tables are written in both CSV and Parquet form.


00_run_config.json
------------------
Complete Step-14 run configuration including:

    input paths;
    script version;
    resolved transition columns;
    frozen representation;
    selected core bins;
    core abundance range;
    dt values;
    interval-position permutation design;
    multiple-testing policy;
    anchor bootstrap design;
    available and matched-all-dt subject modes;
    validated Step-11 and Step-13 provenance.


01_input_manifest.csv / .parquet
--------------------------------
Input and design provenance including:

    Step-5 transition path;
    Step-11 fluctuation directory;
    Step-13 temporal-scaling directory;
    Step-11 shared-bin table;
    Step-13 tested-bin table;
    filtered transition count;
    subject count;
    selected core bins;
    resolved pattern and effective-dt source.


02_core_bin_definitions.csv / .parquet
--------------------------------------
The Step-11 shared absolute xstar bins retained by the Step-13 longitudinal
core.


03_interval_composition_by_pair.csv / .parquet
----------------------------------------------
Structural composition by:

    dt x t0 x t1.

Contains:

    n_transitions
    n_subjects
    n_unique_clonotypes_global
    n_subject_clonotype_pairs
    n_core_bins_present
    mean/sd effective dt
    mean/sd xstar
    interval label.


04_subject_interval_metrics.csv / .parquet
------------------------------------------
Subject x interval x abundance-bin metrics used by both permutation and
anchored analyses.


05_cohort_interval_profiles.csv / .parquet
------------------------------------------
Transition-weighted and equal-subject cohort interval profiles.


06_interval_position_permutation_tests.csv / .parquet
------------------------------------------------------
Full per:

    metric x dt x abundance-bin

interval-position permutation results, including observed statistics,
permutation-null summaries, empirical P-values and BH-adjusted q-values.


07_interval_position_test_summary.csv / .parquet
-------------------------------------------------
Compact metric-level summary of the calendar-position permutation tests.


08_anchored_temporal_profiles.csv / .parquet
--------------------------------------------
Common-start and common-end temporal profiles for:

    available
    matched_all_dt

subject-set modes.


09_anchored_signed_slope_bootstrap.csv / .parquet
--------------------------------------------------
Observed and joint-subject-bootstrap signed lag slopes for both anchors and both
subject-matching modes.


10_analysis_summary.csv / .parquet
----------------------------------
Compact integrated summary combining:

    interval-position permutation evidence

and:

    anchored signed-slope evidence.


11_anchor_subject_matching_diagnostic.csv / .parquet
-----------------------------------------------------
Subject-level structural availability across requested dt values for common
start and common end anchors.


DEBUG_REPORT.txt
----------------
Human-readable Step-14 run summary including:

    input paths;
    selected bins;
    core x range;
    filtered transition count;
    subject count;
    metrics;
    permutation and bootstrap settings;
    matched-subject availability;
    permutation interpretation;
    compact analysis summaries.


Primary run
-----------
A standard Step-14 run is:

    python3 14-interval_position_structure.py \
        --transitions \
        ./results/5-latent_transition_construction/latent_transitions.parquet \
        --step11-dir \
        ./results/11-longitudinal_fluctuation_dynamics \
        --step13-dir \
        ./results/13-temporal_fluctuation_scaling \
        --outdir \
        ./results/14-interval_position_structure \
        --dt-values 1,2,3,4,5 \
        --metrics var_dx,msd \
        --n-perm 2000 \
        --n-boot-anchor 2000 \
        --min-interval-positions 3 \
        --min-subjects-per-interval 3 \
        --min-anchor-dt-points 4 \
        --min-matched-anchor-subjects 3 \
        --seed 123 \
        --streaming


Interpretation
--------------
Step 14 asks:

    "Could the apparent temporal-lag structure of longitudinal clonotype
     fluctuations reflect which calendar intervals and subjects contribute at
     each lag?"

The answer is evaluated in two ways:

    synchronized calendar-position permutation tests;

and:

    common-start/common-end anchored lag analyses, including strict
    matched-subject sensitivity.

A lack of synchronized calendar-position structure supports the interpretation
that the temporal-lag result is not simply a shared calendar-position artifact.

Persistence of a signed lag trend in both anchored and matched-all-dt analyses
provides an additional composition-control check.

Conversely, material dependence on interval position, anchor choice or subject
matching should be reported as a limitation on the interpretation of the
Step-13 temporal-scaling result.
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
import polars as pl


SCRIPT_VERSION = "v2-step14-interval-position-structure-2026-08-25"
PRIMARY_METRIC = "var_dx"
COMPLEMENTARY_METRIC = "msd"
DEFAULT_METRICS = (PRIMARY_METRIC, COMPLEMENTARY_METRIC)
SUPPORTED_METRICS = (
    "msd",
    "var_dx",
    "mad_dx",
    "mean_dx",
    "median_dx",
    "ppos",
)


# =============================================================================
# CLI
# =============================================================================


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Test calendar interval-position structure and anchored lag trends "
            "using latent xstar-conditioned TT transitions."
        ),
    )
    parser.add_argument("--transitions", required=True, help="Step-5 latent transition table.")
    parser.add_argument(
        "--step11-dir",
        dest="step11_dir",
        required=True,
        help=(
            "Current Step-11 output directory from "
            "11-longitudinal_fluctuation_dynamics.py."
        ),
    )
    parser.add_argument(
        "--step13-dir",
        dest="step13_dir",
        required=True,
        help=(
            "Current Step-13 output directory from "
            "13-temporal_fluctuation_scaling.py."
        ),
    )
    parser.add_argument("--outdir", "--out-dir", dest="out_dir", required=True)
    parser.add_argument("--dataset-label", default="healthy")

    parser.add_argument(
        "--dt-values",
        default="1,2,3,4,5",
        help="Comma-separated nominal lag values retained from the transition table.",
    )
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help=(
            "Comma-separated subject-level metrics used in permutation/anchor analyses. "
            "Primary: var_dx; complementary: msd. Other supported metrics are optional sensitivities."
        ),
    )
    parser.add_argument(
        "--core-bin-ids",
        default="",
        help=(
            "Optional explicit comma-separated Step-11 bin IDs. Empty uses "
            "selected_for_analysis from Step-13 02_tested_abundance_bins."
        ),
    )
    parser.add_argument(
        "--pattern",
        default="auto",
        help=(
            "Transition pattern to retain when a pattern column exists. 'auto' prefers "
            "'all' when present, otherwise the unique available pattern."
        ),
    )

    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--n-boot-anchor", type=int, default=2000)
    parser.add_argument("--ci", type=float, default=95.0)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--min-interval-positions",
        type=int,
        default=3,
        help="Minimum calendar positions required for the primary permutation test.",
    )
    parser.add_argument(
        "--min-subjects-per-interval",
        type=int,
        default=3,
        help="Minimum contributing subjects required for a calendar interval position.",
    )
    parser.add_argument(
        "--min-anchor-dt-points",
        type=int,
        default=4,
        help="Minimum dt values required to fit an anchored signed slope.",
    )
    parser.add_argument(
        "--min-matched-anchor-subjects",
        type=int,
        default=3,
        help=(
            "Minimum complete-case subjects required for inferential bootstrap output "
            "in matched_all_dt anchored sensitivity analyses."
        ),
    )

    parser.add_argument("--streaming", action="store_true")
    parser.add_argument("--write-csv", action="store_true", default=True)
    parser.add_argument("--no-write-csv", action="store_false", dest="write_csv")
    parser.add_argument("--write-parquet", action="store_true", default=True)
    parser.add_argument("--no-write-parquet", action="store_false", dest="write_parquet")
    parser.add_argument("--compression", default="zstd")
    return parser


# =============================================================================
# Generic helpers
# =============================================================================


def parse_csv_strings(text: str) -> List[str]:
    return [value.strip() for value in str(text or "").split(",") if value.strip()]


def parse_csv_ints(text: str) -> List[int]:
    values: List[int] = []
    for raw in parse_csv_strings(text):
        value = float(raw)
        if not value.is_integer():
            raise ValueError("Expected integer value, got {!r}".format(raw))
        values.append(int(value))
    return values


def collect(lf: pl.LazyFrame, streaming: bool) -> pl.DataFrame:
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
    raise ValueError("Unsupported input format: {}".format(path))


def read_table(path: Path) -> pl.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pl.read_parquet(path)
    if suffix == ".csv":
        return pl.read_csv(path, infer_schema_length=10000)
    if suffix in {".tsv", ".txt"}:
        return pl.read_csv(path, separator="\t", infer_schema_length=10000)
    raise ValueError("Unsupported table format: {}".format(path))


def find_table(folder: Path, stem: str, required: bool = True) -> Optional[Path]:
    for suffix in (".parquet", ".pq", ".csv", ".tsv", ".txt"):
        candidate = folder / (stem + suffix)
        if candidate.exists():
            return candidate
    if required:
        raise FileNotFoundError("Could not find '{}' in {}".format(stem, folder))
    return None


def write_table(
    frame: pl.DataFrame,
    stem: Path,
    write_csv: bool,
    write_parquet: bool,
    compression: str,
) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    if write_csv:
        frame.write_csv(stem.with_suffix(".csv"))
    if write_parquet:
        frame.write_parquet(stem.with_suffix(".parquet"), compression=compression)


def load_json(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        obj = json.load(handle)
    return obj if isinstance(obj, dict) else {}


def bool_series_to_mask(series: pl.Series) -> pl.Series:
    if series.dtype == pl.Boolean:
        return series.fill_null(False)
    return (
        series.cast(pl.String, strict=False)
        .str.to_lowercase()
        .is_in(["true", "1", "yes", "y"])
        .fill_null(False)
    )


def validate_args(args: argparse.Namespace) -> None:
    args.dt_values_parsed = parse_csv_ints(args.dt_values)
    args.metrics_parsed = parse_csv_strings(args.metrics)
    args.core_bin_ids_parsed = parse_csv_ints(args.core_bin_ids) if args.core_bin_ids else []

    if not args.dt_values_parsed:
        raise ValueError("--dt-values must contain at least one lag")
    unknown_metrics = sorted(set(args.metrics_parsed) - set(SUPPORTED_METRICS))
    if unknown_metrics:
        raise ValueError("Unsupported metrics: {}".format(unknown_metrics))
    if PRIMARY_METRIC not in args.metrics_parsed:
        warnings.warn(
            "Primary metric 'var_dx' is not included in --metrics; this is a non-primary Step-14 run."
        )
    if args.n_perm <= 0:
        raise ValueError("--n-perm must be > 0")
    if args.n_boot_anchor <= 0:
        raise ValueError("--n-boot-anchor must be > 0")
    if not (0 < args.ci < 100):
        raise ValueError("--ci must be between 0 and 100")
    if args.min_interval_positions < 2:
        raise ValueError("--min-interval-positions must be >= 2")
    if args.min_subjects_per_interval < 2:
        raise ValueError("--min-subjects-per-interval must be >= 2")
    if args.min_anchor_dt_points < 3:
        raise ValueError("--min-anchor-dt-points must be >= 3")
    if args.min_matched_anchor_subjects < 2:
        raise ValueError("--min-matched-anchor-subjects must be >= 2")


# =============================================================================
# Step-11 / Step-13 design reuse
# =============================================================================


def validate_step11_metadata(step11_dir: Path, requested_dts: Sequence[int]) -> Dict[str, object]:
    metadata = load_json(step11_dir / "00_analysis_metadata.json")
    if not metadata:
        raise ValueError(
            "Step-11 00_analysis_metadata.json is required for frozen-representation "
            "validation and was not found or is empty."
        )

    representation = str(metadata.get("representation", "")).lower()
    conditioning = str(metadata.get("conditioning", "")).lower()
    mode = str(metadata.get("mode", "")).upper()
    if representation and representation != "latent":
        raise ValueError("Step-11 representation is {!r}, expected 'latent'.".format(representation))
    if conditioning and conditioning != "xstar":
        raise ValueError("Step-11 conditioning is {!r}, expected 'xstar'.".format(conditioning))
    if mode and mode != "TT":
        raise ValueError("Step-11 mode is {!r}, expected 'TT'.".format(mode))

    resolved = metadata.get("resolved_columns", {})
    if not isinstance(resolved, dict):
        raise ValueError("Step-11 metadata have no valid resolved_columns dictionary.")
    dx_source = str(resolved.get("dx", ""))
    x_source = str(resolved.get("x_condition", ""))
    if dx_source not in {"dx_latent", "__computed_x1_minus_x0"}:
        raise ValueError(
            "Step-11 resolved dx source is {!r}, expected 'dx_latent'.".format(dx_source)
        )
    if x_source != "xstar_latent":
        raise ValueError(
            "Step-11 resolved conditioning source is {!r}, expected 'xstar_latent'.".format(
                x_source
            )
        )

    analysis_layer = str(metadata.get("analysis_layer", ""))
    if analysis_layer and analysis_layer != "step11_longitudinal_fluctuation_dynamics":
        raise ValueError(
            "Step-11 analysis_layer is {!r}, expected 'step11_longitudinal_fluctuation_dynamics'.".format(
                analysis_layer
            )
        )
    if metadata.get("is_primary_configuration") is False:
        raise ValueError(
            "Step-11 metadata indicate a non-primary representation/conditioning/mode configuration."
        )

    selected_dt = metadata.get("selected_dt")
    if isinstance(selected_dt, list):
        available = sorted(int(value) for value in selected_dt)
        missing = sorted(set(requested_dts) - set(available))
        if missing:
            raise ValueError(
                "Requested dt values {} are absent from Step-11 metadata {}.".format(
                    missing, available
                )
            )
    return metadata


def validate_step13_config(step13_dir: Path) -> Dict[str, object]:
    config = load_json(step13_dir / "00_run_config.json")
    if not config:
        raise ValueError(
            "Step-13 00_run_config.json is required for frozen-representation "
            "validation and was not found or is empty."
        )

    expected = {
        "upstream_representation_validated": True,
        "upstream_required_representation": "latent",
        "upstream_required_conditioning": "xstar",
        "upstream_required_x_condition_source": "xstar_latent",
        "upstream_required_mode": "TT",
    }
    mismatches = []
    for key, wanted in expected.items():
        observed = config.get(key)
        if observed != wanted:
            mismatches.append("{}={!r} (expected {!r})".format(key, observed, wanted))
    if mismatches:
        raise ValueError(
            "Step-13 configuration is not compatible with the frozen xstar analysis:\n  - "
            + "\n  - ".join(mismatches)
        )

    layer = str(config.get("analysis_layer", ""))
    if layer and layer != "step13_temporal_fluctuation_scaling":
        raise ValueError(
            "Step-13 analysis_layer is {!r}, expected 'step13_temporal_fluctuation_scaling'.".format(
                layer
            )
        )
    if str(config.get("primary_metric", "var_dx")) != "var_dx":
        raise ValueError("Step-13 primary_metric is not 'var_dx'.")
    if config.get("pseudo_overlap_used_for_temporal_scaling") is True:
        raise ValueError(
            "Step-13 used pseudo-overlap restriction for temporal scaling; Step 14 requires the longitudinal core."
        )
    return config


def load_design_tables(
    step11_dir: Path,
    step13_dir: Path,
    explicit_core_bins: Sequence[int],
) -> Tuple[pl.DataFrame, pl.DataFrame, List[int], Path, Path]:
    bins_path = find_table(step11_dir, "04_shared_absolute_bin_definitions", required=True)
    tested_path = find_table(step13_dir, "02_tested_abundance_bins", required=True)
    assert bins_path is not None
    assert tested_path is not None

    bins = read_table(bins_path)
    tested = read_table(tested_path)

    required_bin_columns = {"bin_id", "bin_left", "bin_mid", "bin_right"}
    missing_bins = sorted(required_bin_columns - set(bins.columns))
    if missing_bins:
        raise ValueError("Step-11 bin table missing columns: {}".format(missing_bins))

    bins = bins.with_columns(
        [
            pl.col("bin_id").cast(pl.Int32, strict=False),
            pl.col("bin_left").cast(pl.Float64, strict=False),
            pl.col("bin_mid").cast(pl.Float64, strict=False),
            pl.col("bin_right").cast(pl.Float64, strict=False),
        ]
    ).drop_nulls(["bin_id", "bin_left", "bin_mid", "bin_right"]).sort("bin_id")

    if explicit_core_bins:
        core_bins = sorted(set(int(value) for value in explicit_core_bins))
    else:
        if "selected_for_analysis" not in tested.columns:
            raise ValueError(
                "Step-13 tested-bin table has no selected_for_analysis column. "
                "Provide --core-bin-ids explicitly."
            )
        mask = bool_series_to_mask(tested["selected_for_analysis"])
        core_bins = sorted(
            tested.filter(mask)
            .select(pl.col("bin_id").cast(pl.Int32, strict=False))
            .drop_nulls()
            .get_column("bin_id")
            .to_list()
        )

    if not core_bins:
        raise ValueError("No Step-13 core bins were selected.")

    available_ids = set(int(value) for value in bins["bin_id"].to_list())
    missing_ids = sorted(set(core_bins) - available_ids)
    if missing_ids:
        raise ValueError("Core bin IDs absent from Step-11 bin table: {}".format(missing_ids))

    core = bins.filter(pl.col("bin_id").is_in(core_bins)).sort("bin_id")
    return bins, core, core_bins, bins_path, tested_path


def validate_uniform_bins(bin_table: pl.DataFrame) -> Tuple[float, float, int]:
    ordered = bin_table.sort("bin_id")
    ids = [int(value) for value in ordered["bin_id"].to_list()]
    expected = list(range(min(ids), max(ids) + 1))
    if ids != expected:
        raise ValueError("Step-11 bin IDs are not contiguous; cannot reproduce assignment safely.")

    left = ordered["bin_left"].to_numpy().astype(float)
    right = ordered["bin_right"].to_numpy().astype(float)
    widths = right - left
    width = float(np.nanmedian(widths))
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Invalid Step-11 bin width")
    if not np.allclose(widths, width, rtol=1e-8, atol=1e-10):
        raise ValueError("Step-11 bins are not uniformly spaced; explicit interval join is required.")
    return float(left[0]), width, int(ids[0])


# =============================================================================
# Transition preparation
# =============================================================================


def normalize_bool_expr(column: str) -> pl.Expr:
    text = pl.col(column).cast(pl.String, strict=False).str.to_lowercase()
    return text.is_in(["true", "1", "t", "yes", "y"])


def resolve_pattern(
    lf: pl.LazyFrame,
    schema_names: Sequence[str],
    requested: str,
    dt_values: Sequence[int],
    streaming: bool,
) -> Optional[str]:
    if "pattern" not in schema_names:
        return None
    requested = str(requested or "auto").strip()
    if requested.lower() != "auto":
        return requested

    dt_expr = (
        pl.col("dt").cast(pl.Int32, strict=False)
        if "dt" in schema_names
        else pl.col("t1").cast(pl.Int32, strict=False) - pl.col("t0").cast(pl.Int32, strict=False)
    )
    patterns = collect(
        lf.with_columns(dt_expr.alias("_dt"))
        .filter(pl.col("_dt").is_in(list(dt_values)))
        .select(pl.col("pattern").cast(pl.String).unique()),
        streaming,
    ).get_column("pattern").drop_nulls().to_list()
    patterns = sorted(set(str(value) for value in patterns))
    if not patterns:
        return None
    if "all" in patterns:
        return "all"
    if len(patterns) == 1:
        return patterns[0]
    raise ValueError(
        "Multiple transition patterns are present {} and none is 'all'. "
        "Specify --pattern explicitly.".format(patterns)
    )


def prepare_filtered_transitions(
    lf: pl.LazyFrame,
    schema_names: Sequence[str],
    args: argparse.Namespace,
    full_bins: pl.DataFrame,
    core_bins: Sequence[int],
) -> Tuple[pl.DataFrame, Dict[str, object]]:
    required = ["subject", "aaSeqCDR3", "t0", "t1", "xstar_latent", "dx_latent"]
    missing = [name for name in required if name not in schema_names]
    if missing:
        raise ValueError("Transition table missing required columns: {}".format(missing))

    if "dt" in schema_names:
        dt_expr = pl.col("dt").cast(pl.Int32, strict=False)
    else:
        dt_expr = (
            pl.col("t1").cast(pl.Int32, strict=False)
            - pl.col("t0").cast(pl.Int32, strict=False)
        )

    if "obs_class" in schema_names:
        obs_class_expr = pl.col("obs_class").cast(pl.String).str.to_uppercase()
    elif "obs0" in schema_names and "obs1" in schema_names:
        obs0 = normalize_bool_expr("obs0")
        obs1 = normalize_bool_expr("obs1")
        obs_class_expr = (
            pl.when(obs0 & obs1)
            .then(pl.lit("TT"))
            .when(obs0 & ~obs1)
            .then(pl.lit("TF"))
            .when(~obs0 & obs1)
            .then(pl.lit("FT"))
            .otherwise(pl.lit("FF"))
        )
    else:
        raise ValueError("Transition input requires obs_class or both obs0 and obs1.")

    effective_dt_source = None
    for candidate in ("effective_dt", "dt_effective", "delta_t_effective"):
        if candidate in schema_names:
            effective_dt_source = candidate
            break
    effective_dt_expr = (
        pl.col(effective_dt_source).cast(pl.Float64, strict=False)
        if effective_dt_source is not None
        else dt_expr.cast(pl.Float64)
    )

    pattern_value = resolve_pattern(
        lf=lf,
        schema_names=schema_names,
        requested=args.pattern,
        dt_values=args.dt_values_parsed,
        streaming=args.streaming,
    )

    select_exprs = [
        pl.col("subject").cast(pl.String).alias("subject"),
        pl.col("aaSeqCDR3").cast(pl.String).alias("aaSeqCDR3"),
        pl.col("t0").cast(pl.Int32, strict=False).alias("t0"),
        pl.col("t1").cast(pl.Int32, strict=False).alias("t1"),
        dt_expr.alias("dt"),
        effective_dt_expr.alias("effective_dt"),
        pl.col("xstar_latent").cast(pl.Float64, strict=False).alias("xstar_latent"),
        pl.col("dx_latent").cast(pl.Float64, strict=False).alias("dx_latent"),
        obs_class_expr.alias("obs_class"),
    ]
    if "pattern" in schema_names:
        select_exprs.append(pl.col("pattern").cast(pl.String).alias("pattern"))

    selected = lf.select(select_exprs).filter(
        pl.col("dt").is_in(args.dt_values_parsed)
        & (pl.col("obs_class") == "TT")
        & pl.col("xstar_latent").is_finite()
        & pl.col("dx_latent").is_finite()
    )
    if pattern_value is not None and "pattern" in schema_names:
        selected = selected.filter(pl.col("pattern") == pattern_value)

    x_min, bin_width, first_bin_id = validate_uniform_bins(full_bins)
    n_bins = full_bins.height
    max_right = float(full_bins["bin_right"].max())

    selected = (
        selected.filter(
            (pl.col("xstar_latent") >= x_min)
            & (pl.col("xstar_latent") <= max_right)
        )
        .with_columns(
            (
                ((pl.col("xstar_latent") - x_min) / bin_width)
                .floor()
                .clip(0, n_bins - 1)
                .cast(pl.Int32)
                + first_bin_id
            ).alias("bin_id")
        )
        .filter(pl.col("bin_id").is_in(list(core_bins)))
    )

    data = collect(selected, args.streaming)
    if data.height == 0:
        raise ValueError("No transitions remain after TT/dt/core-bin filtering.")

    # Guard against accidental duplicate all+adjacent representations or repeated rows.
    duplicate_keys = ["subject", "aaSeqCDR3", "t0", "t1", "bin_id"]
    duplicates = (
        data.group_by(duplicate_keys)
        .len()
        .filter(pl.col("len") > 1)
    )
    if duplicates.height > 0:
        raise ValueError(
            "Filtered transition table contains {:,} duplicated subject-clone-interval-bin keys. "
            "Check --pattern or the Step-5 transition table.".format(duplicates.height)
        )

    resolved = {
        "representation": "latent",
        "conditioning": "xstar",
        "mode": "TT",
        "x_condition": "xstar_latent",
        "dx": "dx_latent",
        "pattern": pattern_value if pattern_value is not None else "__not_present",
        "effective_dt": effective_dt_source if effective_dt_source is not None else "__nominal_dt",
    }
    return data, resolved


# =============================================================================
# Interval summaries
# =============================================================================


def add_bin_edges(frame: pl.DataFrame, core_bins: pl.DataFrame) -> pl.DataFrame:
    return frame.join(
        core_bins.select(["bin_id", "bin_left", "bin_mid", "bin_right"]),
        on="bin_id",
        how="left",
    )


def build_interval_composition(data: pl.DataFrame) -> pl.DataFrame:
    return (
        data.group_by(["dt", "t0", "t1"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes_global"),
                pl.struct(["subject", "aaSeqCDR3"]).n_unique().alias("n_subject_clonotype_pairs"),
                pl.col("bin_id").n_unique().alias("n_core_bins_present"),
                pl.col("effective_dt").mean().alias("mean_effective_dt"),
                pl.col("effective_dt").std(ddof=1).alias("sd_effective_dt"),
                pl.col("xstar_latent").mean().alias("mean_xstar_latent"),
                pl.col("xstar_latent").std(ddof=1).alias("sd_xstar_latent"),
            ]
        )
        .with_columns(
            pl.concat_str(
                [pl.col("t0").cast(pl.String), pl.lit("->"), pl.col("t1").cast(pl.String)]
            ).alias("interval_label")
        )
        .sort(["dt", "t0", "t1"])
    )


def build_subject_interval_metrics(
    data: pl.DataFrame,
    core_bins: pl.DataFrame,
) -> pl.DataFrame:
    keys = ["subject", "dt", "t0", "t1", "bin_id"]
    base = (
        data.group_by(keys)
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes"),
                pl.col("xstar_latent").mean().alias("mean_xstar_latent"),
                pl.col("xstar_latent").median().alias("median_xstar_latent"),
                pl.col("effective_dt").mean().alias("mean_effective_dt"),
                pl.col("dx_latent").sum().alias("sum_dx"),
                (pl.col("dx_latent") ** 2).sum().alias("sum_dx2"),
                pl.col("dx_latent").mean().alias("mean_dx"),
                pl.col("dx_latent").median().alias("median_dx"),
                (pl.col("dx_latent") > 0).mean().alias("ppos"),
                (pl.col("dx_latent") ** 2).mean().alias("msd"),
                pl.col("dx_latent").var(ddof=1).alias("var_dx"),
            ]
        )
        .sort(keys)
    )

    medians = base.select(keys + ["median_dx"])
    mad = (
        data.join(medians, on=keys, how="left")
        .with_columns((pl.col("dx_latent") - pl.col("median_dx")).abs().alias("_abs_dev"))
        .group_by(keys)
        .agg(pl.col("_abs_dev").median().alias("mad_dx"))
    )

    out = base.join(mad, on=keys, how="left")
    out = add_bin_edges(out, core_bins)
    return (
        out.with_columns(
            [
                pl.concat_str(
                    [pl.col("t0").cast(pl.String), pl.lit("->"), pl.col("t1").cast(pl.String)]
                ).alias("interval_label"),
                pl.col("t0").alias("interval_position"),
            ]
        )
        .sort(["dt", "bin_id", "t0", "subject"])
    )


def build_pooled_interval_metrics(
    data: pl.DataFrame,
    core_bins: pl.DataFrame,
) -> pl.DataFrame:
    keys = ["dt", "t0", "t1", "bin_id"]
    base = (
        data.group_by(keys)
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_subjects"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_clonotypes_global"),
                pl.struct(["subject", "aaSeqCDR3"]).n_unique().alias("n_subject_clonotype_pairs"),
                pl.col("xstar_latent").mean().alias("mean_xstar_latent"),
                pl.col("xstar_latent").median().alias("median_xstar_latent"),
                pl.col("dx_latent").mean().alias("mean_dx_transition_weighted"),
                pl.col("dx_latent").median().alias("median_dx_transition_weighted"),
                (pl.col("dx_latent") > 0).mean().alias("ppos_transition_weighted"),
                (pl.col("dx_latent") ** 2).mean().alias("msd_transition_weighted"),
                pl.col("dx_latent").var(ddof=1).alias("var_dx_transition_weighted"),
            ]
        )
        .sort(keys)
    )
    medians = base.select(keys + [pl.col("median_dx_transition_weighted").alias("_median")])
    mad = (
        data.join(medians, on=keys, how="left")
        .with_columns((pl.col("dx_latent") - pl.col("_median")).abs().alias("_abs_dev"))
        .group_by(keys)
        .agg(pl.col("_abs_dev").median().alias("mad_dx_transition_weighted"))
    )
    return add_bin_edges(base.join(mad, on=keys, how="left"), core_bins)


def build_cohort_interval_profiles(
    subject_metrics: pl.DataFrame,
    pooled: pl.DataFrame,
) -> pl.DataFrame:
    keys = ["dt", "t0", "t1", "bin_id", "bin_left", "bin_mid", "bin_right"]
    equal = (
        subject_metrics.group_by(keys)
        .agg(
            [
                pl.col("subject").n_unique().alias("n_subjects_equal"),
                pl.col("n_transitions").sum().alias("n_transitions_subject_sum"),
                pl.col("mean_dx").mean().alias("mean_dx_equal_subject"),
                pl.col("median_dx").median().alias("median_dx_equal_subject"),
                pl.col("ppos").mean().alias("ppos_equal_subject"),
                pl.col("msd").mean().alias("msd_equal_subject"),
                pl.col("var_dx").mean().alias("var_dx_equal_subject"),
                pl.col("mad_dx").mean().alias("mad_dx_equal_subject"),
            ]
        )
        .sort(keys)
    )

    out = pooled.join(equal, on=keys, how="left")
    return (
        out.with_columns(
            [
                pl.concat_str(
                    [pl.col("t0").cast(pl.String), pl.lit("->"), pl.col("t1").cast(pl.String)]
                ).alias("interval_label"),
                pl.col("t0").alias("interval_position"),
            ]
        )
        .sort(["dt", "bin_id", "t0"])
    )


# =============================================================================
# Permutation tests for synchronized calendar-position structure
# =============================================================================


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if int(valid.sum()) < 3:
        return np.nan
    xv = x[valid]
    yv = y[valid]
    if np.nanstd(xv) <= 0 or np.nanstd(yv) <= 0:
        return 0.0
    return float(np.corrcoef(xv, yv)[0, 1])


def profile_statistics(positions: np.ndarray, profile: np.ndarray) -> Dict[str, float]:
    valid = np.isfinite(positions) & np.isfinite(profile)
    x = positions[valid].astype(float)
    y = profile[valid].astype(float)
    if y.size < 2:
        return {
            "rms_position_deviation": np.nan,
            "position_range": np.nan,
            "calendar_slope": np.nan,
            "alternation_score": np.nan,
            "n_centered_sign_changes": np.nan,
        }
    centered = y - float(np.mean(y))
    rms = float(np.sqrt(np.mean(centered ** 2)))
    value_range = float(np.max(y) - np.min(y))
    if np.unique(x).size >= 2:
        slope = float(np.polyfit(x, y, deg=1)[0])
    else:
        slope = np.nan

    alternation = np.nan
    sign_changes = np.nan
    if y.size >= 4:
        alt = np.where(np.arange(y.size) % 2 == 0, 1.0, -1.0)
        alt = alt - np.mean(alt)
        alternation = abs(safe_corr(centered, alt))

        nonzero = np.sign(centered)
        # Ignore exact-zero positions in the sign-change count.
        nonzero = nonzero[nonzero != 0]
        if nonzero.size >= 2:
            sign_changes = float(np.sum(nonzero[1:] != nonzero[:-1]))
        else:
            sign_changes = 0.0

    return {
        "rms_position_deviation": rms,
        "position_range": value_range,
        "calendar_slope": slope,
        "alternation_score": alternation,
        "n_centered_sign_changes": sign_changes,
    }


def empirical_p_greater_equal(null_values: np.ndarray, observed: float) -> float:
    finite = null_values[np.isfinite(null_values)]
    if finite.size == 0 or not np.isfinite(observed):
        return np.nan
    return float((1.0 + np.sum(finite >= observed)) / (finite.size + 1.0))


def quantile(values: np.ndarray, q: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan
    return float(np.quantile(finite, q))


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    valid_idx = np.flatnonzero(np.isfinite(p))
    if valid_idx.size == 0:
        return q
    pv = p[valid_idx]
    order = np.argsort(pv)
    ranked = pv[order]
    m = float(len(ranked))
    adjusted = ranked * m / np.arange(1, len(ranked) + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    q[valid_idx] = adjusted[inverse]
    return q


def group_to_subject_position_matrix(
    group: pl.DataFrame,
    metric: str,
    min_subjects_per_interval: int,
) -> Tuple[np.ndarray, np.ndarray, List[str], np.ndarray]:
    positions_all = sorted(int(value) for value in group["t0"].unique().to_list())
    position_counts = (
        group.group_by("t0")
        .agg(pl.col("subject").n_unique().alias("n_subjects"))
        .sort("t0")
    )
    valid_positions = [
        int(row["t0"])
        for row in position_counts.iter_rows(named=True)
        if int(row["n_subjects"]) >= min_subjects_per_interval
    ]
    positions = np.asarray([p for p in positions_all if p in valid_positions], dtype=float)
    if positions.size == 0:
        return positions, np.empty((0, 0), dtype=float), [], np.array([], dtype=int)

    subjects = sorted(str(value) for value in group["subject"].unique().to_list())
    subject_index = {subject: i for i, subject in enumerate(subjects)}
    position_index = {int(position): j for j, position in enumerate(positions.astype(int))}
    matrix = np.full((len(subjects), len(positions)), np.nan, dtype=float)

    for row in group.select(["subject", "t0", metric]).iter_rows(named=True):
        subject = str(row["subject"])
        position = int(row["t0"])
        value = row[metric]
        if position not in position_index or value is None:
            continue
        value_float = float(value)
        if np.isfinite(value_float):
            matrix[subject_index[subject], position_index[position]] = value_float

    counts = np.sum(np.isfinite(matrix), axis=0).astype(int)
    return positions, matrix, subjects, counts


def permute_subject_rows(matrix: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    permuted = matrix.copy()
    for row_index in range(matrix.shape[0]):
        finite_idx = np.flatnonzero(np.isfinite(matrix[row_index]))
        if finite_idx.size <= 1:
            continue
        values = matrix[row_index, finite_idx].copy()
        rng.shuffle(values)
        permuted[row_index, finite_idx] = values
    return permuted


def run_interval_position_permutations(
    subject_metrics: pl.DataFrame,
    metrics: Sequence[str],
    n_perm: int,
    seed: int,
    ci: float,
    min_positions: int,
    min_subjects_per_interval: int,
) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    rows: List[Dict[str, object]] = []
    alpha = (100.0 - ci) / 200.0

    for metric in metrics:
        if metric not in subject_metrics.columns:
            continue
        metric_frame = subject_metrics.filter(pl.col(metric).is_finite())
        for group in metric_frame.partition_by(["dt", "bin_id"], maintain_order=True):
            dt = int(group["dt"][0])
            bin_id = int(group["bin_id"][0])
            bin_left = float(group["bin_left"][0])
            bin_mid = float(group["bin_mid"][0])
            bin_right = float(group["bin_right"][0])

            positions, matrix, subjects, counts = group_to_subject_position_matrix(
                group=group,
                metric=metric,
                min_subjects_per_interval=min_subjects_per_interval,
            )
            n_positions = int(positions.size)
            base: Dict[str, object] = {
                "metric": metric,
                "dt": dt,
                "bin_id": bin_id,
                "bin_left": bin_left,
                "bin_mid": bin_mid,
                "bin_right": bin_right,
                "n_interval_positions": n_positions,
                "interval_positions": ",".join(str(int(value)) for value in positions),
                "min_subjects_across_positions": int(counts.min()) if counts.size else 0,
                "max_subjects_across_positions": int(counts.max()) if counts.size else 0,
                "n_unique_subjects": int(len(subjects)),
                "n_permutations_requested": int(n_perm),
            }

            if n_positions < min_positions:
                base.update(
                    {
                        "test_status": "insufficient_interval_positions",
                        "n_permutations_valid": 0,
                    }
                )
                rows.append(base)
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                profile = np.nanmean(matrix, axis=0)
            observed = profile_statistics(positions, profile)
            base.update(observed)

            null_rms = np.full(n_perm, np.nan, dtype=float)
            null_abs_slope = np.full(n_perm, np.nan, dtype=float)
            null_alt = np.full(n_perm, np.nan, dtype=float)

            for permutation_id in range(n_perm):
                permuted = permute_subject_rows(matrix, rng)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    perm_profile = np.nanmean(permuted, axis=0)
                stats = profile_statistics(positions, perm_profile)
                null_rms[permutation_id] = stats["rms_position_deviation"]
                null_abs_slope[permutation_id] = abs(stats["calendar_slope"])
                null_alt[permutation_id] = stats["alternation_score"]

            observed_abs_slope = abs(float(observed["calendar_slope"])) if np.isfinite(observed["calendar_slope"]) else np.nan
            p_rms = empirical_p_greater_equal(null_rms, float(observed["rms_position_deviation"]))
            p_slope = empirical_p_greater_equal(null_abs_slope, observed_abs_slope)
            p_alt = empirical_p_greater_equal(null_alt, float(observed["alternation_score"]))

            base.update(
                {
                    "test_status": "ok",
                    "n_permutations_valid": int(np.isfinite(null_rms).sum()),
                    "p_rms_position_structure": p_rms,
                    "p_abs_calendar_slope": p_slope,
                    "p_alternation": p_alt,
                    "null_rms_median": quantile(null_rms, 0.5),
                    "null_rms_q_low": quantile(null_rms, alpha),
                    "null_rms_q_high": quantile(null_rms, 1.0 - alpha),
                    "null_abs_slope_median": quantile(null_abs_slope, 0.5),
                    "null_abs_slope_q_low": quantile(null_abs_slope, alpha),
                    "null_abs_slope_q_high": quantile(null_abs_slope, 1.0 - alpha),
                    "null_alternation_median": quantile(null_alt, 0.5),
                    "null_alternation_q_low": quantile(null_alt, alpha),
                    "null_alternation_q_high": quantile(null_alt, 1.0 - alpha),
                }
            )
            rows.append(base)

    result = pl.DataFrame(rows) if rows else pl.DataFrame()
    if result.height == 0:
        return result

    # BH correction is applied across all dt x bin tests within each metric.
    result_rows = result.to_dicts()
    for metric in sorted(set(str(row.get("metric")) for row in result_rows)):
        indices = [i for i, row in enumerate(result_rows) if str(row.get("metric")) == metric]
        for p_column, q_column in (
            ("p_rms_position_structure", "q_bh_rms_position_structure"),
            ("p_abs_calendar_slope", "q_bh_abs_calendar_slope"),
            ("p_alternation", "q_bh_alternation"),
        ):
            p_values = [
                float(result_rows[i].get(p_column, np.nan))
                if result_rows[i].get(p_column) is not None
                else np.nan
                for i in indices
            ]
            q_values = benjamini_hochberg(p_values)
            for local_index, row_index in enumerate(indices):
                result_rows[row_index][q_column] = float(q_values[local_index]) if np.isfinite(q_values[local_index]) else None

    return pl.DataFrame(result_rows).sort(["metric", "dt", "bin_id"])


def summarize_permutation_tests(tests: pl.DataFrame) -> pl.DataFrame:
    if tests.height == 0:
        return pl.DataFrame()
    valid = tests.filter(pl.col("test_status") == "ok")
    if valid.height == 0:
        return pl.DataFrame()

    rows: List[Dict[str, object]] = []
    for group in valid.partition_by("metric", maintain_order=True):
        metric = str(group["metric"][0])
        row: Dict[str, object] = {
            "metric": metric,
            "n_tests": int(group.height),
            "n_dt_tested": int(group["dt"].n_unique()),
            "n_bins_tested": int(group["bin_id"].n_unique()),
        }
        for column, label in (
            ("q_bh_rms_position_structure", "rms_position_structure"),
            ("q_bh_abs_calendar_slope", "calendar_slope"),
            ("q_bh_alternation", "alternation"),
        ):
            if column in group.columns:
                values = group[column].cast(pl.Float64, strict=False)
                finite = values.drop_nulls()
                row["n_{}_q_lt_0_05".format(label)] = int((finite < 0.05).sum()) if len(finite) else 0
                row["min_{}_q".format(label)] = float(finite.min()) if len(finite) else np.nan
        rows.append(row)
    return pl.DataFrame(rows).sort("metric")


# =============================================================================
# Anchored temporal profiles and joint subject bootstrap
# =============================================================================


def metric_equal_subject_aggregate(metric: str) -> pl.Expr:
    if metric == "median_dx":
        return pl.col(metric).median()
    return pl.col(metric).mean()


def complete_case_subjects(
    subjects: Sequence[str],
    matrix: np.ndarray,
) -> Tuple[List[str], np.ndarray]:
    """Return subjects with finite values at every requested dt."""
    if matrix.ndim != 2 or matrix.shape[0] != len(subjects):
        raise ValueError("Subject/matrix shape mismatch in complete-case matching.")
    if matrix.shape[1] == 0:
        return [], np.zeros(len(subjects), dtype=bool)
    mask = np.all(np.isfinite(matrix), axis=1)
    matched = [str(subjects[i]) for i in np.flatnonzero(mask)]
    return matched, mask


def aggregate_subject_matrix(metric: str, matrix: np.ndarray) -> np.ndarray:
    """Aggregate a subject x dt matrix using the Step-14 equal-subject diagnostic rule."""
    if matrix.size == 0:
        return np.array([], dtype=float)
    if metric == "median_dx":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmedian(matrix, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(matrix, axis=0)


def build_anchor_subject_matching_diagnostic(
    subject_metrics: pl.DataFrame,
    dt_values: Sequence[int],
) -> pl.DataFrame:
    """
    Structural anchor availability independent of metric and abundance bin.

    A subject is marked complete when at least one subject-interval summary row is
    present at every requested dt for the relevant common-start/common-end anchor.
    This table is descriptive; strict inferential matching is additionally applied
    within each metric x abundance bin, where finite metric values are required at
    every requested dt.
    """
    min_t0 = int(subject_metrics["t0"].min())
    max_t1 = int(subject_metrics["t1"].max())
    requested = [int(value) for value in dt_values]
    rows: List[Dict[str, object]] = []

    for anchor_type, anchor_value, condition_column in (
        ("common_start", min_t0, "t0"),
        ("common_end", max_t1, "t1"),
    ):
        anchored = subject_metrics.filter(pl.col(condition_column) == anchor_value)
        subjects = sorted(str(value) for value in subject_metrics["subject"].unique().to_list())
        for subject in subjects:
            g = anchored.filter(pl.col("subject").cast(pl.String) == subject)
            available = sorted(
                set(int(value) for value in g["dt"].unique().to_list()) & set(requested)
            ) if g.height else []
            missing = [value for value in requested if value not in available]
            rows.append(
                {
                    "anchor_type": anchor_type,
                    "anchor_value": int(anchor_value),
                    "subject": subject,
                    "n_requested_dt": int(len(requested)),
                    "n_available_dt": int(len(available)),
                    "available_dt_values": ",".join(str(value) for value in available),
                    "missing_dt_values": ",".join(str(value) for value in missing),
                    "complete_all_requested_dt": bool(len(missing) == 0),
                }
            )

    out = pl.DataFrame(rows) if rows else pl.DataFrame()
    if out.height:
        out = out.sort(["anchor_type", "subject"])
    return out


def build_anchored_profiles(
    subject_metrics: pl.DataFrame,
    metrics: Sequence[str],
    dt_values: Sequence[int],
) -> Tuple[pl.DataFrame, int, int]:
    """
    Build both original and strict matched-subject anchored profiles.

    subject_matching='available' reproduces the original availability-based Step-14 profile, where
    each dt uses all finite subject values available at that lag.

    subject_matching='matched_all_dt' retains, separately within each anchor x
    metric x abundance bin, only subjects with a finite metric at every requested
    dt. The contributing subject set is therefore identical across all lags.
    """
    min_t0 = int(subject_metrics["t0"].min())
    max_t1 = int(subject_metrics["t1"].max())
    requested_dt = [int(value) for value in dt_values]
    rows: List[Dict[str, object]] = []

    for anchor_type, anchor_value, condition_column in (
        ("common_start", min_t0, "t0"),
        ("common_end", max_t1, "t1"),
    ):
        anchored = subject_metrics.filter(pl.col(condition_column) == anchor_value)
        bin_ids = sorted(int(value) for value in anchored["bin_id"].unique().to_list())

        for metric in metrics:
            if metric not in anchored.columns:
                continue

            # Original available-subject profile.
            summary = (
                anchored.filter(pl.col(metric).is_finite())
                .group_by(["dt", "bin_id", "bin_left", "bin_mid", "bin_right"])
                .agg(
                    [
                        pl.col("subject").n_unique().alias("n_subjects"),
                        pl.col("n_transitions").sum().alias("n_transitions"),
                        metric_equal_subject_aggregate(metric).alias("equal_subject_value"),
                    ]
                )
                .sort(["bin_id", "dt"])
            )
            for item in summary.iter_rows(named=True):
                row = dict(item)
                row.update(
                    {
                        "anchor_type": anchor_type,
                        "anchor_value": int(anchor_value),
                        "metric": metric,
                        "subject_matching": "available",
                        "n_subjects_union": None,
                        "matched_subject_ids": "",
                    }
                )
                rows.append(row)

            # Strict complete-case matched-subject profile.
            for bin_id in bin_ids:
                frame = anchored.filter(
                    (pl.col("bin_id") == bin_id) & pl.col(metric).is_finite()
                )
                if frame.height == 0:
                    continue
                subjects, dt_array, matrix = build_anchor_matrix(frame, metric, requested_dt)
                matched_subjects, complete_mask = complete_case_subjects(subjects, matrix)
                if not matched_subjects:
                    continue

                matched_matrix = matrix[complete_mask, :]
                matched_profile = aggregate_subject_matrix(metric, matched_matrix)
                bin_row = frame.select(["bin_left", "bin_mid", "bin_right"]).row(0, named=True)

                matched_set = set(matched_subjects)
                matched_frame = frame.filter(
                    pl.col("subject").cast(pl.String).is_in(sorted(matched_set))
                )
                transition_counts = {
                    int(row["dt"]): int(row["n_transitions"])
                    for row in (
                        matched_frame.group_by("dt")
                        .agg(pl.col("n_transitions").sum().alias("n_transitions"))
                        .iter_rows(named=True)
                    )
                }

                for j, dt in enumerate(requested_dt):
                    value = float(matched_profile[j]) if j < matched_profile.size else np.nan
                    rows.append(
                        {
                            "dt": int(dt),
                            "bin_id": int(bin_id),
                            "bin_left": float(bin_row["bin_left"]),
                            "bin_mid": float(bin_row["bin_mid"]),
                            "bin_right": float(bin_row["bin_right"]),
                            "n_subjects": int(len(matched_subjects)),
                            "n_transitions": int(transition_counts.get(int(dt), 0)),
                            "equal_subject_value": value,
                            "anchor_type": anchor_type,
                            "anchor_value": int(anchor_value),
                            "metric": metric,
                            "subject_matching": "matched_all_dt",
                            "n_subjects_union": int(len(subjects)),
                            "matched_subject_ids": ",".join(matched_subjects),
                        }
                    )

    out = pl.DataFrame(rows) if rows else pl.DataFrame()
    if out.height:
        out = out.sort(["subject_matching", "anchor_type", "metric", "bin_id", "dt"])
    return out, min_t0, max_t1


def signed_linear_fit(dt: np.ndarray, values: np.ndarray, min_points: int) -> Tuple[float, float, int]:
    valid = np.isfinite(dt) & np.isfinite(values)
    x = dt[valid].astype(float)
    y = values[valid].astype(float)
    n = int(y.size)
    if n < min_points or np.unique(x).size < 2:
        return np.nan, np.nan, n
    design = np.column_stack([np.ones(n, dtype=float), x - 1.0])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return float(coefficients[0]), float(coefficients[1]), n


def build_anchor_matrix(
    frame: pl.DataFrame,
    metric: str,
    dt_values: Sequence[int],
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    subjects = sorted(str(value) for value in frame["subject"].unique().to_list())
    subject_index = {subject: i for i, subject in enumerate(subjects)}
    dt_index = {int(dt): j for j, dt in enumerate(dt_values)}
    matrix = np.full((len(subjects), len(dt_values)), np.nan, dtype=float)

    for row in frame.select(["subject", "dt", metric]).iter_rows(named=True):
        subject = str(row["subject"])
        dt = int(row["dt"])
        value = row[metric]
        if dt not in dt_index or value is None:
            continue
        value_float = float(value)
        if np.isfinite(value_float):
            matrix[subject_index[subject], dt_index[dt]] = value_float
    return subjects, np.asarray(dt_values, dtype=float), matrix


def weighted_subject_mean(matrix: np.ndarray, multiplicity: np.ndarray) -> np.ndarray:
    valid = np.isfinite(matrix)
    weighted_values = np.nan_to_num(matrix, nan=0.0) * multiplicity[:, None]
    numerator = np.sum(weighted_values, axis=0)
    denominator = np.sum(valid * multiplicity[:, None], axis=0)
    return np.divide(
        numerator,
        denominator,
        out=np.full(matrix.shape[1], np.nan, dtype=float),
        where=denominator > 0,
    )


def bootstrap_anchor_matrix(
    matrix: np.ndarray,
    dt_array: np.ndarray,
    metric: str,
    n_boot: int,
    rng: np.random.Generator,
    min_points: int,
) -> Tuple[np.ndarray, np.ndarray, float, float, int]:
    multiplicity = np.ones(matrix.shape[0], dtype=float)
    if metric == "median_dx":
        observed_profile = aggregate_subject_matrix(metric, matrix)
    else:
        observed_profile = weighted_subject_mean(matrix, multiplicity)
    observed_intercept, observed_slope, n_observed_dt = signed_linear_fit(
        dt_array, observed_profile, min_points
    )

    boot_slopes = np.full(n_boot, np.nan, dtype=float)
    if matrix.shape[0] >= 2:
        for bootstrap_id in range(n_boot):
            sampled = rng.integers(0, matrix.shape[0], size=matrix.shape[0])
            if metric == "median_dx":
                sampled_matrix = matrix[sampled, :]
                profile = aggregate_subject_matrix(metric, sampled_matrix)
            else:
                mult = np.bincount(sampled, minlength=matrix.shape[0]).astype(float)
                profile = weighted_subject_mean(matrix, mult)
            _, slope, _ = signed_linear_fit(dt_array, profile, min_points)
            boot_slopes[bootstrap_id] = slope
    return observed_profile, boot_slopes, observed_intercept, observed_slope, n_observed_dt


def anchored_slope_bootstrap(
    subject_metrics: pl.DataFrame,
    metrics: Sequence[str],
    dt_values: Sequence[int],
    n_boot: int,
    seed: int,
    ci: float,
    min_points: int,
    min_matched_subjects: int,
) -> pl.DataFrame:
    """Run original and strict matched-subject anchored signed-slope bootstraps."""
    rng = np.random.default_rng(seed + 7919)
    alpha = (100.0 - ci) / 200.0
    min_t0 = int(subject_metrics["t0"].min())
    max_t1 = int(subject_metrics["t1"].max())
    rows: List[Dict[str, object]] = []

    for anchor_type, anchor_value, condition_column in (
        ("common_start", min_t0, "t0"),
        ("common_end", max_t1, "t1"),
    ):
        anchored = subject_metrics.filter(pl.col(condition_column) == anchor_value)
        bin_ids = sorted(int(value) for value in anchored["bin_id"].unique().to_list())
        for metric in metrics:
            if metric not in anchored.columns:
                continue
            for bin_id in bin_ids:
                frame = anchored.filter((pl.col("bin_id") == bin_id) & pl.col(metric).is_finite())
                if frame.height == 0:
                    continue
                subjects, dt_array, matrix = build_anchor_matrix(frame, metric, dt_values)
                bin_row = frame.select(["bin_left", "bin_mid", "bin_right"]).row(0, named=True)

                for subject_matching in ("available", "matched_all_dt"):
                    if subject_matching == "matched_all_dt":
                        matched_subjects, complete_mask = complete_case_subjects(subjects, matrix)
                        analysis_subjects = matched_subjects
                        analysis_matrix = matrix[complete_mask, :]
                    else:
                        analysis_subjects = list(subjects)
                        analysis_matrix = matrix

                    counts_by_dt = np.sum(np.isfinite(analysis_matrix), axis=0).astype(int) if analysis_matrix.size else np.zeros(len(dt_values), dtype=int)
                    n_analysis_subjects = int(len(analysis_subjects))
                    test_status = "ok"
                    if subject_matching == "matched_all_dt" and n_analysis_subjects < min_matched_subjects:
                        test_status = "insufficient_matched_subjects"
                    elif n_analysis_subjects < 2:
                        test_status = "insufficient_subjects"

                    observed_profile, boot_slopes, observed_intercept, observed_slope, n_observed_dt = bootstrap_anchor_matrix(
                        matrix=analysis_matrix,
                        dt_array=dt_array,
                        metric=metric,
                        n_boot=(n_boot if test_status == "ok" else 0),
                        rng=rng,
                        min_points=min_points,
                    )

                    finite = boot_slopes[np.isfinite(boot_slopes)]
                    if finite.size:
                        slope_median = float(np.quantile(finite, 0.5))
                        slope_low = float(np.quantile(finite, alpha))
                        slope_high = float(np.quantile(finite, 1.0 - alpha))
                        p_negative = float(np.mean(finite < 0))
                        p_positive = float(np.mean(finite > 0))
                    else:
                        slope_median = slope_low = slope_high = p_negative = p_positive = np.nan

                    rows.append(
                        {
                            "anchor_type": anchor_type,
                            "anchor_value": int(anchor_value),
                            "subject_matching": subject_matching,
                            "metric": metric,
                            "bin_id": int(bin_id),
                            "bin_left": float(bin_row["bin_left"]),
                            "bin_mid": float(bin_row["bin_mid"]),
                            "bin_right": float(bin_row["bin_right"]),
                            "test_status": test_status,
                            "n_subjects_union": int(len(subjects)),
                            "n_subjects_analysis": n_analysis_subjects,
                            "min_subjects_across_dt": int(counts_by_dt.min()) if counts_by_dt.size else 0,
                            "n_observed_dt": int(n_observed_dt),
                            "observed_intercept_at_dt1": observed_intercept,
                            "observed_signed_slope": observed_slope,
                            "n_bootstrap_requested": int(n_boot),
                            "n_valid_bootstrap": int(finite.size),
                            "signed_slope_median": slope_median,
                            "signed_slope_q_low": slope_low,
                            "signed_slope_q_high": slope_high,
                            "signed_slope_negative_fraction": p_negative,
                            "signed_slope_positive_fraction": p_positive,
                            "signed_slope_ci_below_zero": bool(np.isfinite(slope_high) and slope_high < 0),
                            "signed_slope_ci_above_zero": bool(np.isfinite(slope_low) and slope_low > 0),
                            "dt_values": ",".join(str(int(value)) for value in dt_values),
                            "subject_counts_by_dt": ",".join(str(int(value)) for value in counts_by_dt),
                            "analysis_subject_ids": ",".join(analysis_subjects),
                            "observed_profile_values": ",".join(
                                "{:.12g}".format(value) if np.isfinite(value) else "nan"
                                for value in observed_profile
                            ),
                        }
                    )

    out = pl.DataFrame(rows) if rows else pl.DataFrame()
    if out.height:
        out = out.sort(["subject_matching", "anchor_type", "metric", "bin_id"])
    return out


# =============================================================================
# Compact overall summary
# =============================================================================


def build_analysis_summary(
    permutation_summary: pl.DataFrame,
    anchor_slopes: pl.DataFrame,
) -> pl.DataFrame:
    metrics = set()
    if permutation_summary.height:
        metrics.update(str(value) for value in permutation_summary["metric"].to_list())
    if anchor_slopes.height:
        metrics.update(str(value) for value in anchor_slopes["metric"].to_list())

    rows: List[Dict[str, object]] = []
    for metric in sorted(metrics):
        row: Dict[str, object] = {"metric": metric}
        if permutation_summary.height:
            p = permutation_summary.filter(pl.col("metric") == metric)
            if p.height:
                for column in p.columns:
                    if column == "metric":
                        continue
                    row["permutation_{}".format(column)] = p[column][0]

        if anchor_slopes.height:
            a_metric = anchor_slopes.filter(pl.col("metric") == metric)
            matching_modes = (
                ["available", "matched_all_dt"]
                if "subject_matching" in a_metric.columns
                else ["available"]
            )
            for subject_matching in matching_modes:
                a = (
                    a_metric.filter(pl.col("subject_matching") == subject_matching)
                    if "subject_matching" in a_metric.columns
                    else a_metric
                )
                prefix = "" if subject_matching == "available" else "matched_"
                for anchor_type in ("common_start", "common_end"):
                    g = a.filter(pl.col("anchor_type") == anchor_type)
                    if "test_status" in g.columns:
                        g_ok = g.filter(pl.col("test_status") == "ok")
                    else:
                        g_ok = g
                    if g.height:
                        row["{}{}_n_bins_total".format(prefix, anchor_type)] = int(g.height)
                        row["{}{}_n_bins_valid".format(prefix, anchor_type)] = int(g_ok.height)
                        if "n_subjects_analysis" in g.columns:
                            finite_subjects = g["n_subjects_analysis"].cast(pl.Int64, strict=False).drop_nulls()
                            row["{}{}_median_n_subjects".format(prefix, anchor_type)] = (
                                float(finite_subjects.median()) if len(finite_subjects) else np.nan
                            )
                            row["{}{}_min_n_subjects".format(prefix, anchor_type)] = (
                                int(finite_subjects.min()) if len(finite_subjects) else 0
                            )
                    if g_ok.height:
                        row["{}{}_n_negative_ci".format(prefix, anchor_type)] = int(
                            g_ok["signed_slope_ci_below_zero"].cast(pl.Boolean).sum()
                        )
                        row["{}{}_n_positive_ci".format(prefix, anchor_type)] = int(
                            g_ok["signed_slope_ci_above_zero"].cast(pl.Boolean).sum()
                        )
                        finite = g_ok["signed_slope_median"].cast(pl.Float64, strict=False).drop_nulls()
                        row["{}{}_median_signed_slope".format(prefix, anchor_type)] = (
                            float(finite.median()) if len(finite) else np.nan
                        )

                start = a.filter(pl.col("anchor_type") == "common_start")
                end = a.filter(pl.col("anchor_type") == "common_end")
                if "test_status" in start.columns:
                    start = start.filter(pl.col("test_status") == "ok")
                if "test_status" in end.columns:
                    end = end.filter(pl.col("test_status") == "ok")
                start = start.select(
                    ["bin_id", pl.col("signed_slope_ci_below_zero").alias("start_negative")]
                )
                end = end.select(
                    ["bin_id", pl.col("signed_slope_ci_below_zero").alias("end_negative")]
                )
                both = start.join(end, on="bin_id", how="inner")
                if both.height:
                    row["{}n_bins_negative_in_both_anchors".format(prefix)] = int(
                        (pl.Series(both["start_negative"]) & pl.Series(both["end_negative"])).sum()
                    )

        rows.append(row)
    return pl.DataFrame(rows) if rows else pl.DataFrame()


# =============================================================================
# Main
# =============================================================================


def main() -> int:
    args = build_argparser().parse_args()
    validate_args(args)

    transitions_path = Path(args.transitions).expanduser()
    step11_dir = Path(args.step11_dir).expanduser()
    step13_dir = Path(args.step13_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not transitions_path.exists():
        raise FileNotFoundError("Transition table not found: {}".format(transitions_path))
    if not step11_dir.exists():
        raise FileNotFoundError("Step-11 directory not found: {}".format(step11_dir))
    if not step13_dir.exists():
        raise FileNotFoundError("Step-13 directory not found: {}".format(step13_dir))

    step11_metadata = validate_step11_metadata(step11_dir, args.dt_values_parsed)
    step13_config = validate_step13_config(step13_dir)
    full_bins, core_bins, core_bin_ids, bins_path, tested_path = load_design_tables(
        step11_dir=step11_dir,
        step13_dir=step13_dir,
        explicit_core_bins=args.core_bin_ids_parsed,
    )

    if not args.core_bin_ids_parsed and isinstance(step13_config.get("tested_bins"), list):
        configured_bins = sorted(int(value) for value in step13_config["tested_bins"])
        if configured_bins != core_bin_ids:
            raise ValueError(
                "Step-13 tested_bins in 00_run_config.json do not match "
                "selected_for_analysis in 02_tested_abundance_bins: {} vs {}.".format(
                    configured_bins, core_bin_ids
                )
            )

    lf = scan_table(transitions_path)
    schema = lf.collect_schema()
    schema_names = schema.names()

    print("[INFO] ClonoDynamics Step 14 | {}".format(SCRIPT_VERSION))
    print("[READ] {}".format(transitions_path))
    print("[DESIGN] core bins={}".format(core_bin_ids))
    print(
        "[DESIGN] core x range [{:.6f}, {:.6f}]".format(
            float(core_bins["bin_left"].min()),
            float(core_bins["bin_right"].max()),
        )
    )

    data, resolved = prepare_filtered_transitions(
        lf=lf,
        schema_names=schema_names,
        args=args,
        full_bins=full_bins,
        core_bins=core_bin_ids,
    )
    print("[FILTER] transitions={:,}".format(data.height))

    manifest = pl.DataFrame(
        [
            {
                "dataset_label": args.dataset_label,
                "transitions_path": str(transitions_path),
                "step11_dir": str(step11_dir),
                "step13_dir": str(step13_dir),
                "step11_bin_table": str(bins_path),
                "step13_tested_bin_table": str(tested_path),
                "n_transition_rows_filtered": int(data.height),
                "n_subjects": int(data["subject"].n_unique()),
                "n_core_bins": int(len(core_bin_ids)),
                "core_bin_ids": ",".join(str(value) for value in core_bin_ids),
                "core_x_min": float(core_bins["bin_left"].min()),
                "core_x_max": float(core_bins["bin_right"].max()),
                "dt_values": ",".join(str(value) for value in args.dt_values_parsed),
                "resolved_pattern": str(resolved["pattern"]),
                "effective_dt_source": str(resolved["effective_dt"]),
                "conditioning": str(resolved["conditioning"]),
                "x_condition_source": str(resolved["x_condition"]),
                "dx_source": str(resolved["dx"]),
                "step11_metadata_validated": True,
                "step13_config_validated": True,
            }
        ]
    )
    write_table(manifest, out_dir / "01_input_manifest", args.write_csv, args.write_parquet, args.compression)
    write_table(core_bins, out_dir / "02_core_bin_definitions", args.write_csv, args.write_parquet, args.compression)

    composition = build_interval_composition(data)
    write_table(composition, out_dir / "03_interval_composition_by_pair", args.write_csv, args.write_parquet, args.compression)

    print("[SUMMARY] Building subject x interval metrics")
    subject_metrics = build_subject_interval_metrics(data, core_bins)
    write_table(subject_metrics, out_dir / "04_subject_interval_metrics", args.write_csv, args.write_parquet, args.compression)

    pooled = build_pooled_interval_metrics(data, core_bins)
    cohort_profiles = build_cohort_interval_profiles(subject_metrics, pooled)
    write_table(cohort_profiles, out_dir / "05_cohort_interval_profiles", args.write_csv, args.write_parquet, args.compression)

    print("[PERM] {} permutations".format(args.n_perm))
    permutation_tests = run_interval_position_permutations(
        subject_metrics=subject_metrics,
        metrics=args.metrics_parsed,
        n_perm=args.n_perm,
        seed=args.seed,
        ci=args.ci,
        min_positions=args.min_interval_positions,
        min_subjects_per_interval=args.min_subjects_per_interval,
    )
    write_table(permutation_tests, out_dir / "06_interval_position_permutation_tests", args.write_csv, args.write_parquet, args.compression)

    permutation_summary = summarize_permutation_tests(permutation_tests)
    write_table(permutation_summary, out_dir / "07_interval_position_test_summary", args.write_csv, args.write_parquet, args.compression)

    anchor_matching_diagnostic = build_anchor_subject_matching_diagnostic(
        subject_metrics=subject_metrics,
        dt_values=args.dt_values_parsed,
    )
    write_table(
        anchor_matching_diagnostic,
        out_dir / "11_anchor_subject_matching_diagnostic",
        args.write_csv,
        args.write_parquet,
        args.compression,
    )

    anchored_profiles, min_t0, max_t1 = build_anchored_profiles(
        subject_metrics=subject_metrics,
        metrics=args.metrics_parsed,
        dt_values=args.dt_values_parsed,
    )
    write_table(anchored_profiles, out_dir / "08_anchored_temporal_profiles", args.write_csv, args.write_parquet, args.compression)

    print(
        "[BOOT] Anchored signed slopes: {} subject bootstraps; available + matched_all_dt".format(
            args.n_boot_anchor
        )
    )
    anchor_slopes = anchored_slope_bootstrap(
        subject_metrics=subject_metrics,
        metrics=args.metrics_parsed,
        dt_values=args.dt_values_parsed,
        n_boot=args.n_boot_anchor,
        seed=args.seed,
        ci=args.ci,
        min_points=args.min_anchor_dt_points,
        min_matched_subjects=args.min_matched_anchor_subjects,
    )
    write_table(anchor_slopes, out_dir / "09_anchored_signed_slope_bootstrap", args.write_csv, args.write_parquet, args.compression)

    analysis_summary = build_analysis_summary(permutation_summary, anchor_slopes)
    write_table(analysis_summary, out_dir / "10_analysis_summary", args.write_csv, args.write_parquet, args.compression)

    config = vars(args).copy()
    for key in ("dt_values_parsed", "metrics_parsed", "core_bin_ids_parsed"):
        config[key] = list(config[key])
    config.update(
        {
            "step": 14,
            "script": "14-interval_position_structure.py",
            "script_version": SCRIPT_VERSION,
            "analysis_layer": "step14_interval_position_structure",
            "upstream_steps": [5, 11, 13],
            "step12_pseudo_baseline_used": False,
            "generated": datetime.now().isoformat(timespec="seconds"),
            "resolved_columns": resolved,
            "frozen_representation": "latent dx_latent conditioned on xstar_latent",
            "primary_metric": PRIMARY_METRIC,
            "complementary_metric": COMPLEMENTARY_METRIC,
            "analysis_role": "calendar-position and anchored composition-control sensitivity downstream of Step 13; does not replace Step 13 temporal-model inference",
            "pseudo_overlap_used": False,
            "representation": "latent",
            "conditioning": "xstar",
            "x_condition_source": "xstar_latent",
            "dx_source": "dx_latent",
            "mode": "TT",
            "core_bin_ids_resolved": core_bin_ids,
            "core_x_min": float(core_bins["bin_left"].min()),
            "core_x_max": float(core_bins["bin_right"].max()),
            "earliest_timepoint": int(min_t0),
            "latest_timepoint": int(max_t1),
            "permutation_unit": "independent interval-label permutation within subject x dt x abundance bin",
            "interval_position_estimator": "equal-subject cohort profile",
            "permutation_primary_statistic": "RMS deviation of equal-subject cohort interval profile",
            "permutation_multiple_testing": "Benjamini-Hochberg within metric across dt x abundance-bin tests",
            "alternation_is_exploratory": True,
            "sinusoidal_model_fitted": False,
            "anchor_bootstrap_unit": "joint subject cluster across dt",
            "anchor_estimator": "equal-subject aggregation of subject-level interval metrics",
            "anchor_subject_modes": ["available", "matched_all_dt"],
            "matched_anchor_definition": (
                "strict complete-case subjects within anchor x metric x abundance bin; "
                "finite values required at every requested dt"
            ),
            "min_matched_anchor_subjects": int(args.min_matched_anchor_subjects),
            "step11_metadata_loaded": bool(step11_metadata),
            "step11_metadata_validated": True,
            "step13_config_loaded": bool(step13_config),
            "step13_config_validated": True,
        }
    )
    with open(out_dir / "00_run_config.json", "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)

    with open(out_dir / "DEBUG_REPORT.txt", "w", encoding="utf-8") as handle:
        handle.write("STEP 14 — INTERVAL-POSITION / ANCHORED SENSITIVITY\n")
        handle.write("======================================\n\n")
        handle.write("Dataset: {}\n".format(args.dataset_label))
        handle.write("Transitions: {}\n".format(transitions_path))
        handle.write("Representation: latent / xstar_latent conditioning / TT\n")
        handle.write("dt values: {}\n".format(args.dt_values_parsed))
        handle.write("Core bins: {}\n".format(core_bin_ids))
        handle.write(
            "Core x range: [{:.8g}, {:.8g}]\n".format(
                float(core_bins["bin_left"].min()),
                float(core_bins["bin_right"].max()),
            )
        )
        handle.write("Filtered transitions: {:,}\n".format(data.height))
        handle.write("Subjects: {}\n".format(data["subject"].n_unique()))
        handle.write("Metrics: {}\n".format(args.metrics_parsed))
        handle.write("Primary metric: var_dx; complementary metric: msd\n")
        handle.write("Step-13 pseudo-overlap restriction reused: no\n")
        handle.write("Permutations: {}\n".format(args.n_perm))
        handle.write("Anchor bootstraps: {}\n".format(args.n_boot_anchor))
        handle.write("Anchor subject modes: available, matched_all_dt\n")
        handle.write(
            "Minimum matched subjects for inference: {}\n".format(
                args.min_matched_anchor_subjects
            )
        )
        handle.write("Earliest t0 / latest t1: {} / {}\n\n".format(min_t0, max_t1))
        if anchor_matching_diagnostic.height:
            handle.write("Structural matched-subject anchor availability:\n")
            for anchor_type in ("common_start", "common_end"):
                g = anchor_matching_diagnostic.filter(pl.col("anchor_type") == anchor_type)
                n_complete = int(g["complete_all_requested_dt"].cast(pl.Boolean).sum()) if g.height else 0
                complete_ids = (
                    g.filter(pl.col("complete_all_requested_dt") == True)["subject"].to_list()
                    if g.height else []
                )
                handle.write(
                    "  {}: {} complete subjects [{}]\n".format(
                        anchor_type, n_complete, ",".join(str(value) for value in complete_ids)
                    )
                )
            handle.write("\n")
        handle.write("Permutation interpretation:\n")
        handle.write(
            "  Tests synchronized calendar-position structure after independent within-subject interval-label shuffling.\n"
        )
        handle.write(
            "  The alternation score is exploratory and is not a sinusoidal/periodic model test.\n\n"
        )
        if permutation_summary.height:
            handle.write("Permutation summary:\n")
            handle.write(str(permutation_summary))
            handle.write("\n\n")
        if analysis_summary.height:
            handle.write("Overall summary:\n")
            handle.write(str(analysis_summary))
            handle.write("\n")

    print("[DONE] {}".format(out_dir))
    if permutation_summary.height:
        print(permutation_summary)
    if analysis_summary.height:
        print(analysis_summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise