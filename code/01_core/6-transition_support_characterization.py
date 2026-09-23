#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
6-transition_support_characterization.py
========================================

Descriptive characterization and production-readiness audit of the generic
longitudinal transition table used by the replicate-resolved ClonoDynamics
framework.

Overview
--------
This script implements Step 6 of the production ClonoDynamics pipeline.

Canonical upstream scripts are:

    Step 4
        4-multirepresentation_trajectory_assembly.py

    Step 5
        5-longitudinal_transition_assembly.py

Step 6 does NOT estimate longitudinal dynamics.  It characterizes the
transition universe assembled in Step 5 and verifies that the fields required
for the downstream estimands are structurally available on their intended
support.

The production sequence is therefore:

    multirepresentation clonotype states
                |
                v
    generic endpoint pairing / transitions
                |
                v
    transition + support characterization       <- this script
                |
                v
    estimand-specific downstream analyses

The central distinction from earlier latent-centric versions is that the
transition table contains three complementary information layers:

    1. model-based consensus / conditioning
           x0_latent, x1_latent, xmid_latent, optional xstar_latent;

    2. replicate-resolved observed measurements
           replicate-specific counts, depths, frequencies and log-frequencies;

    3. observation annotations
           realized read positivity, p_value, reference-threshold observability,
           operational TT/TF/FT/FF classes, and model-based detectability.

Step 6 reports the structure and completeness of these layers without choosing
or fitting a dynamical model.

Primary downstream estimands supported by this table
-----------------------------------------------------
Forward drift is defined downstream from observed technical replicates by the
replicate-decoupled AB/BA construction:

    AB:
        conditioning = x_observed_rep1(t0)
        displacement = dx_observed_rep2

    BA:
        conditioning = x_observed_rep2(t0)
        displacement = dx_observed_rep1.

Step 5 stores the corresponding eligibility flags:

    forward_ab_eligible
    forward_ba_eligible.

Replicate-consistent fluctuations are defined downstream from the
cross-replicate covariance of observed displacements:

    Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt)

on:

    common4 = both technical replicates read-positive at both endpoints.

Step 6 therefore treats AB support, BA support and common4 support as core
transition-dataset properties.

Operational observation classes
-------------------------------
Step 5 preserves the continuous empirical observation statistic at both
endpoints:

    p_value_t0
    p_value_t1

and the reference-threshold annotation:

    observable_reference_t0
    observable_reference_t1
    obs_class_reference.

At the reference alpha,

    observable_reference = 1[p_value < reference_alpha].

The four classes:

    TT, TF, FT, FF

are OPERATIONAL endpoint-observation classes.  They are not interpreted as
biological persistence, contraction, expansion or absence.

The continuous p_value fields are retained so later sensitivity analyses can
reconstruct the operational observation domain at alternative alpha values
without rebuilding the transition table.

What this script does
---------------------
The script produces five kinds of descriptive/audit information:

    1. transition-universe size and temporal coverage;

    2. coverage by subject, subject-specific interval, nominal interval and dt;

    3. reference-threshold TT/TF/FT/FF composition;

    4. estimand support:
           forward_ab_eligible
           forward_ba_eligible
           common4;

    5. production-readiness checks:
           core-field completeness;
           support-conditional completeness of observed coordinates;
           consistency of stored positivity/support/observability flags;
           optional model-field completeness.

It intentionally does NOT:

    - refit Step-2 latent states;
    - redefine Step-3 diagnostics;
    - apply abundance or posterior-quality filters;
    - remove TT/TF/FT/FF transitions;
    - treat TT as a quality class;
    - estimate forward drift;
    - estimate cross-replicate covariance;
    - estimate temporal scaling;
    - perform operational-threshold sensitivity analyses;
    - count pseudo configurations as additional biological subjects.

Input
-----
The required input is the Step-5 production transition table:

    longitudinal_transitions.parquet

specified with:

    --transitions

The primary empirical analysis ordinarily has no extra grouping columns.
Pseudo/null/experimental datasets may preserve additional grouping variables,
for example:

    --extra-group-cols configuration_id

These variables are included in analysis-unit keys but are never counted as
biological subjects.

Core required columns
---------------------
Structural identifiers:

    subject
    aaSeqCDR3
    t0
    t1
    dt

Model-based conditioning layer:

    x0_latent
    x1_latent
    xmid_latent

Replicate-resolved observed layer:

    count_rep1_t0
    count_rep1_t1
    count_rep2_t0
    count_rep2_t1

    depth_rep1_t0
    depth_rep1_t1
    depth_rep2_t0
    depth_rep2_t1

    f_observed_rep1_t0
    f_observed_rep1_t1
    f_observed_rep2_t0
    f_observed_rep2_t1

    x_observed_rep1_t0
    x_observed_rep1_t1
    x_observed_rep2_t0
    x_observed_rep2_t1

    dx_observed_rep1
    dx_observed_rep2

Realized positivity and estimand support:

    read_positive_rep1_t0
    read_positive_rep1_t1
    read_positive_rep2_t0
    read_positive_rep2_t1

    forward_ab_eligible
    forward_ba_eligible
    common4

Operational observation layer:

    p_value_t0
    p_value_t1
    reference_alpha
    observable_reference_t0
    observable_reference_t1
    obs_class_reference

The script stops if a core field is absent.  This is intentional: Step 6 is a
production contract audit, not a legacy-field fallback layer.

Support-conditional completeness
--------------------------------
Observed log-frequency is undefined for a zero-count replicate because the
production Step 4 does not introduce a pseudocount.  Consequently, requiring
x_observed_rep* to be non-null on every transition would be incorrect.

Instead Step 6 tests the fields on the support where they are mathematically
required.

Examples:

    x_observed_rep1_t0
        must be finite where read_positive_rep1_t0 is True;

    dx_observed_rep1
        must be finite where replicate 1 is positive at both endpoints;

    AB readiness
        requires finite x_observed_rep1_t0 and dx_observed_rep2 on
        forward_ab_eligible rows;

    BA readiness
        requires finite x_observed_rep2_t0 and dx_observed_rep1 on
        forward_ba_eligible rows;

    common4 fluctuation readiness
        requires both observed displacements and xmid_latent to be finite on
        common4 rows.

This distinguishes legitimate observation-boundary missingness from actual
field incompleteness.

Structural consistency audit
----------------------------
Step 6 independently verifies that stored flags agree with their definitions:

    dt == t1 - t0;

    read_positive_repr == (count_repr > 0);

    forward_ab_eligible ==
        rep1 positive at t0 AND rep2 positive at t0 and t1;

    forward_ba_eligible ==
        rep2 positive at t0 AND rep1 positive at t0 and t1;

    common4 ==
        both replicates positive at both endpoints;

    observable_reference_tk ==
        (p_value_tk < reference_alpha), k in {0,1};

    obs_class_reference ==
        concatenation of the two reference endpoint labels.

These checks audit construction consistency only.  They do not validate any
biological model.

Optional model-based fields
---------------------------
The script reports presence/completeness, but does not require, quantities such
as:

    dx_latent
    xstar_latent
    x0_latent_sd
    x1_latent_sd
    dx_latent_sd_approx

    p_detect_state_t0 / p_detect_state_t1
    p_dropout_state_t0 / p_dropout_state_t1

    full-posterior transition summaries
        dx_post_mean
        dx_post_median
        dx_post_q025
        dx_post_q975
        dx_post_sd
        p_dx_gt0
        p_dx_lt0
        posterior_success.

Full-posterior propagation is therefore an optional model-based uncertainty
layer rather than a requirement for construction of the primary
replicate-resolved estimands.

Outputs
-------
The standard CSV outputs are:

    01_transition_dataset_summary.csv
        one-row transition-universe summary;

    02_subject_interval_coverage.csv
        coverage and support for each subject x t0 x t1 x dt block;

    03_subject_dt_coverage.csv
        coverage and support by subject x dt;

    04_coverage_summary_by_dt.csv
        dataset-wide coverage by dt;

    05_interval_pair_coverage.csv
        nominal t0 x t1 x dt coverage across subjects;

    06_reference_observation_class_composition_by_dt.csv
        reference-threshold TT/TF/FT/FF composition by dt;

    07_reference_observation_class_composition_overall.csv
        overall reference-threshold class composition;

    08_estimand_support_by_dt.csv
        AB, BA and common4 support by dt;

    09_estimand_support_overall.csv
        overall AB, BA and common4 support;

    10_core_field_completeness.csv
        unconditional completeness of fields expected on every row;

    11_estimand_readiness_audit.csv
        support-conditional availability of observed coordinates required by
        the primary downstream estimands;

    12_structural_consistency_audit.csv
        consistency of stored positivity, eligibility, dt and reference
        observability definitions;

    13_optional_model_field_completeness.csv
        presence/completeness of optional model-based fields;

    14_dataset_overview.csv
        compact long-format manuscript/supplementary summary.

A Markdown report is also written:

    transition_support_summary.md

and a machine-readable manifest:

    transition_support_manifest.json.

When --write-parquet is supplied, each tabular output is also written as
Parquet.

Optional descriptive figures
----------------------------
With --figures, Step 6 writes separate descriptive panels only:

    01_transition_coverage_by_dt
        transition counts and contributing subject-interval units by dt;

    02_estimand_support_by_dt
        fractions eligible for AB, BA and common4 by dt;

    03_subject_interval_coverage_heatmap
        log10 transition counts across subject x nominal interval blocks.

These are coverage/support figures, not dynamical estimates.

Typical usage
-------------
Primary longitudinal production analysis:

    python3 code/dynamics/6-transition_support_characterization.py \
        --transitions \
        results/5-longitudinal_transition_assembly/p_05/longitudinal_transitions.parquet \
        --outdir \
        results/6-transition_support_characterization/p_05 \
        --streaming \
        --write-parquet

Pseudo/null analysis preserving configuration identity:

    python3 code/dynamics/6-transition_support_characterization.py \
        --transitions pseudo_longitudinal_transitions.parquet \
        --outdir results/6-transition_support_characterization/pseudo \
        --extra-group-cols configuration_id \
        --streaming

Interpretive boundary
---------------------
Step 6 answers:

    "What transition universe was assembled, what support is available for
     each downstream estimand, and is the table structurally ready for those
     analyses?"

It does not answer:

    "What are the longitudinal drift or fluctuation dynamics?"

Those quantities are estimated only in later estimand-specific steps.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import polars as pl


SCRIPT_VERSION = "v4-replicate-resolved-transition-support-2026-09-09"
OBS_CLASSES = ("TT", "TF", "FT", "FF")

CORE_REQUIRED_COLUMNS = [
    "subject",
    "aaSeqCDR3",
    "t0",
    "t1",
    "dt",
    "x0_latent",
    "x1_latent",
    "xmid_latent",
    "count_rep1_t0",
    "count_rep1_t1",
    "count_rep2_t0",
    "count_rep2_t1",
    "depth_rep1_t0",
    "depth_rep1_t1",
    "depth_rep2_t0",
    "depth_rep2_t1",
    "f_observed_rep1_t0",
    "f_observed_rep1_t1",
    "f_observed_rep2_t0",
    "f_observed_rep2_t1",
    "x_observed_rep1_t0",
    "x_observed_rep1_t1",
    "x_observed_rep2_t0",
    "x_observed_rep2_t1",
    "dx_observed_rep1",
    "dx_observed_rep2",
    "read_positive_rep1_t0",
    "read_positive_rep1_t1",
    "read_positive_rep2_t0",
    "read_positive_rep2_t1",
    "forward_ab_eligible",
    "forward_ba_eligible",
    "common4",
    "p_value_t0",
    "p_value_t1",
    "reference_alpha",
    "observable_reference_t0",
    "observable_reference_t1",
    "obs_class_reference",
]

# Fields expected to be defined on every transition row.  Observed log-frequency
# fields are intentionally excluded because they are null by design at zero-count
# replicate states and are audited conditionally in the estimand-readiness table.
UNCONDITIONAL_CORE_FIELDS = [
    "subject",
    "aaSeqCDR3",
    "t0",
    "t1",
    "dt",
    "x0_latent",
    "x1_latent",
    "xmid_latent",
    "count_rep1_t0",
    "count_rep1_t1",
    "count_rep2_t0",
    "count_rep2_t1",
    "depth_rep1_t0",
    "depth_rep1_t1",
    "depth_rep2_t0",
    "depth_rep2_t1",
    "f_observed_rep1_t0",
    "f_observed_rep1_t1",
    "f_observed_rep2_t0",
    "f_observed_rep2_t1",
    "read_positive_rep1_t0",
    "read_positive_rep1_t1",
    "read_positive_rep2_t0",
    "read_positive_rep2_t1",
    "forward_ab_eligible",
    "forward_ba_eligible",
    "common4",
    "p_value_t0",
    "p_value_t1",
    "reference_alpha",
    "observable_reference_t0",
    "observable_reference_t1",
    "obs_class_reference",
]

OPTIONAL_MODEL_FIELDS = [
    "f0_latent",
    "f1_latent",
    "dx_latent",
    "x0_condition_latent",
    "xstar_latent",
    "xstar_source",
    "x0_latent_sd",
    "x1_latent_sd",
    "dx_latent_sd_approx",
    "f0_latent_q025",
    "f0_latent_q975",
    "f1_latent_q025",
    "f1_latent_q975",
    "x0_latent_q025",
    "x0_latent_q975",
    "x1_latent_q025",
    "x1_latent_q975",
    "p_detect_state_t0",
    "p_detect_state_t1",
    "p_dropout_state_t0",
    "p_dropout_state_t1",
    "dx_post_mean",
    "dx_post_median",
    "dx_post_q025",
    "dx_post_q975",
    "dx_post_sd",
    "p_dx_gt0",
    "p_dx_lt0",
    "posterior_success",
    "posterior_message",
    "x_positive_current_t0",
    "x_positive_current_t1",
    "dx_positive_current",
    "replicate_posterior_bhattacharyya_t0",
    "replicate_posterior_bhattacharyya_t1",
]

TRUE_STRINGS = {"true", "t", "1", "yes", "y"}


def split_csv(text: Optional[str]) -> List[str]:
    if text is None:
        return []
    return [part.strip() for part in str(text).split(",") if part.strip()]


def collect(lf: pl.LazyFrame, streaming: bool) -> pl.DataFrame:
    if streaming:
        try:
            return lf.collect(engine="streaming")
        except TypeError:
            return lf.collect(streaming=True)
    return lf.collect()


def require_columns(names: Sequence[str], required: Sequence[str], label: str) -> None:
    missing = [column for column in required if column not in names]
    if missing:
        raise ValueError(
            f"{label} is missing required production columns: {missing}. "
            "Re-run the current Step 4 and Step 5 production scripts rather "
            "than using legacy fallbacks."
        )


def bool_expr(column: str) -> pl.Expr:
    """Normalize a bool-like field to Boolean without changing row membership."""
    return (
        pl.col(column)
        .cast(pl.String, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
        .is_in(sorted(TRUE_STRINGS))
        .fill_null(False)
    )


def finite_expr(column: str) -> pl.Expr:
    """True for finite numeric values after permissive Float64 conversion."""
    x = pl.col(column).cast(pl.Float64, strict=False)
    return x.is_not_null() & x.is_finite()


def nonnull_expr(column: str) -> pl.Expr:
    return pl.col(column).is_not_null()


def safe_fraction(num: int, den: int) -> Optional[float]:
    if den <= 0:
        return None
    return float(num / den)


def write_table(
    df: pl.DataFrame,
    outdir: Path,
    stem: str,
    write_parquet: bool,
    compression: str,
) -> None:
    df.write_csv(outdir / f"{stem}.csv")
    if write_parquet:
        df.write_parquet(outdir / f"{stem}.parquet", compression=compression)


def add_overview_row(
    rows: List[Dict[str, Any]],
    metric: str,
    value: Any,
    unit: str = "",
    description: str = "",
) -> None:
    rows.append(
        {
            "metric": metric,
            "value": value,
            "unit": unit,
            "description": description,
        }
    )


def assert_single_reference_alpha(
    base: pl.LazyFrame,
    streaming: bool,
) -> float:
    values = collect(
        base.select(
            pl.col("reference_alpha")
            .cast(pl.Float64, strict=False)
            .drop_nulls()
            .unique()
            .sort()
        ),
        streaming,
    ).get_column("reference_alpha").to_list()

    values = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if len(values) != 1:
        raise ValueError(
            "Step-6 production input must contain exactly one finite "
            f"reference_alpha value; found {values}."
        )
    return float(values[0])


def build_core_completeness(
    base: pl.LazyFrame,
    schema_names: Sequence[str],
    total_rows: int,
    streaming: bool,
) -> pl.DataFrame:
    numeric_like = {
        "subject", "t0", "t1", "dt",
        "x0_latent", "x1_latent", "xmid_latent",
        "count_rep1_t0", "count_rep1_t1", "count_rep2_t0", "count_rep2_t1",
        "depth_rep1_t0", "depth_rep1_t1", "depth_rep2_t0", "depth_rep2_t1",
        "f_observed_rep1_t0", "f_observed_rep1_t1",
        "f_observed_rep2_t0", "f_observed_rep2_t1",
        "p_value_t0", "p_value_t1", "reference_alpha",
    }

    expressions: List[pl.Expr] = []
    for column in UNCONDITIONAL_CORE_FIELDS:
        valid = finite_expr(column) if column in numeric_like else nonnull_expr(column)
        expressions.extend(
            [
                valid.sum().alias(f"{column}__n_valid"),
                valid.mean().alias(f"{column}__fraction_valid"),
            ]
        )

    result = collect(base.select(expressions), streaming).row(0, named=True)
    rows: List[Dict[str, Any]] = []
    for column in UNCONDITIONAL_CORE_FIELDS:
        n_valid = int(result[f"{column}__n_valid"] or 0)
        rows.append(
            {
                "field": column,
                "present": column in set(schema_names),
                "n_expected": int(total_rows),
                "n_valid": n_valid,
                "n_missing_or_nonfinite": int(total_rows - n_valid),
                "fraction_valid": result[f"{column}__fraction_valid"],
                "expectation": "defined on every transition row",
            }
        )
    return pl.DataFrame(rows).sort("field")


def _conditional_readiness_row(
    base: pl.LazyFrame,
    streaming: bool,
    name: str,
    support: pl.Expr,
    required_valid: pl.Expr,
    description: str,
) -> Dict[str, Any]:
    result = collect(
        base.select(
            [
                support.sum().alias("n_support"),
                (support & required_valid).sum().alias("n_ready"),
            ]
        ),
        streaming,
    ).row(0, named=True)

    n_support = int(result["n_support"] or 0)
    n_ready = int(result["n_ready"] or 0)
    return {
        "audit": name,
        "n_support": n_support,
        "n_ready": n_ready,
        "n_not_ready": int(n_support - n_ready),
        "fraction_ready_within_support": safe_fraction(n_ready, n_support),
        "definition": description,
    }


def build_estimand_readiness(
    base: pl.LazyFrame,
    streaming: bool,
) -> pl.DataFrame:
    rp1_t0 = bool_expr("read_positive_rep1_t0")
    rp1_t1 = bool_expr("read_positive_rep1_t1")
    rp2_t0 = bool_expr("read_positive_rep2_t0")
    rp2_t1 = bool_expr("read_positive_rep2_t1")

    ab = bool_expr("forward_ab_eligible")
    ba = bool_expr("forward_ba_eligible")
    common4 = bool_expr("common4")

    rows = [
        _conditional_readiness_row(
            base,
            streaming,
            "rep1_x_t0_on_positive_support",
            rp1_t0,
            finite_expr("x_observed_rep1_t0"),
            "x_observed_rep1_t0 finite whenever replicate 1 is read-positive at t0",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "rep1_x_t1_on_positive_support",
            rp1_t1,
            finite_expr("x_observed_rep1_t1"),
            "x_observed_rep1_t1 finite whenever replicate 1 is read-positive at t1",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "rep2_x_t0_on_positive_support",
            rp2_t0,
            finite_expr("x_observed_rep2_t0"),
            "x_observed_rep2_t0 finite whenever replicate 2 is read-positive at t0",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "rep2_x_t1_on_positive_support",
            rp2_t1,
            finite_expr("x_observed_rep2_t1"),
            "x_observed_rep2_t1 finite whenever replicate 2 is read-positive at t1",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "rep1_displacement_on_two_endpoint_support",
            rp1_t0 & rp1_t1,
            finite_expr("dx_observed_rep1"),
            "dx_observed_rep1 finite when replicate 1 is positive at both endpoints",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "rep2_displacement_on_two_endpoint_support",
            rp2_t0 & rp2_t1,
            finite_expr("dx_observed_rep2"),
            "dx_observed_rep2 finite when replicate 2 is positive at both endpoints",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "forward_AB_estimand_ready",
            ab,
            finite_expr("x_observed_rep1_t0") & finite_expr("dx_observed_rep2"),
            "AB: finite rep1 initial conditioning and rep2 displacement on forward_ab_eligible",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "forward_BA_estimand_ready",
            ba,
            finite_expr("x_observed_rep2_t0") & finite_expr("dx_observed_rep1"),
            "BA: finite rep2 initial conditioning and rep1 displacement on forward_ba_eligible",
        ),
        _conditional_readiness_row(
            base,
            streaming,
            "common4_fluctuation_estimand_ready",
            common4,
            finite_expr("dx_observed_rep1")
            & finite_expr("dx_observed_rep2")
            & finite_expr("xmid_latent"),
            "common4: both observed displacements and latent midpoint conditioning finite",
        ),
    ]
    return pl.DataFrame(rows)


def _consistency_row(
    base: pl.LazyFrame,
    streaming: bool,
    name: str,
    comparison: pl.Expr,
    description: str,
) -> Dict[str, Any]:
    result = collect(
        base.select(
            [
                pl.len().alias("n_checked"),
                comparison.fill_null(False).sum().alias("n_match"),
            ]
        ),
        streaming,
    ).row(0, named=True)
    n_checked = int(result["n_checked"] or 0)
    n_match = int(result["n_match"] or 0)
    return {
        "audit": name,
        "n_checked": n_checked,
        "n_match": n_match,
        "n_mismatch": int(n_checked - n_match),
        "fraction_match": safe_fraction(n_match, n_checked),
        "definition": description,
    }


def build_structural_consistency(
    base: pl.LazyFrame,
    streaming: bool,
) -> pl.DataFrame:
    c1_0 = pl.col("count_rep1_t0").cast(pl.Float64, strict=False)
    c1_1 = pl.col("count_rep1_t1").cast(pl.Float64, strict=False)
    c2_0 = pl.col("count_rep2_t0").cast(pl.Float64, strict=False)
    c2_1 = pl.col("count_rep2_t1").cast(pl.Float64, strict=False)

    rp1_0 = bool_expr("read_positive_rep1_t0")
    rp1_1 = bool_expr("read_positive_rep1_t1")
    rp2_0 = bool_expr("read_positive_rep2_t0")
    rp2_1 = bool_expr("read_positive_rep2_t1")

    obs0 = bool_expr("observable_reference_t0")
    obs1 = bool_expr("observable_reference_t1")
    alpha = pl.col("reference_alpha").cast(pl.Float64, strict=False)
    p0 = pl.col("p_value_t0").cast(pl.Float64, strict=False)
    p1 = pl.col("p_value_t1").cast(pl.Float64, strict=False)

    expected_obs_class = pl.concat_str(
        [
            pl.when(obs0).then(pl.lit("T")).otherwise(pl.lit("F")),
            pl.when(obs1).then(pl.lit("T")).otherwise(pl.lit("F")),
        ]
    )

    rows = [
        _consistency_row(
            base,
            streaming,
            "dt_equals_t1_minus_t0",
            pl.col("dt").cast(pl.Int64, strict=False)
            == (
                pl.col("t1").cast(pl.Int64, strict=False)
                - pl.col("t0").cast(pl.Int64, strict=False)
            ),
            "stored dt equals t1 - t0",
        ),
        _consistency_row(
            base,
            streaming,
            "rep1_t0_read_positive_matches_count",
            rp1_0 == (c1_0 > 0),
            "read_positive_rep1_t0 equals count_rep1_t0 > 0",
        ),
        _consistency_row(
            base,
            streaming,
            "rep1_t1_read_positive_matches_count",
            rp1_1 == (c1_1 > 0),
            "read_positive_rep1_t1 equals count_rep1_t1 > 0",
        ),
        _consistency_row(
            base,
            streaming,
            "rep2_t0_read_positive_matches_count",
            rp2_0 == (c2_0 > 0),
            "read_positive_rep2_t0 equals count_rep2_t0 > 0",
        ),
        _consistency_row(
            base,
            streaming,
            "rep2_t1_read_positive_matches_count",
            rp2_1 == (c2_1 > 0),
            "read_positive_rep2_t1 equals count_rep2_t1 > 0",
        ),
        _consistency_row(
            base,
            streaming,
            "forward_ab_flag_definition",
            bool_expr("forward_ab_eligible") == (rp1_0 & rp2_0 & rp2_1),
            "forward_ab_eligible equals rep1 positive at t0 AND rep2 positive at t0 and t1",
        ),
        _consistency_row(
            base,
            streaming,
            "forward_ba_flag_definition",
            bool_expr("forward_ba_eligible") == (rp2_0 & rp1_0 & rp1_1),
            "forward_ba_eligible equals rep2 positive at t0 AND rep1 positive at t0 and t1",
        ),
        _consistency_row(
            base,
            streaming,
            "common4_flag_definition",
            bool_expr("common4") == (rp1_0 & rp1_1 & rp2_0 & rp2_1),
            "common4 equals both replicates positive at both endpoints",
        ),
        _consistency_row(
            base,
            streaming,
            "reference_observable_t0_definition",
            obs0 == (p0 < alpha),
            "observable_reference_t0 equals p_value_t0 < reference_alpha",
        ),
        _consistency_row(
            base,
            streaming,
            "reference_observable_t1_definition",
            obs1 == (p1 < alpha),
            "observable_reference_t1 equals p_value_t1 < reference_alpha",
        ),
        _consistency_row(
            base,
            streaming,
            "reference_obs_class_definition",
            pl.col("obs_class_reference").cast(pl.String, strict=False)
            == expected_obs_class,
            "obs_class_reference concatenates reference endpoint T/F labels",
        ),
        _consistency_row(
            base,
            streaming,
            "reference_obs_class_allowed_values",
            pl.col("obs_class_reference").cast(pl.String, strict=False).is_in(list(OBS_CLASSES)),
            "obs_class_reference is one of TT/TF/FT/FF",
        ),
        _consistency_row(
            base,
            streaming,
            "p_value_t0_in_unit_interval",
            finite_expr("p_value_t0") & (p0 >= 0.0) & (p0 <= 1.0),
            "p_value_t0 is finite and within [0,1]",
        ),
        _consistency_row(
            base,
            streaming,
            "p_value_t1_in_unit_interval",
            finite_expr("p_value_t1") & (p1 >= 0.0) & (p1 <= 1.0),
            "p_value_t1 is finite and within [0,1]",
        ),
        _consistency_row(
            base,
            streaming,
            "reference_alpha_in_unit_interval",
            finite_expr("reference_alpha") & (alpha > 0.0) & (alpha < 1.0),
            "reference_alpha is finite and strictly between 0 and 1",
        ),
        _consistency_row(
            base,
            streaming,
            "replicate_depths_positive",
            (pl.col("depth_rep1_t0").cast(pl.Float64, strict=False) > 0)
            & (pl.col("depth_rep1_t1").cast(pl.Float64, strict=False) > 0)
            & (pl.col("depth_rep2_t0").cast(pl.Float64, strict=False) > 0)
            & (pl.col("depth_rep2_t1").cast(pl.Float64, strict=False) > 0),
            "all replicate sequencing depths are positive",
        ),
    ]
    return pl.DataFrame(rows)


def build_optional_model_completeness(
    base: pl.LazyFrame,
    schema_names: Sequence[str],
    total_rows: int,
    streaming: bool,
) -> pl.DataFrame:
    names = set(schema_names)
    present_fields = [field for field in OPTIONAL_MODEL_FIELDS if field in names]

    values: Dict[str, Any] = {}
    if present_fields:
        expressions: List[pl.Expr] = []
        for field in present_fields:
            # posterior_message and xstar_source are text; all others are numeric/bool.
            if field in {"posterior_message", "xstar_source"}:
                valid = nonnull_expr(field)
            elif field == "posterior_success":
                valid = nonnull_expr(field)
            else:
                valid = finite_expr(field)
            expressions.extend(
                [
                    valid.sum().alias(f"{field}__n_valid"),
                    valid.mean().alias(f"{field}__fraction_valid"),
                ]
            )
        values = collect(base.select(expressions), streaming).row(0, named=True)

    rows: List[Dict[str, Any]] = []
    for field in OPTIONAL_MODEL_FIELDS:
        present = field in names
        n_valid = int(values.get(f"{field}__n_valid", 0) or 0) if present else 0
        fraction = values.get(f"{field}__fraction_valid") if present else None
        rows.append(
            {
                "field": field,
                "present": bool(present),
                "n_total": int(total_rows),
                "n_valid": n_valid,
                "fraction_valid": fraction,
                "role": "optional model-based / sensitivity / uncertainty field",
            }
        )
    return pl.DataFrame(rows).sort("field")


def make_figures(
    dt_coverage: pl.DataFrame,
    subject_interval_coverage: pl.DataFrame,
    support_by_dt: pl.DataFrame,
    figure_dir: Path,
    dpi: int,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    # 1. Coverage by dt.
    if dt_coverage.height > 0:
        x = dt_coverage["dt"].to_list()
        transitions = dt_coverage["n_transitions"].to_list()
        intervals = dt_coverage["n_group_subject_interval_pairs"].to_list()

        fig, ax1 = plt.subplots(figsize=(4.8, 3.2))
        line1 = ax1.plot(x, transitions, marker="o", linewidth=1.1, label="Transitions")[0]
        ax1.set_xlabel("Temporal lag, Δt")
        ax1.set_ylabel("Transitions")
        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax1.spines["top"].set_visible(False)

        ax2 = ax1.twinx()
        line2 = ax2.plot(
            x,
            intervals,
            marker="s",
            linewidth=1.0,
            linestyle="--",
            label="Subject-interval units",
        )[0]
        ax2.set_ylabel("Subject-interval units")
        ax2.spines["top"].set_visible(False)
        ax2.legend(
            [line1, line2],
            ["Transitions", "Subject-interval units"],
            frameon=False,
            loc="best",
        )
        fig.tight_layout()
        fig.savefig(figure_dir / "01_transition_coverage_by_dt.pdf")
        fig.savefig(figure_dir / "01_transition_coverage_by_dt.png", dpi=dpi)
        plt.close(fig)

    # 2. Estimand support fractions by dt.
    if support_by_dt.height > 0:
        x = support_by_dt["dt"].to_list()
        fig, ax = plt.subplots(figsize=(4.8, 3.2))
        ax.plot(
            x,
            support_by_dt["fraction_forward_ab_eligible"].to_list(),
            marker="o",
            linewidth=1.0,
            label="AB eligible",
        )
        ax.plot(
            x,
            support_by_dt["fraction_forward_ba_eligible"].to_list(),
            marker="s",
            linewidth=1.0,
            label="BA eligible",
        )
        ax.plot(
            x,
            support_by_dt["fraction_common4"].to_list(),
            marker="^",
            linewidth=1.0,
            label="common4",
        )
        ax.set_xlabel("Temporal lag, Δt")
        ax.set_ylabel("Fraction of transition universe")
        ax.set_ylim(0, 1)
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(figure_dir / "02_estimand_support_by_dt.pdf")
        fig.savefig(figure_dir / "02_estimand_support_by_dt.png", dpi=dpi)
        plt.close(fig)

    # 3. Subject x nominal interval coverage heatmap.
    if subject_interval_coverage.height > 0:
        intervals = (
            subject_interval_coverage.select(["t0", "t1", "dt"])
            .unique()
            .sort(["dt", "t0", "t1"])
        )
        subjects = subject_interval_coverage.select("subject").unique().sort("subject")

        interval_keys = [
            (row["t0"], row["t1"], row["dt"])
            for row in intervals.iter_rows(named=True)
        ]
        subject_values = subjects["subject"].to_list()
        interval_index = {key: idx for idx, key in enumerate(interval_keys)}
        subject_index = {value: idx for idx, value in enumerate(subject_values)}

        matrix = np.full((len(subject_values), len(interval_keys)), np.nan, dtype=float)
        # When extra grouping variables are used, multiple rows can exist for the
        # same biological subject/nominal interval. Sum them for this descriptive
        # biological-subject heatmap; extra-group-aware tables remain available.
        heatmap_rows = (
            subject_interval_coverage.group_by(["subject", "t0", "t1", "dt"])
            .agg(pl.col("n_transitions").sum().alias("n_transitions"))
        )
        for row in heatmap_rows.iter_rows(named=True):
            i = subject_index[row["subject"]]
            j = interval_index[(row["t0"], row["t1"], row["dt"])]
            matrix[i, j] = float(row["n_transitions"])

        display_matrix = np.log10(matrix + 1.0)
        width = max(5.5, 0.42 * len(interval_keys) + 1.8)
        height = max(3.4, 0.30 * len(subject_values) + 1.2)
        fig, ax = plt.subplots(figsize=(width, height))
        image = ax.imshow(display_matrix, aspect="auto", interpolation="nearest")
        labels = [f"{t0}→{t1}\n(Δt={dt})" for t0, t1, dt in interval_keys]
        ax.set_xticks(range(len(interval_keys)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks(range(len(subject_values)))
        ax.set_yticklabels([str(value) for value in subject_values])
        ax.set_xlabel("Nominal temporal interval")
        ax.set_ylabel("Biological subject")
        cbar = fig.colorbar(image, ax=ax)
        cbar.set_label("log10(transitions + 1)")
        fig.tight_layout()
        fig.savefig(figure_dir / "03_subject_interval_coverage_heatmap.pdf")
        fig.savefig(figure_dir / "03_subject_interval_coverage_heatmap.png", dpi=dpi)
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Characterize the Step-5 replicate-resolved transition universe and "
            "audit support/readiness for downstream estimands."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--transitions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--extra-group-cols",
        default="",
        help=(
            "Comma-separated grouping variables preserved upstream, e.g. "
            "configuration_id. These extend analysis-unit keys but never "
            "biological subject N."
        ),
    )
    parser.add_argument("--streaming", action="store_true")
    parser.add_argument("--write-parquet", action="store_true")
    parser.add_argument("--compression", default="zstd")
    parser.add_argument("--figures", action="store_true")
    parser.add_argument("--figure-dpi", type=int, default=300)
    parser.add_argument("--report-md", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    transitions_path = args.transitions.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    extra_group_cols = split_csv(args.extra_group_cols)

    if not transitions_path.exists():
        raise FileNotFoundError(f"Transition parquet not found: {transitions_path}")
    if transitions_path.suffix.lower() not in {".parquet", ".pq"}:
        raise ValueError("Production Step 6 expects a Parquet transition table.")

    print(f"[STEP6] {SCRIPT_VERSION}")
    print(f"[STEP6] input={transitions_path}")
    print(f"[STEP6] extra_group_cols={extra_group_cols or '<none>'}")

    lf = pl.scan_parquet(str(transitions_path))
    schema_obj = lf.collect_schema()
    schema_names = list(schema_obj.names())
    schema_set = set(schema_names)

    require_columns(schema_names, CORE_REQUIRED_COLUMNS, "Step-5 transition table")
    require_columns(schema_names, extra_group_cols, "Step-5 transition table")

    reference_alpha = assert_single_reference_alpha(lf, args.streaming)
    print(f"[STEP6] reference_alpha={reference_alpha:g}")

    # Normalize only the bool-like columns used for descriptive grouping/audits.
    base = lf.with_columns(
        [
            bool_expr("read_positive_rep1_t0").alias("_rp1_t0"),
            bool_expr("read_positive_rep1_t1").alias("_rp1_t1"),
            bool_expr("read_positive_rep2_t0").alias("_rp2_t0"),
            bool_expr("read_positive_rep2_t1").alias("_rp2_t1"),
            bool_expr("forward_ab_eligible").alias("_ab"),
            bool_expr("forward_ba_eligible").alias("_ba"),
            bool_expr("common4").alias("_common4"),
            bool_expr("observable_reference_t0").alias("_obs_ref_t0"),
            bool_expr("observable_reference_t1").alias("_obs_ref_t1"),
        ]
    )

    total_rows = int(collect(base.select(pl.len().alias("n")), args.streaming).item())
    if total_rows <= 0:
        raise ValueError("Transition table contains zero rows.")

    group_subject_key = extra_group_cols + ["subject"]
    group_subject_clone_key = extra_group_cols + ["subject", "aaSeqCDR3"]
    group_subject_interval_key = extra_group_cols + ["subject", "t0", "t1", "dt"]
    nominal_interval_key = ["t0", "t1", "dt"]

    # ------------------------------------------------------------------
    # 1. Overall transition universe
    # ------------------------------------------------------------------
    overall = collect(
        base.select(
            [
                pl.len().alias("n_transition_rows"),
                pl.col("subject").n_unique().alias("n_biological_subjects"),
                pl.struct(group_subject_key).n_unique().alias("n_group_subject_units"),
                pl.col("aaSeqCDR3").n_unique().alias("n_unique_cdr3_unstratified"),
                pl.struct(group_subject_clone_key)
                .n_unique()
                .alias("n_group_subject_clonotype_pairs"),
                pl.struct(group_subject_interval_key)
                .n_unique()
                .alias("n_group_subject_interval_pairs"),
                pl.struct(nominal_interval_key)
                .n_unique()
                .alias("n_distinct_nominal_interval_pairs"),
                pl.col("dt").min().alias("min_dt"),
                pl.col("dt").max().alias("max_dt"),
                pl.col("_ab").cast(pl.Float64).mean().alias("fraction_forward_ab_eligible"),
                pl.col("_ba").cast(pl.Float64).mean().alias("fraction_forward_ba_eligible"),
                pl.col("_common4").cast(pl.Float64).mean().alias("fraction_common4"),
            ]
        ),
        args.streaming,
    ).with_columns(pl.lit(reference_alpha).alias("reference_alpha"))

    write_table(
        overall,
        outdir,
        "01_transition_dataset_summary",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 2. Coverage by group/subject/actual interval
    # ------------------------------------------------------------------
    subject_interval_coverage = collect(
        base.group_by(group_subject_interval_key)
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("aaSeqCDR3").n_unique().alias("n_clonotypes"),
                pl.col("_ab").sum().alias("n_forward_ab_eligible"),
                pl.col("_ba").sum().alias("n_forward_ba_eligible"),
                pl.col("_common4").sum().alias("n_common4"),
                pl.col("_ab").cast(pl.Float64).mean().alias("fraction_forward_ab_eligible"),
                pl.col("_ba").cast(pl.Float64).mean().alias("fraction_forward_ba_eligible"),
                pl.col("_common4").cast(pl.Float64).mean().alias("fraction_common4"),
                *[
                    (pl.col("obs_class_reference") == cls)
                    .sum()
                    .alias(f"n_{cls}")
                    for cls in OBS_CLASSES
                ],
            ]
        )
        .with_columns(
            [
                (pl.col(f"n_{cls}") / pl.col("n_transitions")).alias(
                    f"fraction_{cls}"
                )
                for cls in OBS_CLASSES
            ]
        )
        .sort(group_subject_key + ["dt", "t0", "t1"]),
        args.streaming,
    )
    write_table(
        subject_interval_coverage,
        outdir,
        "02_subject_interval_coverage",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 3. Coverage by group/subject/dt
    # ------------------------------------------------------------------
    subject_dt_key = extra_group_cols + ["subject", "dt"]
    subject_dt_coverage = collect(
        base.group_by(subject_dt_key)
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("aaSeqCDR3").n_unique().alias("n_clonotypes"),
                pl.struct(["t0", "t1"]).n_unique().alias("n_interval_pairs"),
                pl.col("_ab").sum().alias("n_forward_ab_eligible"),
                pl.col("_ba").sum().alias("n_forward_ba_eligible"),
                pl.col("_common4").sum().alias("n_common4"),
                pl.col("_ab").cast(pl.Float64).mean().alias("fraction_forward_ab_eligible"),
                pl.col("_ba").cast(pl.Float64).mean().alias("fraction_forward_ba_eligible"),
                pl.col("_common4").cast(pl.Float64).mean().alias("fraction_common4"),
            ]
        )
        .sort(subject_dt_key),
        args.streaming,
    )
    write_table(
        subject_dt_coverage,
        outdir,
        "03_subject_dt_coverage",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 4. Dataset-wide coverage by dt
    # ------------------------------------------------------------------
    dt_coverage = collect(
        base.group_by("dt")
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_biological_subjects"),
                pl.struct(group_subject_key).n_unique().alias("n_group_subject_units"),
                pl.struct(group_subject_interval_key)
                .n_unique()
                .alias("n_group_subject_interval_pairs"),
                pl.struct(["t0", "t1"])
                .n_unique()
                .alias("n_distinct_nominal_interval_pairs"),
                pl.struct(group_subject_clone_key)
                .n_unique()
                .alias("n_group_subject_clonotype_pairs"),
                pl.col("_ab").cast(pl.Float64).mean().alias("fraction_forward_ab_eligible"),
                pl.col("_ba").cast(pl.Float64).mean().alias("fraction_forward_ba_eligible"),
                pl.col("_common4").cast(pl.Float64).mean().alias("fraction_common4"),
            ]
        )
        .sort("dt"),
        args.streaming,
    )

    interval_distribution = (
        subject_interval_coverage.group_by("dt")
        .agg(
            [
                pl.col("n_transitions")
                .median()
                .alias("median_transitions_per_group_subject_interval"),
                pl.col("n_transitions")
                .quantile(0.25)
                .alias("q25_transitions_per_group_subject_interval"),
                pl.col("n_transitions")
                .quantile(0.75)
                .alias("q75_transitions_per_group_subject_interval"),
                pl.col("n_clonotypes")
                .median()
                .alias("median_clonotypes_per_group_subject_interval"),
            ]
        )
        .sort("dt")
    )
    dt_coverage = dt_coverage.join(interval_distribution, on="dt", how="left")
    write_table(
        dt_coverage,
        outdir,
        "04_coverage_summary_by_dt",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 5. Nominal interval-pair coverage
    # ------------------------------------------------------------------
    interval_pair_coverage = collect(
        base.group_by(["t0", "t1", "dt"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_biological_subjects"),
                pl.struct(group_subject_key).n_unique().alias("n_group_subject_units"),
                pl.struct(group_subject_clone_key)
                .n_unique()
                .alias("n_group_subject_clonotype_pairs"),
                pl.col("_ab").cast(pl.Float64).mean().alias("fraction_forward_ab_eligible"),
                pl.col("_ba").cast(pl.Float64).mean().alias("fraction_forward_ba_eligible"),
                pl.col("_common4").cast(pl.Float64).mean().alias("fraction_common4"),
            ]
        )
        .sort(["dt", "t0", "t1"]),
        args.streaming,
    )
    write_table(
        interval_pair_coverage,
        outdir,
        "05_interval_pair_coverage",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 6–7. Reference observation-class composition
    # ------------------------------------------------------------------
    class_by_dt = collect(
        base.group_by(["dt", "obs_class_reference"])
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_biological_subjects"),
                pl.struct(group_subject_interval_key)
                .n_unique()
                .alias("n_group_subject_interval_pairs"),
            ]
        )
        .sort(["dt", "obs_class_reference"]),
        args.streaming,
    )
    totals_by_dt = class_by_dt.group_by("dt").agg(
        pl.col("n_transitions").sum().alias("_dt_total")
    )
    class_by_dt = (
        class_by_dt.join(totals_by_dt, on="dt", how="left")
        .with_columns(
            (pl.col("n_transitions") / pl.col("_dt_total")).alias(
                "transition_fraction"
            )
        )
        .drop("_dt_total")
        .with_columns(pl.lit(reference_alpha).alias("reference_alpha"))
    )
    write_table(
        class_by_dt,
        outdir,
        "06_reference_observation_class_composition_by_dt",
        args.write_parquet,
        args.compression,
    )

    class_overall = collect(
        base.group_by("obs_class_reference")
        .agg(
            [
                pl.len().alias("n_transitions"),
                pl.col("subject").n_unique().alias("n_biological_subjects"),
            ]
        )
        .sort("obs_class_reference"),
        args.streaming,
    ).with_columns(
        [
            (pl.col("n_transitions") / pl.lit(total_rows)).alias(
                "transition_fraction"
            ),
            pl.lit(reference_alpha).alias("reference_alpha"),
        ]
    )
    write_table(
        class_overall,
        outdir,
        "07_reference_observation_class_composition_overall",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 8–9. Primary-estimand support
    # ------------------------------------------------------------------
    support_by_dt = dt_coverage.select(
        [
            "dt",
            "n_transitions",
            "fraction_forward_ab_eligible",
            "fraction_forward_ba_eligible",
            "fraction_common4",
        ]
    )
    write_table(
        support_by_dt,
        outdir,
        "08_estimand_support_by_dt",
        args.write_parquet,
        args.compression,
    )

    ov = overall.row(0, named=True)
    support_overall = pl.DataFrame(
        [
            {
                "support": "forward_ab_eligible",
                "n_transitions": int(round(total_rows * float(ov["fraction_forward_ab_eligible"]))),
                "fraction_of_transition_universe": float(ov["fraction_forward_ab_eligible"]),
                "downstream_role": "AB: rep1 initial conditioning, rep2 displacement",
            },
            {
                "support": "forward_ba_eligible",
                "n_transitions": int(round(total_rows * float(ov["fraction_forward_ba_eligible"]))),
                "fraction_of_transition_universe": float(ov["fraction_forward_ba_eligible"]),
                "downstream_role": "BA: rep2 initial conditioning, rep1 displacement",
            },
            {
                "support": "common4",
                "n_transitions": int(round(total_rows * float(ov["fraction_common4"]))),
                "fraction_of_transition_universe": float(ov["fraction_common4"]),
                "downstream_role": "cross-replicate fluctuation covariance",
            },
        ]
    )
    # Replace rounded counts with exact sums from the source table.
    exact_support = collect(
        base.select(
            [
                pl.col("_ab").sum().alias("ab"),
                pl.col("_ba").sum().alias("ba"),
                pl.col("_common4").sum().alias("common4"),
            ]
        ),
        args.streaming,
    ).row(0, named=True)
    support_overall = support_overall.with_columns(
        pl.Series(
            "n_transitions",
            [
                int(exact_support["ab"] or 0),
                int(exact_support["ba"] or 0),
                int(exact_support["common4"] or 0),
            ],
        )
    )
    write_table(
        support_overall,
        outdir,
        "09_estimand_support_overall",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # 10–13. Production readiness audits
    # ------------------------------------------------------------------
    core_completeness = build_core_completeness(
        base,
        schema_names,
        total_rows,
        args.streaming,
    )
    write_table(
        core_completeness,
        outdir,
        "10_core_field_completeness",
        args.write_parquet,
        args.compression,
    )

    estimand_readiness = build_estimand_readiness(base, args.streaming)
    write_table(
        estimand_readiness,
        outdir,
        "11_estimand_readiness_audit",
        args.write_parquet,
        args.compression,
    )

    structural_consistency = build_structural_consistency(base, args.streaming)
    write_table(
        structural_consistency,
        outdir,
        "12_structural_consistency_audit",
        args.write_parquet,
        args.compression,
    )

    optional_completeness = build_optional_model_completeness(
        base,
        schema_names,
        total_rows,
        args.streaming,
    )
    write_table(
        optional_completeness,
        outdir,
        "13_optional_model_field_completeness",
        args.write_parquet,
        args.compression,
    )

    # Production readiness flags used in compact overview/report.
    min_core_fraction = float(
        core_completeness.get_column("fraction_valid").min()
    )
    readiness_fractions = [
        value
        for value in estimand_readiness.get_column(
            "fraction_ready_within_support"
        ).to_list()
        if value is not None and math.isfinite(float(value))
    ]
    min_estimand_readiness = (
        float(min(readiness_fractions)) if readiness_fractions else float("nan")
    )
    max_structural_mismatch = int(
        structural_consistency.get_column("n_mismatch").max() or 0
    )

    # ------------------------------------------------------------------
    # 14. Compact manuscript/supplementary overview
    # ------------------------------------------------------------------
    class_fraction_map = {
        row["obs_class_reference"]: float(row["transition_fraction"])
        for row in class_overall.iter_rows(named=True)
    }

    overview_rows: List[Dict[str, Any]] = []
    add_overview_row(
        overview_rows,
        "biological_subjects",
        ov["n_biological_subjects"],
        "count",
        "Unique biological subjects; extra grouping configurations are not counted as subjects.",
    )
    add_overview_row(
        overview_rows,
        "group_subject_units",
        ov["n_group_subject_units"],
        "count",
        "Unique extra-group × subject units; equals biological subjects when no extra groups are supplied.",
    )
    add_overview_row(
        overview_rows,
        "transitions",
        ov["n_transition_rows"],
        "count",
        "Total generic longitudinal transition rows.",
    )
    add_overview_row(
        overview_rows,
        "group_subject_clonotype_pairs",
        ov["n_group_subject_clonotype_pairs"],
        "count",
        "Unique extra-group × subject × clonotype combinations.",
    )
    add_overview_row(
        overview_rows,
        "group_subject_interval_pairs",
        ov["n_group_subject_interval_pairs"],
        "count",
        "Unique extra-group × subject × temporal interval units.",
    )
    add_overview_row(
        overview_rows,
        "distinct_nominal_interval_pairs",
        ov["n_distinct_nominal_interval_pairs"],
        "count",
        "Distinct nominal t0→t1 combinations.",
    )
    add_overview_row(overview_rows, "minimum_dt", ov["min_dt"], "weeks")
    add_overview_row(overview_rows, "maximum_dt", ov["max_dt"], "weeks")
    add_overview_row(
        overview_rows,
        "reference_alpha",
        reference_alpha,
        "probability threshold",
        "Reference operational observation threshold; alternative alpha values are evaluated downstream.",
    )
    add_overview_row(
        overview_rows,
        "forward_ab_eligible_fraction",
        ov["fraction_forward_ab_eligible"],
        "fraction",
    )
    add_overview_row(
        overview_rows,
        "forward_ba_eligible_fraction",
        ov["fraction_forward_ba_eligible"],
        "fraction",
    )
    add_overview_row(
        overview_rows,
        "common4_fraction",
        ov["fraction_common4"],
        "fraction",
    )
    for cls in OBS_CLASSES:
        add_overview_row(
            overview_rows,
            f"reference_{cls}_fraction",
            class_fraction_map.get(cls),
            "fraction",
            "Reference-threshold operational observation class.",
        )
    add_overview_row(
        overview_rows,
        "minimum_core_field_completeness",
        min_core_fraction,
        "fraction",
    )
    add_overview_row(
        overview_rows,
        "minimum_estimand_readiness_within_support",
        min_estimand_readiness,
        "fraction",
    )
    add_overview_row(
        overview_rows,
        "maximum_structural_consistency_mismatches",
        max_structural_mismatch,
        "count",
    )

    # Optional posterior propagation summary, if available.
    if "posterior_success" in schema_set:
        posterior_success = collect(
            base.select(
                bool_expr("posterior_success")
                .cast(pl.Float64)
                .mean()
                .alias("fraction")
            ),
            args.streaming,
        ).item()
        add_overview_row(
            overview_rows,
            "optional_posterior_propagation_success_fraction",
            posterior_success,
            "fraction",
            "Optional model-based transition-posterior propagation diagnostic.",
        )

    dataset_overview = pl.DataFrame(overview_rows)
    write_table(
        dataset_overview,
        outdir,
        "14_dataset_overview",
        args.write_parquet,
        args.compression,
    )

    # ------------------------------------------------------------------
    # Markdown report and manifest
    # ------------------------------------------------------------------
    report_path = args.report_md or (outdir / "transition_support_summary.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)

    report_lines = [
        "# Transition and estimand-support characterization",
        "",
        f"- generated: {datetime.now().isoformat(timespec='seconds')}",
        f"- script_version: {SCRIPT_VERSION}",
        f"- input: `{transitions_path}`",
        f"- transitions: {int(ov['n_transition_rows']):,}",
        f"- biological subjects: {int(ov['n_biological_subjects'])}",
        f"- group-subject units: {int(ov['n_group_subject_units'])}",
        f"- group-subject-clonotype pairs: {int(ov['n_group_subject_clonotype_pairs']):,}",
        f"- group-subject-interval pairs: {int(ov['n_group_subject_interval_pairs']):,}",
        f"- distinct nominal interval pairs: {int(ov['n_distinct_nominal_interval_pairs'])}",
        f"- temporal lag range: {ov['min_dt']}–{ov['max_dt']} weeks",
        f"- reference alpha: {reference_alpha:g}",
        f"- AB eligible fraction: {float(ov['fraction_forward_ab_eligible']):.4%}",
        f"- BA eligible fraction: {float(ov['fraction_forward_ba_eligible']):.4%}",
        f"- common4 fraction: {float(ov['fraction_common4']):.4%}",
        f"- minimum core-field completeness: {min_core_fraction:.4%}",
        f"- minimum estimand readiness within intended support: {min_estimand_readiness:.4%}",
        f"- maximum structural-consistency mismatches: {max_structural_mismatch}",
        "",
        "## Coverage by temporal lag",
        "",
    ]
    for row in dt_coverage.iter_rows(named=True):
        report_lines.append(
            "- Δt={dt}: {n_transitions:,} transitions; "
            "{n_group_subject_interval_pairs} group-subject interval units; "
            "AB={fraction_forward_ab_eligible:.2%}, "
            "BA={fraction_forward_ba_eligible:.2%}, "
            "common4={fraction_common4:.2%}.".format(**row)
        )

    report_lines.extend(["", "## Reference observation classes", ""])
    for row in class_overall.iter_rows(named=True):
        report_lines.append(
            "- {obs_class_reference}: {n_transitions:,} transitions "
            "({transition_fraction:.2%}).".format(**row)
        )

    report_lines.extend(
        [
            "",
            "## Scope",
            "",
            "Step 6 characterizes the generic Step-5 transition universe and "
            "its support for downstream replicate-resolved estimands. It does "
            "not estimate drift, cross-replicate covariance, temporal scaling, "
            "or threshold sensitivity.",
            "",
            "TT/TF/FT/FF are reference-threshold operational observation "
            "classes and are not biological presence/absence states.",
            "",
            "## Output tables",
            "",
        ]
    )
    for path in sorted(outdir.glob("*.csv")):
        report_lines.append(f"- `{path.name}`")
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "canonical_upstream": "5-longitudinal_transition_assembly.py",
        "input": str(transitions_path),
        "extra_group_cols": extra_group_cols,
        "biological_n_definition": "unique subject only",
        "reference_alpha": reference_alpha,
        "reference_observation_classes_are_biological_states": False,
        "step6_applies_observation_class_filtering": False,
        "step6_estimates_dynamics": False,
        "core_support": {
            "forward_ab_eligible": "rep1 positive at t0; rep2 positive at t0 and t1",
            "forward_ba_eligible": "rep2 positive at t0; rep1 positive at t0 and t1",
            "common4": "both replicates positive at both endpoints",
        },
        "primary_downstream_estimands": {
            "forward": "observed replicate-decoupled AB/BA",
            "fluctuations": "Cov(dx_observed_rep1, dx_observed_rep2) on common4",
            "fluctuation_conditioning": "xmid_latent",
        },
        "continuous_observation_statistics_preserved": ["p_value_t0", "p_value_t1"],
        "optional_model_based_fields_are_required": False,
        "production_readiness": {
            "minimum_core_field_completeness": min_core_fraction,
            "minimum_estimand_readiness_within_support": min_estimand_readiness,
            "maximum_structural_consistency_mismatches": max_structural_mismatch,
        },
    }
    manifest_path = outdir / "transition_support_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if args.figures:
        make_figures(
            dt_coverage=dt_coverage,
            subject_interval_coverage=subject_interval_coverage,
            support_by_dt=support_by_dt,
            figure_dir=outdir / "figures",
            dpi=int(args.figure_dpi),
        )

    print(f"[STEP6] Wrote outputs: {outdir}")
    print(f"[STEP6] Wrote report: {report_path}")
    print(f"[STEP6] Wrote manifest: {manifest_path}")
    print(
        "[STEP6] production_readiness: "
        f"core={min_core_fraction:.6f}; "
        f"estimand={min_estimand_readiness:.6f}; "
        f"max_mismatch={max_structural_mismatch}"
    )
    print("[DONE]")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
