#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4-multirepresentation_trajectory_assembly.py
============================================

Production assembly of longitudinal clonotype-state trajectories for the
replicate-resolved ClonoDynamics framework.

Overview
--------
This script implements Step 4 of the current ClonoDynamics pipeline.

Step 2 preserves three conceptually distinct state layers for every
subject-clonotype-time state:

    1. model-based latent consensus

       A paired-replicate posterior summary of clonotype abundance under the
       Negative-Binomial observation model. The canonical point estimate is
       `x_latent_median` (and `f_latent_median`), with posterior uncertainty
       retained separately. This layer is used downstream for conditioning
       geometry, uncertainty characterization and latent-representation
       sensitivity analyses. It is not treated as ground-truth abundance.

    2. replicate-resolved observed measurements

       The two technical replicates are kept separate. Canonical observed
       relative frequencies are reconstructed directly from replicate counts
       and replicate sequencing depths:

           f_observed_rep1 = count_rep1 / depth_rep1
           f_observed_rep2 = count_rep2 / depth_rep2

       and, only for positive counts,

           x_observed_rep1 = log(f_observed_rep1)
           x_observed_rep2 = log(f_observed_rep2).

       No pseudocount is added. Therefore an observed zero remains a realized
       zero-count observation and its observed log-frequency is null rather
       than an arbitrarily imputed value.

       These replicate-resolved observed measurements provide the measurement
       layer used by the primary downstream replicate-decoupled forward-drift
       and cross-replicate fluctuation estimands.

    3. observation annotations

       Realized replicate positivity, empirical p-value-based operational
       observability, and posterior-predictive detectability/dropout are
       retained as separate quantities. They are not interpreted as
       interchangeable measures of biological presence.

Step 4 assembles these layers into one standardized long-format state table
with one row for each

    subject x aaSeqCDR3 x time

combination, optionally extended by explicit grouping variables for pseudo,
null or experimental workflows.

Step 4 is an ASSEMBLY step. It does not refit the observation model, does not
alter posterior estimates, does not estimate dynamics, and does not apply
abundance-, posterior-quality-, detectability- or observation-class filters.

Position in the pipeline
------------------------

    Step 1
    1-repertoire_characterization.py
        -> descriptive observed repertoire architecture
        -> no parameter transfer to Step 2

    Step 2
    2-multirepresentation_clonotype_state_inference.py
        -> model-based latent consensus + uncertainty
        -> replicate-resolved observed counts/depths
        -> p_value / operational observability
        -> posterior-predictive detectability

    Step 3
    3-latent_state_and_observation_model_diagnostics.py
        -> diagnostic characterization only
        -> no state filtering and no transformation of Step-2 states

    Step 4
    4-multirepresentation_trajectory_assembly.py
        -> standardized multirepresentation longitudinal state table

    Step 5
        -> generic endpoint pairing / transition assembly

    Step 6
        -> transition-support and eligibility characterization

    Step 7+
        -> estimand-specific pseudo benchmarks and longitudinal dynamics

Canonical input
---------------
The canonical production input is the combined Step-2 table:

    per_clone_latent_subject.parquet

from:

    2-multirepresentation_clonotype_state_inference.py

The script accepts either:

    --per-clone <single table>

or:

    --per-clone-dir <directory of table parts>.

Supported formats are Parquet/PQ, Feather/Arrow, CSV, TSV and TXT. Parquet is
recommended for production-scale data.

Required production fields
--------------------------
Identifiers:

    subject
    time
    aaSeqCDR3

Model-based consensus:

    f_latent_median
    x_latent_median

Replicate-resolved observed measurements:

    count_rep1
    count_rep2
    depth_rep1
    depth_rep2

Observation annotations:

    p_value
    observable
    p_detect_state
    p_dropout_state

Column names can be overridden explicitly from the CLI. There is no automatic
cross-representation fallback: an observed abundance column is never silently
used as a latent abundance column, and vice versa.

State-key uniqueness
--------------------
The canonical Step-2 production table must contain exactly one row for each
state key:

    subject x aaSeqCDR3 x time

or, when extra grouping columns are supplied:

    <extra groups> x subject x aaSeqCDR3 x time.

Duplicate state keys are treated as an upstream data-contract violation and
cause Step 4 to stop. Values from duplicate rows are NOT averaged, summed or
otherwise collapsed in production.

This differs intentionally from older Step-4 implementations that silently
aggregated duplicate state rows.

Model-based latent consensus
----------------------------
The canonical latent point estimate is propagated without re-estimation:

    f_latent = f_latent_median
    x_latent = x_latent_median.

When present, Step 4 also carries forward:

    f_latent_q025
    f_latent_q975
    x_latent_q025
    x_latent_q975
    x_latent_sd
    posterior_entropy
    logP.

For backward compatibility with the immediately preceding transition builder,
`freq_latent` is written as an alias of `f_latent`. Generic aliases such as
`freq` and `log_freq` are deliberately not created because they obscure which
representation is being used.

Replicate-resolved observed representation
-------------------------------------------
Observed relative frequencies are reconstructed from the raw replicate
measurement opportunity:

    f_observed_rep1 = count_rep1 / depth_rep1
    f_observed_rep2 = count_rep2 / depth_rep2.

The corresponding observed log-frequencies are defined only for positive
counts:

    x_observed_rep1 = log(f_observed_rep1),  if count_rep1 > 0
    x_observed_rep2 = log(f_observed_rep2),  if count_rep2 > 0.

For a zero count, the corresponding `x_observed_rep*` value is null. Step 4
never adds a pseudocount to create an observed log-frequency.

The table also contains:

    read_positive_rep1
    read_positive_rep2
    both_replicates_positive
    n_positive_replicates
    positive_state_class
    positive_replicate_source.

The technical-detection support classes are defined as:

    TWO_POSITIVE    both replicate counts > 0
    ONE_POSITIVE    exactly one replicate count > 0
    ZERO_POSITIVE   neither replicate count > 0.

`ZERO_POSITIVE` is not expected in the canonical Step-2 union of observed
clonotypes, but is represented explicitly if encountered.

Auxiliary positive-only representation
--------------------------------------
For compatibility with Step-3 representation diagnostics, Step 4 also
constructs the auxiliary positive-only state directly from the replicate
observations:

    TWO_POSITIVE:
        x_positive_current = (x_observed_rep1 + x_observed_rep2) / 2
        f_positive_current = geometric mean of the two observed frequencies

    ONE_POSITIVE:
        x_positive_current = observed log-frequency of the positive replicate
        f_positive_current = observed frequency of the positive replicate.

This auxiliary representation is not the primary downstream fluctuation or
forward-drift measurement layer.

Operational observability
-------------------------
The empirical `p_value` is preserved for every state so downstream threshold
sensitivity can reconstruct:

    observable(alpha) = 1[p_value < alpha]

on exactly the same state and transition universe.

Step 4 writes a canonical reference annotation:

    observable_reference

using the explicit CLI value:

    --reference-alpha

(default 0.05).

The legacy output name `observable` is retained as an alias of
`observable_reference` for compatibility. When upstream `observable` is
present, Step 4 validates it against `p_value < reference_alpha` by default and
stops if the two disagree. This check can be disabled explicitly with
`--no-validate-observable` for legacy reconstruction only.

Model-based detectability
-------------------------
The posterior-predictive state-level quantities:

    p_detect_state
    p_dropout_state

are propagated unchanged. Replicate-specific posterior-predictive quantities
are also retained when available:

    p_detect_rep1
    p_detect_rep2
    p_dropout_rep1
    p_dropout_rep2.

These are model-based observation quantities and remain distinct from realized
read positivity and p-value-based operational observability.

Replicate-specific latent summaries
------------------------------------
Single-replicate latent posterior summaries emitted by `noiseK_latent.py` are
retained when present. They are diagnostic/sensitivity representations. They
are not the primary observed measurement layer used for replicate-decoupled
forward drift or cross-replicate fluctuation covariance.

Filtering policy
----------------
No state is removed on the basis of:

    abundance;
    posterior width;
    posterior entropy;
    replicate agreement;
    ONE_POSITIVE/TWO_POSITIVE status;
    p_detect_state;
    p_value;
    operational observability.

Instead, the script validates the canonical state contract. Invalid identifiers,
negative counts, non-positive sequencing depths, invalid latent point estimates,
or p-values outside [0,1] cause an explicit error rather than silent filtering.

Canonical output
----------------
The production output is:

    multirepresentation_trajectories_long.parquet

with one row per standardized state key.

Core fields include:

    subject
    aaSeqCDR3
    time

    f_latent
    freq_latent
    x_latent

    count_rep1
    count_rep2
    depth_rep1
    depth_rep2

    f_observed_rep1
    f_observed_rep2
    x_observed_rep1
    x_observed_rep2

    read_positive_rep1
    read_positive_rep2
    both_replicates_positive
    n_positive_replicates
    positive_state_class
    positive_replicate_source

    f_positive_current
    x_positive_current

    p_value
    reference_alpha
    observable_reference
    observable

    p_detect_state
    p_dropout_state.

Optional posterior and replicate-specific latent fields are appended when
available.

Additional outputs
------------------

    included_subject_clones.parquet
        unique subject x aaSeqCDR3 pairs, extended by extra grouping columns;

    trajectory_assembly_report.md
        concise human-readable construction report;

    trajectory_assembly_manifest.json
        machine-readable state-representation contract and provenance.

Role in downstream inference
----------------------------
Step 4 does not define any dynamical estimand.

The current downstream contract is:

    primary forward drift
        observed replicate-resolved AB/BA estimator;

    primary fluctuations
        cross-replicate covariance of observed displacements on common4;

    primary fluctuation conditioning
        xmid_latent;

    latent single-replicate / xstar / positive-only representations
        diagnostic or sensitivity roles;

    operational observation classes
        reconstructed from p_value at the requested threshold.

Thus Step 4 preserves the independent measurement channels required by the
replicate-resolved framework instead of collapsing them into one universal
abundance representation.

Compatible with Python >= 3.9 and Polars 1.36.x.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence

try:
    import polars as pl
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "ERROR: polars is required. Install it in the active environment with:\n"
        "  pip install polars pyarrow\n"
    ) from exc


SCRIPT_VERSION = "v12-multirepresentation-trajectory-contract-2026-09-09"

PARQUET_EXTS = {".parquet", ".pq"}
CSV_EXTS = {".csv", ".tsv", ".txt"}
FEATHER_EXTS = {".feather", ".arrow"}
TRUE_STRINGS = {"true", "t", "1", "yes", "y", "observable", "observed"}
FALSE_STRINGS = {"false", "f", "0", "no", "n", "unobservable", "not_observed"}

SINGLE_REP_LATENT_COLS = [
    "rep1_f_latent_mean", "rep1_f_latent_median", "rep1_f_latent_mode",
    "rep1_f_latent_q025", "rep1_f_latent_q975",
    "rep1_x_latent_mean", "rep1_x_latent_median", "rep1_x_latent_mode",
    "rep1_x_latent_q025", "rep1_x_latent_q975", "rep1_x_latent_sd",
    "rep1_posterior_entropy",
    "rep2_f_latent_mean", "rep2_f_latent_median", "rep2_f_latent_mode",
    "rep2_f_latent_q025", "rep2_f_latent_q975",
    "rep2_x_latent_mean", "rep2_x_latent_median", "rep2_x_latent_mode",
    "rep2_x_latent_q025", "rep2_x_latent_q975", "rep2_x_latent_sd",
    "rep2_posterior_entropy",
]

SINGLE_REP_ALIAS_COLS = [
    "f_rep1_latent_median", "f_rep2_latent_median",
    "x_rep1_latent_median", "x_rep2_latent_median",
    "x_rep1_latent_q025", "x_rep1_latent_q975",
    "x_rep2_latent_q025", "x_rep2_latent_q975",
]

OPTIONAL_MODEL_COLUMNS = [
    "f_latent_q025",
    "f_latent_q975",
    "x_latent_q025",
    "x_latent_q975",
    "x_latent_sd",
    "posterior_entropy",
    "logP",
    "p_detect_rep1",
    "p_detect_rep2",
    "p_dropout_rep1",
    "p_dropout_rep2",
    "replicate_posterior_bhattacharyya",
]


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def split_csv(text: Optional[str]) -> list[str]:
    if text is None:
        return []
    return [item.strip() for item in str(text).split(",") if item.strip()]


def collect_streaming(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)


def collect_schema_names(lf: pl.LazyFrame) -> list[str]:
    try:
        return list(lf.collect_schema().names())
    except Exception:
        return list(lf.schema.keys())


def file_kind(path: Path) -> str:
    suffix = Path(path).suffix.lower()
    if suffix in PARQUET_EXTS:
        return "parquet"
    if suffix in FEATHER_EXTS:
        return "feather"
    if suffix in CSV_EXTS:
        return "csv"
    raise ValueError(f"Unsupported table extension: {path}")


def list_table_files(directory: Path) -> list[Path]:
    directory = Path(directory)
    files: list[Path] = []
    for pattern in (
        "*.parquet", "*.pq", "*.feather", "*.arrow", "*.csv", "*.tsv", "*.txt"
    ):
        files.extend(sorted(directory.glob(pattern)))
    return sorted(set(files))


def scan_one(path: Path, sep: str = "\t") -> pl.LazyFrame:
    path = Path(path)
    kind = file_kind(path)
    if kind == "parquet":
        return pl.scan_parquet(str(path))
    if kind == "feather":
        return pl.read_ipc(str(path)).lazy()
    delimiter = "," if path.suffix.lower() == ".csv" else sep
    return pl.scan_csv(
        str(path),
        separator=delimiter,
        infer_schema_length=10000,
    )


def scan_input(
    per_clone: Optional[Path],
    per_clone_dir: Optional[Path],
    sep: str,
) -> tuple[pl.LazyFrame, str, list[Path]]:
    if per_clone is None and per_clone_dir is None:
        raise ValueError("Provide either --per-clone or --per-clone-dir")
    if per_clone is not None and per_clone_dir is not None:
        raise ValueError("Use only one input mode: --per-clone OR --per-clone-dir")

    if per_clone is not None:
        path = Path(per_clone)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")
        return scan_one(path, sep=sep), str(path), [path]

    assert per_clone_dir is not None
    directory = Path(per_clone_dir)
    if not directory.exists():
        raise FileNotFoundError(f"Input directory not found: {directory}")

    files = list_table_files(directory)
    if not files:
        raise FileNotFoundError(f"No supported table files found in: {directory}")

    kinds = {file_kind(path) for path in files}
    if kinds == {"parquet"}:
        return (
            pl.scan_parquet([str(path) for path in files]),
            str(directory),
            files,
        )

    lfs = [scan_one(path, sep=sep) for path in files]
    return pl.concat(lfs, how="diagonal_relaxed"), str(directory), files


def require_columns(columns: Sequence[str], required: Sequence[str]) -> None:
    names = set(columns)
    missing = [column for column in required if column not in names]
    if missing:
        raise ValueError(
            "Input Step-2 table is missing required production columns: "
            + ", ".join(missing)
        )


def bool_from_any_expr(column: str) -> pl.Expr:
    """Normalize an upstream Boolean-like field without silently inventing values."""
    text = (
        pl.col(column)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    return (
        pl.when(text.is_in(sorted(TRUE_STRINGS)))
        .then(pl.lit(True))
        .when(text.is_in(sorted(FALSE_STRINGS)))
        .then(pl.lit(False))
        .otherwise(pl.lit(None, dtype=pl.Boolean))
    )


def validate_core_rows(
    lf: pl.LazyFrame,
    args: argparse.Namespace,
    extra_group_cols: Sequence[str],
) -> None:
    """Fail on invalid production rows instead of silently filtering them."""
    group_null_exprs = [pl.col(column).is_null() for column in extra_group_cols]
    any_group_null = (
        pl.any_horizontal(group_null_exprs)
        if group_null_exprs
        else pl.lit(False)
    )

    summary = collect_streaming(
        lf.select(
            [
                pl.len().alias("n_rows"),
                (
                    pl.col(args.subject_col).is_null()
                    | pl.col(args.time_col).is_null()
                    | pl.col(args.clone_col).is_null()
                    | (
                        pl.col(args.clone_col)
                        .cast(pl.Utf8, strict=False)
                        .str.len_chars()
                        == 0
                    )
                    | any_group_null
                ).sum().alias("invalid_identifier_rows"),
                (
                    (pl.col(args.count_rep1_col).cast(pl.Float64, strict=False) < 0)
                    | (pl.col(args.count_rep2_col).cast(pl.Float64, strict=False) < 0)
                    | pl.col(args.count_rep1_col).cast(pl.Float64, strict=False).is_null()
                    | pl.col(args.count_rep2_col).cast(pl.Float64, strict=False).is_null()
                ).sum().alias("invalid_count_rows"),
                (
                    (pl.col(args.depth_rep1_col).cast(pl.Float64, strict=False) <= 0)
                    | (pl.col(args.depth_rep2_col).cast(pl.Float64, strict=False) <= 0)
                    | pl.col(args.depth_rep1_col).cast(pl.Float64, strict=False).is_null()
                    | pl.col(args.depth_rep2_col).cast(pl.Float64, strict=False).is_null()
                ).sum().alias("invalid_depth_rows"),
                (
                    pl.col(args.latent_f_col).cast(pl.Float64, strict=False).is_null()
                    | ~pl.col(args.latent_f_col).cast(pl.Float64, strict=False).is_finite()
                    | (pl.col(args.latent_f_col).cast(pl.Float64, strict=False) <= 0)
                    | pl.col(args.latent_x_col).cast(pl.Float64, strict=False).is_null()
                    | ~pl.col(args.latent_x_col).cast(pl.Float64, strict=False).is_finite()
                ).sum().alias("invalid_latent_rows"),
                (
                    pl.col(args.p_value_col).cast(pl.Float64, strict=False).is_null()
                    | ~pl.col(args.p_value_col).cast(pl.Float64, strict=False).is_finite()
                    | (pl.col(args.p_value_col).cast(pl.Float64, strict=False) < 0)
                    | (pl.col(args.p_value_col).cast(pl.Float64, strict=False) > 1)
                ).sum().alias("invalid_p_value_rows"),
                (
                    pl.col(args.p_detect_state_col).cast(pl.Float64, strict=False).is_null()
                    | ~pl.col(args.p_detect_state_col).cast(pl.Float64, strict=False).is_finite()
                    | (pl.col(args.p_detect_state_col).cast(pl.Float64, strict=False) < 0)
                    | (pl.col(args.p_detect_state_col).cast(pl.Float64, strict=False) > 1)
                    | pl.col(args.p_dropout_state_col).cast(pl.Float64, strict=False).is_null()
                    | ~pl.col(args.p_dropout_state_col).cast(pl.Float64, strict=False).is_finite()
                    | (pl.col(args.p_dropout_state_col).cast(pl.Float64, strict=False) < 0)
                    | (pl.col(args.p_dropout_state_col).cast(pl.Float64, strict=False) > 1)
                ).sum().alias("invalid_detectability_rows"),
                bool_from_any_expr(args.observable_col)
                .is_null()
                .sum()
                .alias("invalid_observable_rows"),
            ]
        )
    ).row(0, named=True)

    problems = {
        key: int(value or 0)
        for key, value in summary.items()
        if key != "n_rows" and int(value or 0) > 0
    }
    if problems:
        raise ValueError(
            "Step-4 production state validation failed: "
            + ", ".join(f"{key}={value}" for key, value in problems.items())
        )


def validate_unique_state_keys(
    lf: pl.LazyFrame,
    args: argparse.Namespace,
    extra_group_cols: Sequence[str],
) -> None:
    key_exprs = [pl.col(column) for column in extra_group_cols]
    key_exprs.extend(
        [
            pl.col(args.subject_col).alias("subject"),
            pl.col(args.clone_col).alias("aaSeqCDR3"),
            pl.col(args.time_col).alias("time"),
        ]
    )
    key_names = list(extra_group_cols) + ["subject", "aaSeqCDR3", "time"]

    duplicates = collect_streaming(
        lf.select(key_exprs)
        .group_by(key_names)
        .agg(pl.len().alias("n_rows"))
        .filter(pl.col("n_rows") > 1)
        .sort("n_rows", descending=True)
        .head(10)
    )

    if not duplicates.is_empty():
        raise ValueError(
            "Duplicate Step-2 state keys detected. Production Step 4 does not "
            "aggregate duplicate subject-clonotype-time states. Examples:\n"
            + str(duplicates)
        )


def validate_reference_observability(
    lf: pl.LazyFrame,
    args: argparse.Namespace,
) -> None:
    if not args.validate_observable:
        return

    expected = (
        pl.col(args.p_value_col).cast(pl.Float64, strict=False)
        < pl.lit(float(args.reference_alpha))
    )
    upstream = bool_from_any_expr(args.observable_col)

    result = collect_streaming(
        lf.select(
            (expected != upstream)
            .fill_null(True)
            .sum()
            .alias("n_mismatches")
        )
    )
    mismatches = int(result["n_mismatches"][0] or 0)
    if mismatches:
        raise ValueError(
            f"Upstream {args.observable_col!r} disagrees with "
            f"{args.p_value_col} < {args.reference_alpha:g} for "
            f"{mismatches:,} states. Use the correct --reference-alpha or "
            "--no-validate-observable only for explicit legacy reconstruction."
        )


def build_trajectory_lazy(
    lf: pl.LazyFrame,
    columns: Sequence[str],
    args: argparse.Namespace,
    extra_group_cols: Sequence[str],
    sort_output: bool,
) -> tuple[pl.LazyFrame, list[str]]:
    names = set(columns)

    required = [
        args.subject_col,
        args.time_col,
        args.clone_col,
        args.latent_f_col,
        args.latent_x_col,
        args.count_rep1_col,
        args.count_rep2_col,
        args.depth_rep1_col,
        args.depth_rep2_col,
        args.p_value_col,
        args.observable_col,
        args.p_detect_state_col,
        args.p_dropout_state_col,
    ]
    require_columns(columns, required)

    extra_missing = [column for column in extra_group_cols if column not in names]
    if extra_missing:
        raise ValueError(
            "--extra-group-cols not found in input: " + ", ".join(extra_missing)
        )

    validate_core_rows(lf, args, extra_group_cols)
    validate_unique_state_keys(lf, args, extra_group_cols)
    validate_reference_observability(lf, args)

    count1 = pl.col(args.count_rep1_col).cast(pl.Float64, strict=False)
    count2 = pl.col(args.count_rep2_col).cast(pl.Float64, strict=False)
    depth1 = pl.col(args.depth_rep1_col).cast(pl.Float64, strict=False)
    depth2 = pl.col(args.depth_rep2_col).cast(pl.Float64, strict=False)

    f1 = count1 / depth1
    f2 = count2 / depth2
    pos1 = count1 > 0
    pos2 = count2 > 0
    x1 = pl.when(pos1).then(f1.log()).otherwise(pl.lit(None, dtype=pl.Float64))
    x2 = pl.when(pos2).then(f2.log()).otherwise(pl.lit(None, dtype=pl.Float64))

    n_positive = pos1.cast(pl.Int8) + pos2.cast(pl.Int8)
    positive_class = (
        pl.when(n_positive == 2)
        .then(pl.lit("TWO_POSITIVE"))
        .when(n_positive == 1)
        .then(pl.lit("ONE_POSITIVE"))
        .otherwise(pl.lit("ZERO_POSITIVE"))
    )
    positive_source = (
        pl.when(pos1 & pos2)
        .then(pl.lit("BOTH"))
        .when(pos1 & ~pos2)
        .then(pl.lit("REP1_ONLY"))
        .when(~pos1 & pos2)
        .then(pl.lit("REP2_ONLY"))
        .otherwise(pl.lit("NONE"))
    )

    f_positive = (
        pl.when(pos1 & pos2)
        .then((f1 * f2).sqrt())
        .when(pos1 & ~pos2)
        .then(f1)
        .when(~pos1 & pos2)
        .then(f2)
        .otherwise(pl.lit(None, dtype=pl.Float64))
    )
    x_positive = (
        pl.when(pos1 & pos2)
        .then(0.5 * (x1 + x2))
        .when(pos1 & ~pos2)
        .then(x1)
        .when(~pos1 & pos2)
        .then(x2)
        .otherwise(pl.lit(None, dtype=pl.Float64))
    )

    p_value = pl.col(args.p_value_col).cast(pl.Float64, strict=False)
    observable_reference = p_value < pl.lit(float(args.reference_alpha))

    select_exprs: list[pl.Expr] = []
    for column in extra_group_cols:
        select_exprs.append(pl.col(column))

    select_exprs.extend(
        [
            pl.col(args.subject_col).cast(pl.Int64, strict=False).alias("subject"),
            pl.col(args.clone_col).cast(pl.Utf8, strict=False).alias("aaSeqCDR3"),
            pl.col(args.time_col).cast(pl.Int64, strict=False).alias("time"),
            pl.col(args.latent_f_col).cast(pl.Float64, strict=False).alias("f_latent"),
            pl.col(args.latent_f_col).cast(pl.Float64, strict=False).alias("freq_latent"),
            pl.col(args.latent_x_col).cast(pl.Float64, strict=False).alias("x_latent"),
            pl.col(args.count_rep1_col).cast(pl.Int64, strict=False).alias("count_rep1"),
            pl.col(args.count_rep2_col).cast(pl.Int64, strict=False).alias("count_rep2"),
            pl.col(args.depth_rep1_col).cast(pl.Int64, strict=False).alias("depth_rep1"),
            pl.col(args.depth_rep2_col).cast(pl.Int64, strict=False).alias("depth_rep2"),
            f1.alias("f_observed_rep1"),
            f2.alias("f_observed_rep2"),
            x1.alias("x_observed_rep1"),
            x2.alias("x_observed_rep2"),
            pos1.alias("read_positive_rep1"),
            pos2.alias("read_positive_rep2"),
            (pos1 & pos2).alias("both_replicates_positive"),
            n_positive.cast(pl.Int8).alias("n_positive_replicates"),
            positive_class.alias("positive_state_class"),
            positive_source.alias("positive_replicate_source"),
            f_positive.alias("f_positive_current"),
            x_positive.alias("x_positive_current"),
            p_value.alias("p_value"),
            pl.lit(float(args.reference_alpha), dtype=pl.Float64).alias("reference_alpha"),
            observable_reference.alias("observable_reference"),
            observable_reference.alias("observable"),
            pl.col(args.p_detect_state_col)
            .cast(pl.Float64, strict=False)
            .alias("p_detect_state"),
            pl.col(args.p_dropout_state_col)
            .cast(pl.Float64, strict=False)
            .alias("p_dropout_state"),
        ]
    )

    if "pair_id" in names:
        select_exprs.insert(
            len(extra_group_cols) + 3,
            pl.col("pair_id").cast(pl.Utf8, strict=False).alias("pair_id"),
        )
    else:
        pair_id_expr = pl.concat_str(
            [
                pl.col(args.subject_col).cast(pl.Utf8, strict=False),
                pl.col(args.time_col).cast(pl.Utf8, strict=False),
            ],
            separator="_",
        ).alias("pair_id")
        select_exprs.insert(len(extra_group_cols) + 3, pair_id_expr)

    for column in OPTIONAL_MODEL_COLUMNS:
        if column in names:
            if column == "logP":
                select_exprs.append(
                    pl.col(column).cast(pl.Float64, strict=False).alias(column)
                )
            else:
                select_exprs.append(
                    pl.col(column).cast(pl.Float64, strict=False).alias(column)
                )

    for column in SINGLE_REP_LATENT_COLS + SINGLE_REP_ALIAS_COLS:
        if column in names:
            select_exprs.append(
                pl.col(column).cast(pl.Float64, strict=False).alias(column)
            )

    out = lf.select(select_exprs)

    final_order = list(extra_group_cols) + [
        "subject",
        "aaSeqCDR3",
        "time",
        "pair_id",
        "f_latent",
        "freq_latent",
        "x_latent",
    ]

    optional_consensus = [
        "f_latent_q025",
        "f_latent_q975",
        "x_latent_q025",
        "x_latent_q975",
        "x_latent_sd",
        "posterior_entropy",
        "logP",
    ]
    final_order.extend([column for column in optional_consensus if column in names])

    final_order.extend(
        [
            "count_rep1",
            "count_rep2",
            "depth_rep1",
            "depth_rep2",
            "f_observed_rep1",
            "f_observed_rep2",
            "x_observed_rep1",
            "x_observed_rep2",
            "read_positive_rep1",
            "read_positive_rep2",
            "both_replicates_positive",
            "n_positive_replicates",
            "positive_state_class",
            "positive_replicate_source",
            "f_positive_current",
            "x_positive_current",
            "p_value",
            "reference_alpha",
            "observable_reference",
            "observable",
            "p_detect_state",
            "p_dropout_state",
        ]
    )

    for column in [
        "p_detect_rep1",
        "p_detect_rep2",
        "p_dropout_rep1",
        "p_dropout_rep2",
        "replicate_posterior_bhattacharyya",
    ]:
        if column in names:
            final_order.append(column)

    for column in SINGLE_REP_LATENT_COLS + SINGLE_REP_ALIAS_COLS:
        if column in names and column not in final_order:
            final_order.append(column)

    out_schema = set(collect_schema_names(out))
    final_order = [column for column in final_order if column in out_schema]
    out = out.select(final_order)

    if sort_output:
        out = out.sort(list(extra_group_cols) + ["subject", "aaSeqCDR3", "time"])

    key_cols = list(extra_group_cols) + ["subject", "aaSeqCDR3", "time"]
    return out, key_cols


def write_lazy_or_collect(
    lf: pl.LazyFrame,
    out_path: Path,
    output_format: str,
    parquet_compression: str,
    streaming: bool,
) -> Path:
    ensure_dir(out_path.parent)
    fmt = output_format.lower()

    if fmt == "parquet":
        out = out_path.with_suffix(".parquet")
        if hasattr(lf, "sink_parquet"):
            try:
                lf.sink_parquet(
                    str(out),
                    compression=parquet_compression,
                    statistics=True,
                    mkdir=True,
                )
                return out
            except TypeError:
                try:
                    lf.sink_parquet(str(out), compression=parquet_compression)
                    return out
                except TypeError:
                    lf.sink_parquet(str(out))
                    return out
        df = collect_streaming(lf) if streaming else lf.collect()
        df.write_parquet(str(out), compression=parquet_compression)
        return out

    df = collect_streaming(lf) if streaming else lf.collect()
    if fmt == "csv":
        out = out_path.with_suffix(".csv")
        df.write_csv(str(out))
        return out
    if fmt == "feather":
        out = out_path.with_suffix(".feather")
        df.write_ipc(str(out))
        return out
    raise ValueError(f"Unknown output format: {output_format}")


def write_partitioned(
    lf: pl.LazyFrame,
    out_root: Path,
    partition_cols: Sequence[str],
    output_format: str,
    parquet_compression: str,
    streaming: bool,
) -> tuple[Path, int]:
    """Write one file per partition after trajectory assembly.

    This path collects the already assembled state table. For large production
    datasets prefer one unpartitioned Parquet output unless partitioning is
    required explicitly for pseudo/null workflows.
    """
    out_root = ensure_dir(out_root)
    df = collect_streaming(lf) if streaming else lf.collect()
    if df.is_empty():
        return out_root, 0

    n_partitions = 0
    for key, sub in df.partition_by(list(partition_cols), as_dict=True).items():
        if not isinstance(key, tuple):
            key = (key,)
        name_parts = [
            f"{column}={value}" for column, value in zip(partition_cols, key)
        ]
        subdir = (
            ensure_dir(out_root.joinpath(*name_parts[:-1]))
            if len(name_parts) > 1
            else out_root
        )
        stem = "trajectories_" + "__".join(
            str(item).replace("/", "-") for item in name_parts
        )
        path_base = subdir / stem

        if output_format == "parquet":
            sub.write_parquet(
                str(path_base.with_suffix(".parquet")),
                compression=parquet_compression,
            )
        elif output_format == "csv":
            sub.write_csv(str(path_base.with_suffix(".csv")))
        elif output_format == "feather":
            sub.write_ipc(str(path_base.with_suffix(".feather")))
        else:
            raise ValueError(output_format)
        n_partitions += 1

    return out_root, n_partitions


def write_report(path: Path, payload: dict) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Multirepresentation trajectory assembly report\n\n")
        for key, value in payload.items():
            if isinstance(value, (list, dict, tuple)):
                value = json.dumps(value, ensure_ascii=False)
            handle.write(f"- {key}: {value}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Assemble model-based, replicate-resolved observed and observation-"
            "annotation state layers into long-format ClonoDynamics trajectories."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--per-clone", type=Path, default=None)
    parser.add_argument("--per-clone-dir", type=Path, default=None)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/4-multirepresentation_trajectory_assembly"),
    )
    parser.add_argument(
        "--output-format",
        choices=["parquet", "csv", "feather"],
        default="parquet",
    )
    parser.add_argument("--parquet-compression", default="zstd")
    parser.add_argument("--file-sep", default="\t")
    parser.add_argument("--streaming", action="store_true")

    parser.add_argument("--subject-col", default="subject")
    parser.add_argument("--time-col", default="time")
    parser.add_argument("--clone-col", default="aaSeqCDR3")

    parser.add_argument("--latent-f-col", default="f_latent_median")
    parser.add_argument("--latent-x-col", default="x_latent_median")
    parser.add_argument("--count-rep1-col", default="count_rep1")
    parser.add_argument("--count-rep2-col", default="count_rep2")
    parser.add_argument("--depth-rep1-col", default="depth_rep1")
    parser.add_argument("--depth-rep2-col", default="depth_rep2")
    parser.add_argument("--p-value-col", default="p_value")
    parser.add_argument("--observable-col", default="observable")
    parser.add_argument("--p-detect-state-col", default="p_detect_state")
    parser.add_argument("--p-dropout-state-col", default="p_dropout_state")

    parser.add_argument(
        "--reference-alpha",
        type=float,
        default=0.05,
        help="Reference operational observation threshold used for observable_reference.",
    )
    parser.add_argument(
        "--validate-observable",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Validate upstream observable against p_value < reference_alpha. "
            "Disable only for explicit legacy reconstruction."
        ),
    )

    parser.add_argument(
        "--extra-group-cols",
        default="",
        help="Comma-separated grouping columns, e.g. experiment,configuration_id.",
    )
    parser.add_argument(
        "--partition-by",
        default="",
        help="Comma-separated partition columns; must also be extra group columns.",
    )
    parser.add_argument("--no-sort", action="store_true")
    parser.add_argument("--no-included-clones", action="store_true")
    parser.add_argument("--write-csv-compat", action="store_true")
    parser.add_argument("--skip-summary-stats", action="store_true")
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if not (0.0 < float(args.reference_alpha) < 1.0):
        raise ValueError("--reference-alpha must be strictly between 0 and 1")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    validate_arguments(args)

    outdir = ensure_dir(args.outdir)
    extra_group_cols = split_csv(args.extra_group_cols)
    partition_cols = split_csv(args.partition_by)

    if partition_cols:
        missing_partitions = [
            column for column in partition_cols if column not in extra_group_cols
        ]
        if missing_partitions:
            raise ValueError(
                "--partition-by columns must also be listed in --extra-group-cols: "
                + ", ".join(missing_partitions)
            )

    print(f"[STEP4] {SCRIPT_VERSION}")
    lf0, input_label, input_files = scan_input(
        args.per_clone,
        args.per_clone_dir,
        sep=args.file_sep,
    )
    columns = collect_schema_names(lf0)
    print(f"[STEP4] Input: {input_label}")
    print(f"[STEP4] Input files: {len(input_files)}")
    print(
        "[STEP4] Extra group columns:",
        extra_group_cols if extra_group_cols else "<none>",
    )

    trajectories, key_cols = build_trajectory_lazy(
        lf=lf0,
        columns=columns,
        args=args,
        extra_group_cols=extra_group_cols,
        sort_output=not bool(args.no_sort),
    )

    output_base = outdir / "multirepresentation_trajectories_long"
    n_partitions = 0
    if partition_cols:
        written_path, n_partitions = write_partitioned(
            trajectories,
            outdir / "multirepresentation_trajectory_parts",
            partition_cols=partition_cols,
            output_format=args.output_format,
            parquet_compression=args.parquet_compression,
            streaming=bool(args.streaming),
        )
        print(
            f"[STEP4] Wrote partitioned trajectories: {written_path} "
            f"partitions={n_partitions}"
        )
    else:
        written_path = write_lazy_or_collect(
            trajectories,
            output_base,
            output_format=args.output_format,
            parquet_compression=args.parquet_compression,
            streaming=bool(args.streaming),
        )
        print(f"[STEP4] Wrote trajectories: {written_path}")

    included_path = None
    if not args.no_included_clones:
        clone_key = list(extra_group_cols) + ["subject", "aaSeqCDR3"]
        included = trajectories.select(clone_key).unique()
        if not args.no_sort:
            included = included.sort(clone_key)
        included_path = write_lazy_or_collect(
            included,
            outdir / "included_subject_clones",
            output_format=args.output_format,
            parquet_compression=args.parquet_compression,
            streaming=bool(args.streaming),
        )
        print(f"[STEP4] Wrote included subject-clones: {included_path}")

    if args.write_csv_compat and args.output_format != "csv" and not partition_cols:
        csv_path = write_lazy_or_collect(
            trajectories,
            output_base,
            output_format="csv",
            parquet_compression=args.parquet_compression,
            streaming=bool(args.streaming),
        )
        print(f"[STEP4] Wrote CSV compatibility copy: {csv_path}")

    if args.skip_summary_stats:
        stats = {}
        print("[STEP4] Skipped summary statistics")
    else:
        try:
            stats = collect_streaming(
                trajectories.select(
                    [
                        pl.len().alias("n_state_rows"),
                        pl.col("subject").n_unique().alias("n_subjects"),
                        pl.struct(list(extra_group_cols) + ["subject", "aaSeqCDR3"])
                        .n_unique()
                        .alias("n_unique_subject_clone_pairs"),
                        pl.col("read_positive_rep1")
                        .cast(pl.Float64)
                        .mean()
                        .alias("fraction_rep1_positive"),
                        pl.col("read_positive_rep2")
                        .cast(pl.Float64)
                        .mean()
                        .alias("fraction_rep2_positive"),
                        pl.col("both_replicates_positive")
                        .cast(pl.Float64)
                        .mean()
                        .alias("fraction_two_positive"),
                        (pl.col("positive_state_class") == "ONE_POSITIVE")
                        .cast(pl.Float64)
                        .mean()
                        .alias("fraction_one_positive"),
                    ]
                )
            ).to_dicts()[0]
        except Exception as exc:
            eprint(f"[WARNING] Could not compute summary statistics: {exc}")
            stats = {}

    output_columns = collect_schema_names(trajectories)
    report_payload = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "input": input_label,
        "n_input_files": len(input_files),
        "output": str(written_path),
        "output_format": args.output_format,
        "parquet_compression": args.parquet_compression,
        "state_key": key_cols,
        "extra_group_cols": extra_group_cols,
        "partition_by": partition_cols,
        "n_partitions": n_partitions,
        "reference_alpha": float(args.reference_alpha),
        "observable_validation": bool(args.validate_observable),
        "duplicate_state_policy": "error",
        "cross_representation_fallback": False,
        "pseudocount_for_observed_log_frequency": False,
        "output_columns": output_columns,
        "included_subject_clones_output": (
            str(included_path) if included_path is not None else None
        ),
        **stats,
    }
    report_path = outdir / "trajectory_assembly_report.md"
    write_report(report_path, report_payload)

    manifest = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "canonical_step2_source": "2-multirepresentation_clonotype_state_inference.py",
        "canonical_output": str(written_path),
        "state_key": key_cols,
        "state_layers": {
            "model_based_consensus": [
                "f_latent",
                "x_latent",
                "x_latent_sd",
            ],
            "replicate_resolved_observed": [
                "count_rep1",
                "count_rep2",
                "depth_rep1",
                "depth_rep2",
                "f_observed_rep1",
                "f_observed_rep2",
                "x_observed_rep1",
                "x_observed_rep2",
            ],
            "observation_annotations": [
                "read_positive_rep1",
                "read_positive_rep2",
                "p_value",
                "observable_reference",
                "p_detect_state",
                "p_dropout_state",
            ],
            "auxiliary_positive_only": [
                "f_positive_current",
                "x_positive_current",
                "positive_state_class",
            ],
        },
        "reference_alpha": float(args.reference_alpha),
        "observable_definition": "p_value < reference_alpha",
        "latent_consensus_is_ground_truth": False,
        "dynamic_estimands_defined_here": False,
        "primary_downstream_estimands": {
            "forward_drift": "observed replicate-decoupled AB/BA",
            "fluctuations": "cross-replicate covariance of observed displacements on common4",
            "fluctuation_conditioning": "xmid_latent",
        },
        "filtering": False,
        "duplicate_state_policy": "error",
        "cross_representation_fallback": False,
        "observed_log_frequency_pseudocount": False,
    }
    manifest_path = outdir / "trajectory_assembly_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    print(f"[STEP4] Wrote report: {report_path}")
    print(f"[STEP4] Wrote manifest: {manifest_path}")
    print("[STEP4] Completed. No state-quality filtering was applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
