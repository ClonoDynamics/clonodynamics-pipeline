#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2-noise_aware_latent_inference.py

Noise-aware Bayesian inference of latent clonotype abundance from paired
technical-replicate TCR repertoire sequencing data.

Overview
--------
This script implements Step 2 of the ClonoDynamics pipeline.

Its purpose is to infer, for every clonotype observed in a pair of technical
replicates, a latent underlying relative frequency together with posterior
uncertainty and posterior-predictive detectability metrics.

The analysis is explicitly designed to separate the latent clonotype state
from stochastic variation introduced by finite sequencing depth and technical
sampling. Rather than using the observed frequency in either replicate as the
clonotype abundance estimate, the two replicate count observations are treated
as noisy measurements of a common latent frequency.

For each subject and sampling time point, the script:

    1. identifies the two corresponding technical-replicate repertoire files;

    2. reads clonotype amino-acid CDR3 sequences and read counts;

    3. collapses duplicated aaSeqCDR3 entries within each replicate by summing
       their read counts;

    4. merges the two technical replicates by aaSeqCDR3 using the union of all
       clonotypes detected in either replicate;

    5. represents absence from one replicate as an observed count of zero rather
       than removing the clonotype;

    6. computes observed replicate-specific frequencies and descriptive
       replicate-level abundance coordinates;

    7. fits a Negative Binomial observation model with a power-law prior for
       latent clonotype frequency;

    8. evaluates the posterior distribution of latent frequency on a fine
       log-spaced grid;

    9. reports posterior point estimates, credible intervals, posterior
       uncertainty and replicate-specific posterior diagnostics;

   10. propagates posterior-predictive detectability and dropout probabilities;

   11. calculates pair-level fit and replicate-concordance diagnostics;

   12. combines all successfully fitted subject-timepoint pairs into a common
       clonotype-level table;

   13. derives empirical logP-based observability calls using subject-specific,
       pair-specific or global empirical null pools;

   14. generates compact subject/clonotype/time tables of latent frequency and
       observability for subsequent longitudinal transition reconstruction.

The script performs both latent-state inference and preparation of the
downstream-compatible state tables used by later ClonoDynamics steps.


Input
-----
The script expects a directory containing paired technical-replicate repertoire
files for each subject and sampling time point.

By default, repertoire files are identified from their filename using the
pattern

    <subject>_<time>-<replica>.<extension>

where:

    subject
        Integer subject identifier.

    time
        Integer sampling-time identifier.

    replica
        Technical replicate identifier, expected to be 1 or 2.

For example:

    1_1-1.tsv
    1_1-2.tsv
    1_2-1.tsv
    1_2-2.tsv

represents two sampling time points from subject 1, each measured in duplicate.

Only complete replicate pairs are processed. A subject-timepoint is included
when both replica 1 and replica 2 are available.

The default filename pattern recognizes:

    - TSV / delimited text
    - CSV
    - Parquet
    - Feather

The filename convention can be changed through the `--pattern` argument,
provided that the regular expression contains the named groups:

    subject
    time
    replica

Each repertoire file must contain at least:

    aaSeqCDR3
        Amino-acid sequence of the TCR CDR3 clonotype.

    readCount
        Number of sequencing reads assigned to that clonotype.

Additional columns are allowed but are not required by this script.

Within each replicate:

    - readCount is converted to numeric format;
    - non-finite read counts are removed;
    - non-positive read counts are removed;
    - repeated aaSeqCDR3 entries are collapsed by summing readCount.

For text-based files, the default separator is a tab character (`\t`) and can
be changed through `--file-sep`.


Technical-replicate pairing
---------------------------
For a given subject s and sampling time t, let the two replicate files be

    R_(s,t,1)
    R_(s,t,2).

The two repertoires are merged by aaSeqCDR3 using an outer join.

For clonotype i, this produces the observed count pair

    (c_i1, c_i2),

where a clonotype not detected in one replicate is assigned count zero in that
replicate.

Therefore, a clonotype detected in only one replicate is retained as

    (c_i1 > 0, c_i2 = 0)

or

    (c_i1 = 0, c_i2 > 0),

rather than being discarded.

The sequencing depths of the two replicates are

    N_1 = sum_i c_i1

and

    N_2 = sum_i c_i2.

Observed replicate-specific relative frequencies are

    f_i1^obs = c_i1 / N_1

and

    f_i2^obs = c_i2 / N_2.

The script also computes a descriptive observed-frequency coordinate

                        __________________________________
    f_i^geo =          / f_i1^obs * f_i2^obs

when the clonotype is detected in both replicates.

When the clonotype is detected in only one replicate, the detected replicate
frequency is retained as the descriptive `freq_geo` value.

These observed-frequency quantities are diagnostic/descriptive variables.
They are not used as substitutes for the inferred latent clonotype frequency.


Mathematical model
------------------

1. Latent clonotype frequency

For each clonotype i within one technical-replicate pair, let

    f_i > 0

denote the unknown latent relative frequency underlying both technical
measurements.

The two observed read counts are treated as conditionally independent
measurements of this same latent state.


2. Negative Binomial observation model

For replicate r in {1, 2},

    c_ir | f_i, N_r, kappa
        ~ NegativeBinomial(mean = f_i * N_r, size = kappa),

so that the expected observed read count is

    E[c_ir | f_i] = f_i N_r.

The model therefore accounts explicitly for the different sequencing depths
of the two technical replicates.

The dispersion parameter kappa controls the amount of count variability around
the expected count under the Negative Binomial observation model.

The two-replicate likelihood for clonotype i can be written generically as

    L_i(f_i)
        = P(c_i1 | f_i, N_1, kappa)
          P(c_i2 | f_i, N_2, kappa).

The exact numerical implementation and parameterization of the Negative
Binomial likelihood are provided by `noiseK_latent.py`.


3. Power-law prior on latent frequency

Latent clonotype frequencies are assigned a power-law prior

    p(f_i) ∝ f_i^(-gamma),

over the frequency support used by the latent model.

Here:

    gamma
        Power-law exponent controlling the relative prior weight assigned to
        low- versus high-frequency clonotypes.

The use of a heavy-tailed prior reflects the empirically observed broad
clonotype-frequency architecture of TCR repertoires characterized in Step 1
of the ClonoDynamics pipeline.


4. Posterior distribution

For a clonotype with observed counts (c_i1, c_i2), Bayes' rule gives

    p(f_i | c_i1, c_i2, N_1, N_2, gamma, kappa)

        ∝

        P(c_i1 | f_i, N_1, kappa)
        P(c_i2 | f_i, N_2, kappa)
        p(f_i).

Thus, information from both technical replicates is combined into a single
posterior distribution for the latent clonotype frequency.

The posterior is evaluated numerically on a log-spaced frequency grid rather
than represented by a single observed-frequency proxy.

The default grid contains

    500

frequency points and can be changed with

    --grid-size.


5. Log-frequency representation

For downstream dynamical analyses, latent abundance is also represented in
natural-log frequency space:

    x_i = ln(f_i).

The log-frequency representation is particularly useful because clonotype
frequencies span several orders of magnitude.

Posterior summaries are therefore provided in both frequency space and
log-frequency space whenever available.


6. Posterior point estimates

For each clonotype, the latent posterior generated by `noiseK_latent.py`
provides summaries including estimates such as

    f_latent_mean
    f_latent_median
    f_latent_mode

and the corresponding log-frequency quantities

    x_latent_mean
    x_latent_median
    x_latent_mode.

For downstream ClonoDynamics analyses, the posterior median provides a robust
point estimate of the latent abundance state.


7. Posterior credible interval

Posterior uncertainty is summarized using the 95% credible interval

    [f_i,0.025 , f_i,0.975]

in frequency space and, when available,

    [x_i,0.025 , x_i,0.975]

in log-frequency space.

The log-frequency credible-interval width is

    W_i,95
        = x_i,0.975 - x_i,0.025.

This quantity provides a direct clonotype-specific measure of uncertainty in
latent abundance across the large dynamic range of the repertoire.


8. Posterior dispersion and entropy

Additional posterior-uncertainty descriptors include

    x_latent_sd

and

    posterior_entropy.

Posterior entropy summarizes the spread of probability mass across the
discretized latent-frequency posterior.

These quantities are retained as inferential diagnostics and are not used by
this script to remove clonotypes.


Posterior-predictive detectability
----------------------------------
The latent model also returns posterior-predictive probabilities describing
whether a clonotype with the inferred latent state would be detected or missed
under the sequencing observation process.

The principal state-level quantities propagated by this script are

    p_detect_state
        Posterior-predictive detectability of the inferred paired-replicate
        latent state.

    p_dropout_state
        Corresponding posterior-predictive dropout probability.

The script also retains replicate-specific diagnostic quantities, when
provided by `noiseK_latent.py`, including

    p_detect_rep1
    p_detect_rep2
    p_dropout_rep1
    p_dropout_rep2.

State-level detectability is treated as the primary continuous detectability
measure for downstream longitudinal analyses, whereas replicate-specific
probabilities are retained as diagnostic information.

The exact posterior-predictive construction of these quantities is implemented
in `noiseK_latent.py` and is not redefined independently in this wrapper.


Observed replicate states
-------------------------
For descriptive and diagnostic purposes, the script records

    present_rep1
    present_rep2
    present_both
    present_only_one.

These variables describe actual observed read-count status and must be
distinguished from posterior-predictive detectability.

Observed presence is a realized sequencing outcome.

Posterior-predictive detectability is a probabilistic property inferred from
the latent state and observation model.


Replicate-specific latent diagnostics
-------------------------------------
In addition to the joint latent-frequency posterior, `noiseK_latent.py`
provides replicate-specific posterior quantities used to assess consistency
between the two technical measurements.

The script calculates pair-level summaries including:

    replicate_latent_spearman

        Spearman correlation between replicate-specific latent log-frequency
        estimates among clonotypes detected in both replicates.

    replicate_latent_abs_dx_median

        Median absolute difference between replicate-specific latent
        log-frequency estimates.

    replicate_latent_abs_dx_q90

        90th percentile of that absolute difference.

    overlap_jaccard

        Jaccard overlap of the clonotype sets detected in the two technical
        replicates.

    fraction_present_only_one

        Fraction of union clonotypes detected in only one replicate.

These quantities characterize replicate agreement and model behavior but are
not used as automatic exclusion criteria.


Latent count proxy
------------------
The model directly estimates latent relative frequency rather than a latent
integer read count.

For convenience in downstream comparisons with count-based representations,
the script constructs a latent count-scale proxy using the mean sequencing
depth of the pair:

    N_mean = (N_1 + N_2) / 2

and

    c_i^latent
        = f_i^latent,median * N_mean.

This is stored as

    count_latent_median.

A corresponding log-count proxy is

    x_count_latent_median
        = ln(count_latent_median + 1).

These quantities are derived representations of the inferred latent frequency;
they are not independently fitted latent variables.


Pair-level model fitting
------------------------
The latent model is fitted independently for every complete subject-timepoint
technical-replicate pair.

Thus, if subject s has sampling times

    t = 1, ..., T,

the model is fitted separately to

    (s,1), (s,2), ..., (s,T).

The fitted parameters are therefore pair-specific.

The default initial values are

    gamma_init = 1.6
    k_init     = 50.0.

The fitting procedure can be configured through:

    --gamma-init
    --k-init
    --grid-size
    --chunk-size
    --maxiter
    --optimizer-ftol.

No model parameters from one subject-timepoint pair are imposed on another
pair during the primary fitting stage.


Pair-level QC diagnostics
-------------------------
After latent inference, the script creates a pair-level diagnostic table.

This QC stage is explicitly non-filtering.

Possible diagnostic labels include

    PASS_DIAGNOSTIC
    FAIL_FIT
    REVIEW_NUMERICAL
    REVIEW_MULTIPLE_FLAGS
    REVIEW_SINGLE_FLAG.

Diagnostic flags can reflect:

    - failed optimization;
    - non-finite model parameters;
    - non-positive model parameters;
    - unusually large replicate depth imbalance;
    - unusually broad posterior distributions;
    - unusually high posterior entropy;
    - unusually large replicate latent-frequency differences;
    - unusually low replicate overlap;
    - unusually low replicate latent-frequency concordance.

Robust cohort-relative flags are based on median/MAD-style standardized
deviations.

These flags are descriptive only.

The script does NOT exclude successfully fitted pairs on the basis of these
post-fit QC diagnostics.

Only a pair for which latent-model fitting fails lacks a valid downstream
latent output.


Use of all complete replicate pairs
-----------------------------------
This version of Step 2 intentionally processes every complete technical-
replicate pair.

Upstream PASS/FAIL filtering is not applied.

The compatibility arguments

    --use-qc
    --qc-csv

are retained only for legacy command-line compatibility and are ignored.

This design allows posterior uncertainty, detectability and replicate
concordance to be quantified across the complete set of successfully fitted
technical-replicate pairs rather than conditioning downstream inference on an
upstream QC selection.


Empirical logP pooling
----------------------
After pair-level latent inference, all per-clone outputs are combined and an
empirical distribution of the model-derived `logP` statistic is constructed.

The empirical reference pool is controlled by

    --null

with three available choices:

    clone

        Each subject-timepoint replicate pair forms its own empirical pool.
        The pool identifier is the pair_id.

    subject

        All successfully fitted subject-timepoint pairs from the same subject
        are pooled together.

        If a subject has fewer than `--min-pairs-per-subject` available pairs,
        its rows are assigned to the global pool.

    global

        All clonotype observations are pooled into a single empirical
        distribution.

The default is

    --null subject

with

    --min-pairs-per-subject 2.


Empirical logP p-values
-----------------------
Let

    l_i

denote the `logP` value associated with clonotype observation i, and let

    F_pool(l)

denote the empirical cumulative distribution function of logP within the
assigned pool.

The script calculates

    p_emp_low,i
        = F_pool(l_i)

and

    p_emp_high,i
        = 1 - F_pool(l_i).

The tail used for the observability call is selected with

    --tail low

or

    --tail high.

For the default low-tail analysis,

    p_value_i = p_emp_low,i.

For high-tail analysis,

    p_value_i = p_emp_high,i.


Empirical observability
-----------------------
A clonotype-timepoint observation is classified as empirically observable when

    p_value_i < alpha,

where alpha is controlled by

    --alpha

and defaults to

    alpha = 0.05.

The resulting Boolean field is

    observable.

This empirical observability classification is distinct from the continuous
posterior-predictive detectability quantities

    p_detect_state
    p_dropout_state.

The first is a thresholded classification derived from the empirical
distribution of `logP`.

The second represents continuous posterior-predictive information returned by
the latent observation model.

Both are retained because they answer different questions.


Alpha sensitivity and aggregate-only mode
-----------------------------------------
Latent-model fitting itself does not depend on the empirical observability
threshold alpha.

Therefore, once the computationally intensive per-pair latent posterior files
have been generated, alternative alpha thresholds can be evaluated without
refitting the latent model.

The

    --aggregate-only

mode reads an existing per-clone latent cache and recomputes:

    - empirical logP p-values;
    - empirical cutoffs;
    - observable classifications;
    - downstream trajectory tables.

This is intended for sensitivity analyses such as

    alpha = 0.05
    alpha = 0.02
    alpha = 0.01

while keeping the same underlying latent-frequency inference.


Output
------
The script creates a structured set of pair-level, clonotype-level and
trajectory-level outputs.

For a results directory

    <results-dir>/

the principal outputs are:


1. per_clone_latent_logP/

    One file per successfully fitted subject-timepoint replicate pair:

        latent_per_clone_logP_<subject>_<time>.parquet

    or the equivalent Feather/CSV file according to
    `--intermediate-format`.

    Each row represents one clonotype in the union of the two technical
    replicates and may contain:

        aaSeqCDR3

        count_rep1
        count_rep2

        freq_rep1
        freq_rep2
        freq_geo

        x_obs_rep1
        x_obs_rep2
        x_obs_geo

        present_rep1
        present_rep2
        present_both
        present_only_one

        latent posterior summaries
        posterior uncertainty metrics
        replicate-specific latent diagnostics

        p_detect_state
        p_dropout_state

        replicate-specific detectability/dropout diagnostics

        depth_rep1
        depth_rep2
        depth_mean

        count_latent_median
        x_count_latent_median

        subject
        time
        pair_id

        fitted pair-level model parameters.


2. posterior NPZ files

    For each fitted pair, `noiseK_latent.py` may additionally write

        latent_per_clone_logP_<subject>_<time>.posterior.npz

    containing the numerical posterior representation required for detailed
    posterior diagnostics or downstream posterior propagation.


3. latent_params_by_pair.csv

    One row per attempted technical-replicate pair.

    Contains:

        pair identifiers;
        input filenames;
        sequencing depths;
        fitted model parameters;
        log-likelihood;
        replicate overlap;
        posterior uncertainty summaries;
        detectability summaries;
        replicate concordance metrics;
        fit success/failure information.


4. latent_pair_qc.csv

    Pair-level non-filtering diagnostic table.

    This file reports numerical and cohort-relative warning flags but does not
    remove successfully fitted pairs from downstream analysis.


5. latent_params_aggregated.csv

    Median fitted model parameters aggregated according to the selected null
    structure.

    Under the default

        --null subject

    parameters are summarized by subject.

    For other null choices, a global summary is generated.


6. latent_logP_cutoffs_<null>.csv

    Empirical logP reference thresholds for each pool.

    Contains fields including

        pool_id
        alpha
        tail
        logP_cutoff_low_tail
        logP_cutoff_high_tail
        n_rows_pool.


7. per_clone_latent_<null>.<format>

    Main combined clonotype-level output.

    It contains all per-pair latent results together with

        pool_id
        p_emp_low
        p_emp_high
        p_value
        observable.

    This is the complete downstream-compatible clonotype-level table.


8. latent_trajectory_observability_<null>.csv

    Compact table with one row per

        subject × aaSeqCDR3 × time

    summarizing empirical observability.

    Principal fields include

        subject
        aaSeqCDR3
        time
        observable
        p_emp_low
        p_value
        n_rows.


9. latent_trajectory_frequency_<null>.csv

    Compact latent-state table with one row per

        subject × aaSeqCDR3 × time.

    Depending on the fields provided by `noiseK_latent.py`, it propagates
    quantities including

        count_rep1
        count_rep2
        depth_rep1
        depth_rep2

        freq_rep1
        freq_rep2
        freq_geo

        f_latent_mean
        f_latent_median
        f_latent_mode
        f_latent_q025
        f_latent_q975

        x_latent_mean
        x_latent_median
        x_latent_mode
        x_latent_q025
        x_latent_q975
        x_latent_sd

        p_detect_state
        p_dropout_state

        p_detect_rep1
        p_detect_rep2
        p_dropout_rep1
        p_dropout_rep2

        logP
        p_value
        observable.

    This table provides the principal latent clonotype-state representation
    used for subsequent longitudinal trajectory and transition construction.


Intermediate file formats
-------------------------
Per-clone tables can be written as

    parquet
    feather
    csv

through

    --intermediate-format.

The default is

    parquet.

Parquet compression defaults to

    snappy

and can be changed with

    --parquet-compression.

When a non-CSV intermediate format is selected,

    --write-csv-compat

can optionally generate CSV compatibility copies of the pair-level per-clone
tables.


Parallel execution
------------------
Pair-level model fitting can be parallelized across independent
subject-timepoint replicate pairs using

    --n-jobs.

Because each replicate pair is fitted independently, parallelization does not
change the statistical model.

The script also limits nested numerical-library threading to reduce CPU
oversubscription during multiprocessing.


Software dependency
-------------------
This script requires

    noiseK_latent.py

and specifically imports

    fit_noiseK_latent_powerlaw.

`noiseK_latent.py` must either:

    - be located in the same directory as this script; or
    - be available through PYTHONPATH.

The detailed Negative Binomial likelihood implementation, latent posterior
construction, model optimization and posterior-predictive detectability
calculations are defined in that module.


Typical usage
-------------
Standard fit + aggregate analysis:

    python3 2-noise_aware_latent_inference.py \
        --data-dir ./longitudinali \
        --results-dir ./results/2-noise_aware_latent_inference/p_05 \
        --alpha 0.05 \
        --null subject \
        --tail low \
        --n-jobs 4 \
        --intermediate-format parquet \
        --parquet-compression snappy


Example with automatic macOS-friendly parallelization:

    python3 2-noise_aware_latent_inference.py \
        --data-dir ./longitudinali \
        --results-dir ./results/2-noise_aware_latent_inference/p_05 \
        --alpha 0.05 \
        --null subject \
        --n-jobs 0 \
        --n-jobs-profile balanced


Re-aggregation at a different alpha without refitting:

    python3 2-noise_aware_latent_inference.py \
        --aggregate-only \
        --fit-cache-dir ./results/2-noise_aware_latent_inference/p_05 \
        --aggregate-results-dir ./results/2-noise_aware_latent_inference/p_02 \
        --alpha 0.02 \
        --null subject \
        --tail low


Conceptual role in ClonoDynamics
--------------------------------
Step 1 of ClonoDynamics characterizes the observed repertoire architecture,
including the heavy-tailed distribution of clonotype frequencies.

This Step 2 then converts paired noisy repertoire measurements into
probabilistic latent clonotype states:

    observed replicate counts
                |
                v
    Negative Binomial observation model
                +
        power-law frequency prior
                |
                v
       posterior p(f_i | data)
                |
        +-------+-------+
        |               |
        v               v
 latent abundance    uncertainty
        |               |
        +-------+-------+
                |
                v
    posterior-predictive detectability
                |
                v
 empirical observability annotation
                |
                v
 subject × clonotype × time latent states
                |
                v
 downstream longitudinal transition analysis

The central output is therefore not a denoised observed frequency, but a
posterior representation of the latent clonotype abundance state together
with its inferential uncertainty and detectability.

This distinction is essential for subsequent ClonoDynamics analyses, in which
finite-time clonotype transitions are constructed from latent states rather
than directly from raw replicate counts or single-replicate observed
frequencies.
"""

from __future__ import annotations

import argparse
import gc
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

# Avoid nested BLAS/NumPy oversubscription when using multiprocessing.
for _var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_var, "1")

import numpy as np
import pandas as pd


PATTERN_DEFAULT = (
    r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])"
    r"(?:\.(?:tsv|csv|parquet|feather))?$"
)

PARQUET_EXTS = {".parquet", ".pq"}
FEATHER_EXTS = {".feather", ".arrow"}
CSV_EXTS = {".csv", ".tsv", ".txt"}


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(p: Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_macos_n_jobs(n_jobs: int, profile: str = "balanced") -> int:
    if int(n_jobs) > 0:
        return int(n_jobs)
    cpu = os.cpu_count() or 4
    profile = str(profile or "balanced").lower().strip()
    if profile == "memory":
        return 1
    if profile == "speed":
        return max(1, min(6, cpu - 2))
    return max(1, min(4, cpu - 2))


def configure_arrow_runtime(n_jobs: int) -> None:
    try:
        import pyarrow as pa
        pa.set_cpu_count(max(1, int(n_jobs)))
        try:
            pa.set_io_thread_count(max(1, min(4, int(n_jobs))))
        except Exception:
            pass
    except Exception:
        pass


def format_ext(fmt: str) -> str:
    fmt = str(fmt).lower().strip()
    if fmt == "parquet":
        return ".parquet"
    if fmt == "feather":
        return ".feather"
    if fmt == "csv":
        return ".csv"
    raise ValueError(f"Unknown format: {fmt}")


def write_table(df: pd.DataFrame, path: Path, fmt: str = "parquet", compression: str = "snappy", index: bool = False) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    fmt = str(fmt).lower().strip()
    if fmt == "parquet":
        out = path.with_suffix(".parquet")
        df.to_parquet(out, compression=compression, index=index)
        return out
    if fmt == "feather":
        out = path.with_suffix(".feather")
        if index:
            df = df.reset_index()
        df.to_feather(out)
        return out
    if fmt == "csv":
        out = path.with_suffix(".csv")
        df.to_csv(out, index=index)
        return out
    raise ValueError(f"Unknown format: {fmt}")


def read_table(fp: Path, file_sep: str = "\t", columns: Optional[Sequence[str]] = None) -> pd.DataFrame:
    fp = Path(fp)
    suffix = fp.suffix.lower()
    if suffix in PARQUET_EXTS:
        return pd.read_parquet(fp, columns=list(columns) if columns else None)
    if suffix in FEATHER_EXTS:
        df = pd.read_feather(fp)
    elif suffix == ".csv":
        df = pd.read_csv(fp)
    else:
        df = pd.read_csv(fp, sep=file_sep)
    if columns:
        keep = [c for c in columns if c in df.columns]
        return df[keep].copy()
    return df


def list_perclone_files(perclone_dir: Path) -> List[Path]:
    perclone_dir = Path(perclone_dir)
    for pattern in ("*.parquet", "*.feather", "*.csv"):
        files = sorted(perclone_dir.glob(pattern))
        if files:
            return files
    return []


def empirical_cdf(sorted_arr: np.ndarray, x: float) -> float:
    if sorted_arr.size == 0:
        return float("nan")
    return float(np.searchsorted(sorted_arr, x, side="right") / sorted_arr.size)


# =============================================================================
# Raw repertoire input and pair indexing
# =============================================================================


def read_rep_file(fp: Path, sep: str = "\t") -> pd.DataFrame:
    df = read_table(fp, file_sep=sep)
    required = {"aaSeqCDR3", "readCount"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{fp.name}: missing columns {sorted(missing)}; expected {sorted(required)}")
    out = df[["aaSeqCDR3", "readCount"]].copy()
    out["readCount"] = pd.to_numeric(out["readCount"], errors="coerce")
    out = out[np.isfinite(out["readCount"])].copy()
    out = out[out["readCount"] > 0].copy()
    # Collapse duplicated aaSeqCDR3 rows if present.
    out = out.groupby("aaSeqCDR3", as_index=False)["readCount"].sum()
    out["readCount"] = out["readCount"].astype(np.int64)
    return out


def merge_two_reps(fp1: Path, fp2: Path, sep: str = "\t", key: str = "aaSeqCDR3"):
    d1 = read_rep_file(fp1, sep=sep).rename(columns={"readCount": "count_rep1"})
    d2 = read_rep_file(fp2, sep=sep).rename(columns={"readCount": "count_rep2"})
    merged = d1.merge(d2, on=key, how="outer")
    merged["count_rep1"] = merged["count_rep1"].fillna(0).astype(np.int64)
    merged["count_rep2"] = merged["count_rep2"].fillna(0).astype(np.int64)
    depth_rep1 = int(merged["count_rep1"].sum())
    depth_rep2 = int(merged["count_rep2"].sum())
    if depth_rep1 <= 0 or depth_rep2 <= 0:
        raise ValueError(f"Non-positive sequencing depth: rep1={depth_rep1}, rep2={depth_rep2}")
    depths = [depth_rep1, depth_rep2]
    count_cols = ["count_rep1", "count_rep2"]
    return merged, depths, count_cols, depth_rep1, depth_rep2


def index_pairs(data_dir: Path, pattern: re.Pattern) -> List[Tuple[int, int, Path, Path]]:
    index: Dict[Tuple[int, int, int], Path] = {}
    for fp in Path(data_dir).iterdir():
        if not fp.is_file():
            continue
        m = pattern.match(fp.name)
        if not m:
            continue
        subject = int(m.group("subject"))
        time = int(m.group("time"))
        replica = int(m.group("replica"))
        index[(subject, time, replica)] = fp
    pairs: List[Tuple[int, int, Path, Path]] = []
    for subject, time in sorted({(s, t) for (s, t, _r) in index.keys()}):
        fp1 = index.get((subject, time, 1))
        fp2 = index.get((subject, time, 2))
        if fp1 and fp2:
            pairs.append((subject, time, fp1, fp2))
    return pairs


# =============================================================================
# Optional QC filtering
# =============================================================================


def parse_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    norm = s.astype(str).str.strip().str.lower()
    truthy = {"true", "1", "yes", "y", "pass", "passed"}
    return norm.isin(truthy)


def load_qc_pass_pairs(qc_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(qc_csv)
    required = {"subject", "time", "QC_pass"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"QC CSV missing required columns: {sorted(missing)}")
    out = df[["subject", "time", "QC_pass"]].copy()
    out["subject"] = pd.to_numeric(out["subject"], errors="raise").astype(int)
    out["time"] = pd.to_numeric(out["time"], errors="raise").astype(int)
    out["QC_pass"] = parse_bool_series(out["QC_pass"])
    return out.drop_duplicates(subset=["subject", "time"], keep="first")


def build_allowed_pair_set(qc_df: pd.DataFrame) -> Set[Tuple[int, int]]:
    keep = qc_df.loc[qc_df["QC_pass"] == True, ["subject", "time"]].copy()
    return set(map(tuple, keep.to_records(index=False)))


# =============================================================================
# Latent model fitting
# =============================================================================


def _fit_one_pair_worker(job: dict) -> dict:
    subject = int(job["subject"])
    time = int(job["time"])
    fp1 = Path(job["fp1"])
    fp2 = Path(job["fp2"])
    pair_id = f"{subject}_{time}"
    allowed_pairs_was_used = bool(job.get("allowed_pairs_was_used", False))

    try:
        from noiseK_latent import fit_noiseK_latent_powerlaw

        merged, depths, count_cols, depth_rep1, depth_rep2 = merge_two_reps(
            fp1, fp2, sep=job["file_sep"]
        )

        merged["freq_rep1"] = merged["count_rep1"] / depth_rep1
        merged["freq_rep2"] = merged["count_rep2"] / depth_rep2

        f1 = merged["freq_rep1"].to_numpy(dtype=float)
        f2 = merged["freq_rep2"].to_numpy(dtype=float)
        c1 = merged["count_rep1"].to_numpy(dtype=np.int64)
        c2 = merged["count_rep2"].to_numpy(dtype=np.int64)

        freq_geo = np.zeros(len(merged), dtype=float)
        both = (c1 > 0) & (c2 > 0)
        only1 = (c1 > 0) & (c2 == 0)
        only2 = (c2 > 0) & (c1 == 0)
        freq_geo[both] = np.sqrt(f1[both] * f2[both])
        freq_geo[only1] = f1[only1]
        freq_geo[only2] = f2[only2]
        merged["freq_geo"] = freq_geo
        merged["x_obs_rep1"] = np.where(c1 > 0, np.log(np.maximum(f1, 1e-300)), np.nan)
        merged["x_obs_rep2"] = np.where(c2 > 0, np.log(np.maximum(f2, 1e-300)), np.nan)
        merged["x_obs_geo"] = np.log(np.maximum(freq_geo, 1e-300))
        merged["present_rep1"] = c1 > 0
        merged["present_rep2"] = c2 > 0
        merged["present_both"] = both
        merged["present_only_one"] = only1 | only2

        posterior_npz_path = Path(job["perclone_dir"]) / f"latent_per_clone_logP_{pair_id}.posterior.npz"

        fit = fit_noiseK_latent_powerlaw(
            merged=merged,
            depths=depths,
            count_cols=count_cols,
            gamma_init=float(job["gamma_init"]),
            k_init=float(job["k_init"]),
            grid_size=int(job["grid_size"]),
            chunk_size=int(job["chunk_size"]),
            maxiter=job.get("maxiter"),
            optimizer_ftol=job.get("optimizer_ftol"),
            posterior_npz_path=str(posterior_npz_path),
        )

        per_clone = fit.per_clone_logP.copy()
        per_clone["depth_rep1"] = depth_rep1
        per_clone["depth_rep2"] = depth_rep2

        # The updated noiseK_latent model provides one posterior-predictive
        # detectability for the paired-replicate state, in addition to the
        # replicate-specific diagnostic probabilities.
        state_detect_cols = {"p_detect_state", "p_dropout_state"}
        missing_state_detect = state_detect_cols - set(per_clone.columns)
        if missing_state_detect:
            raise RuntimeError(
                "noiseK_latent output is missing state-level detectability columns: "
                f"{sorted(missing_state_detect)}. Use the updated noiseK_latent.py."
            )

        # Latent count proxy for downstream analyses.
        # The latent model estimates a frequency f; multiplying by the mean
        # sequencing depth gives the expected read-count scale for the pair.
        depth_mean = 0.5 * (float(depth_rep1) + float(depth_rep2))
        per_clone["depth_mean"] = depth_mean
        if "f_latent_median" in per_clone.columns:
            f_lat = pd.to_numeric(per_clone["f_latent_median"], errors="coerce")
            per_clone["count_latent_median"] = f_lat * depth_mean
            per_clone["x_count_latent_median"] = np.log(np.maximum(per_clone["count_latent_median"], 0.0) + 1.0)

        per_clone["subject"] = subject
        per_clone["time"] = time
        per_clone["pair_id"] = pair_id

        params = getattr(fit, "params", {}) if getattr(fit, "params", None) is not None else {}
        if isinstance(params, dict):
            for k, v in params.items():
                # Attach model params to each row for direct downstream use.
                if k not in per_clone.columns:
                    per_clone[k] = v
                else:
                    per_clone[f"param_{k}"] = v

        stem = f"latent_per_clone_logP_{pair_id}"
        out_fp = write_table(
            per_clone,
            Path(job["perclone_dir"]) / stem,
            fmt=job["intermediate_format"],
            compression=job["parquet_compression"],
            index=False,
        )
        if bool(job.get("write_csv_compat", False)) and job["intermediate_format"] != "csv":
            per_clone.to_csv((Path(job["perclone_dir"]) / stem).with_suffix(".csv"), index=False)

        n_clones_rep1 = int((merged["count_rep1"] > 0).sum())
        n_clones_rep2 = int((merged["count_rep2"] > 0).sum())
        inter = int(((merged["count_rep1"] > 0) & (merged["count_rep2"] > 0)).sum())
        union = int(((merged["count_rep1"] > 0) | (merged["count_rep2"] > 0)).sum())
        overlap_jaccard = float(inter / union) if union > 0 else float("nan")
        depth_ratio = float(max(depth_rep1, depth_rep2) / max(1, min(depth_rep1, depth_rep2)))

        # Pair-level latent QC diagnostics. These summarize the fit but do not
        # determine whether the pair is retained.
        posterior_x_width = (
            pd.to_numeric(per_clone.get("x_latent_q975"), errors="coerce")
            - pd.to_numeric(per_clone.get("x_latent_q025"), errors="coerce")
        )
        shared_detected = (per_clone["count_rep1"] > 0) & (per_clone["count_rep2"] > 0)
        rep_corr, n_corr = _safe_spearman(
            per_clone.loc[shared_detected, "x_rep1_latent_median"],
            per_clone.loc[shared_detected, "x_rep2_latent_median"],
            min_n=20,
        )
        abs_dx = pd.to_numeric(per_clone.get("dx_rep_latent_median"), errors="coerce").abs()
        fraction_present_only_one = float(per_clone["present_only_one"].mean()) if len(per_clone) else float("nan")
        fraction_joint_posterior_wide = float((posterior_x_width > 4.0).mean()) if len(per_clone) else float("nan")
        p_detect_state = pd.to_numeric(per_clone.get("p_detect_state"), errors="coerce")
        p_dropout_state = pd.to_numeric(per_clone.get("p_dropout_state"), errors="coerce")

        row = {
            "subject": subject,
            "time": time,
            "pair_id": pair_id,
            "file_rep1": fp1.name,
            "file_rep2": fp2.name,
            "depth_rep1": depth_rep1,
            "depth_rep2": depth_rep2,
            "depth_ratio": depth_ratio,
            "n_clones_rep1": n_clones_rep1,
            "n_clones_rep2": n_clones_rep2,
            "n_clones_union": union,
            "n_clones_intersection": inter,
            "overlap_jaccard": overlap_jaccard,
            "n_clones_qc": int(len(per_clone)),
            "n_shared_detected_qc": int(n_corr),
            "posterior_x_width_median": _safe_numeric_summary(posterior_x_width, "median"),
            "posterior_x_width_q90": _safe_numeric_summary(posterior_x_width, "q90"),
            "posterior_x_sd_median": _safe_numeric_summary(per_clone.get("x_latent_sd"), "median"),
            "posterior_entropy_median": _safe_numeric_summary(per_clone.get("posterior_entropy"), "median"),
            "p_detect_state_median": _safe_numeric_summary(p_detect_state, "median"),
            "p_detect_state_q10": _safe_numeric_summary(p_detect_state, "q10"),
            "p_dropout_state_median": _safe_numeric_summary(p_dropout_state, "median"),
            "replicate_latent_spearman": rep_corr,
            "replicate_latent_abs_dx_median": _safe_numeric_summary(abs_dx, "median"),
            "replicate_latent_abs_dx_q90": _safe_numeric_summary(abs_dx, "q90"),
            "fraction_present_only_one": fraction_present_only_one,
            "fraction_joint_posterior_wide": fraction_joint_posterior_wide,
            "success": bool(getattr(fit, "success", True)),
            "message": str(getattr(fit, "message", "")),
            "log_likelihood": float(getattr(fit, "log_likelihood", np.nan)),
            "QC_pass": True if allowed_pairs_was_used else np.nan,
            "QC_used_for_filtering": allowed_pairs_was_used,
            "perclone_file": str(out_fp),
        }
        if isinstance(params, dict):
            for k, v in params.items():
                row[k] = v

        del merged, per_clone, fit, f1, f2, c1, c2, freq_geo
        gc.collect()
        return row

    except Exception as e:
        return {
            "subject": subject,
            "time": time,
            "pair_id": pair_id,
            "file_rep1": fp1.name,
            "file_rep2": fp2.name,
            "success": False,
            "message": f"ERROR: {type(e).__name__}: {e}",
            "QC_pass": True if allowed_pairs_was_used else np.nan,
            "QC_used_for_filtering": allowed_pairs_was_used,
        }



def _safe_numeric_summary(series: pd.Series, func: str) -> float:
    x = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if x.empty:
        return float("nan")
    if func == "median":
        return float(x.median())
    if func == "mean":
        return float(x.mean())
    if func == "q10":
        return float(x.quantile(0.10))
    if func == "q90":
        return float(x.quantile(0.90))
    raise ValueError(f"Unknown summary function: {func}")


def _safe_spearman(x: pd.Series, y: pd.Series, min_n: int = 20) -> Tuple[float, int]:
    a = pd.to_numeric(x, errors="coerce")
    b = pd.to_numeric(y, errors="coerce")
    keep = np.isfinite(a) & np.isfinite(b)
    n = int(keep.sum())
    if n < int(min_n):
        return float("nan"), n
    return float(a[keep].corr(b[keep], method="spearman")), n


def _robust_z(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    med = x.median(skipna=True)
    mad = (x - med).abs().median(skipna=True)
    if not np.isfinite(mad) or mad <= 0:
        return pd.Series(np.nan, index=x.index, dtype=float)
    return (x - med) / mad


def build_latent_pair_qc(params_csv: Path, results_dir: Path, robust_cutoff: float = 3.5) -> Path:
    """Create a non-filtering post-fit QC table from pair-level latent diagnostics."""
    df = pd.read_csv(params_csv)
    if df.empty:
        raise ValueError(f"Cannot build latent pair QC from empty table: {params_csv}")

    success = parse_bool_series(df.get("success", pd.Series(False, index=df.index)))
    df["fit_failure_flag"] = ~success

    # Parameter validity / numerical validity.
    gamma = pd.to_numeric(df.get("gamma", np.nan), errors="coerce")
    kappa = pd.to_numeric(df.get("k", np.nan), errors="coerce")
    ll = pd.to_numeric(df.get("log_likelihood", np.nan), errors="coerce")
    df["nonfinite_parameter_flag"] = ~(np.isfinite(gamma) & np.isfinite(kappa) & np.isfinite(ll))
    df["nonpositive_parameter_flag"] = (gamma <= 0) | (kappa <= 0)

    # Robust outlier flags are descriptive and cohort-relative; they do not filter data.
    high_metrics = [
        "depth_ratio",
        "posterior_x_width_median",
        "posterior_x_sd_median",
        "posterior_entropy_median",
        "replicate_latent_abs_dx_median",
        "replicate_latent_abs_dx_q90",
    ]
    low_metrics = ["overlap_jaccard", "replicate_latent_spearman"]

    flag_cols = []
    for col in high_metrics:
        if col in df.columns:
            zcol = f"robust_z_{col}"
            fcol = f"outlier_high_{col}"
            df[zcol] = _robust_z(df[col])
            df[fcol] = df[zcol] > float(robust_cutoff)
            flag_cols.append(fcol)
    for col in low_metrics:
        if col in df.columns:
            zcol = f"robust_z_{col}"
            fcol = f"outlier_low_{col}"
            df[zcol] = _robust_z(df[col])
            df[fcol] = df[zcol] < -float(robust_cutoff)
            flag_cols.append(fcol)

    if flag_cols:
        df["n_robust_outlier_flags"] = df[flag_cols].fillna(False).astype(bool).sum(axis=1)
    else:
        df["n_robust_outlier_flags"] = 0

    df["latent_qc_status"] = np.select(
        [
            df["fit_failure_flag"],
            df["nonfinite_parameter_flag"] | df["nonpositive_parameter_flag"],
            df["n_robust_outlier_flags"] >= 2,
            df["n_robust_outlier_flags"] == 1,
        ],
        ["FAIL_FIT", "REVIEW_NUMERICAL", "REVIEW_MULTIPLE_FLAGS", "REVIEW_SINGLE_FLAG"],
        default="PASS_DIAGNOSTIC",
    )
    df["included_in_downstream"] = success
    df["qc_filter_applied"] = False
    df["qc_interpretation"] = (
        "Post-fit diagnostic only; no pair was excluded by latent_pair_qc. "
        "Only pairs whose model fit failed lack downstream latent output."
    )

    preferred = [
        "subject", "time", "pair_id", "file_rep1", "file_rep2",
        "success", "message", "latent_qc_status", "included_in_downstream",
        "gamma", "k", "log_likelihood", "fmin", "grid_size", "chunk_size",
        "depth_rep1", "depth_rep2", "depth_ratio",
        "n_clones_rep1", "n_clones_rep2", "n_clones_union", "n_clones_intersection", "overlap_jaccard",
        "n_clones_qc", "n_shared_detected_qc",
        "posterior_x_width_median", "posterior_x_width_q90",
        "posterior_x_sd_median", "posterior_entropy_median",
        "p_detect_state_median", "p_detect_state_q10", "p_dropout_state_median",
        "replicate_latent_spearman", "replicate_latent_abs_dx_median", "replicate_latent_abs_dx_q90",
        "fraction_present_only_one", "fraction_joint_posterior_wide",
        "n_robust_outlier_flags", "fit_failure_flag", "nonfinite_parameter_flag",
        "nonpositive_parameter_flag", "qc_filter_applied", "qc_interpretation",
    ]
    ordered = [c for c in preferred if c in df.columns]
    remainder = [c for c in df.columns if c not in ordered]
    df = df[ordered + remainder]

    out = Path(results_dir) / "latent_pair_qc.csv"
    df.to_csv(out, index=False)
    print("Saved non-filtering latent pair QC:", out)
    return out

def run_latent_fit(
    data_dir: Path,
    results_dir: Path,
    perclone_dir: Path,
    file_sep: str,
    pattern_str: str,
    gamma_init: float,
    k_init: float,
    grid_size: int,
    chunk_size: int,
    allowed_pairs: Optional[Set[Tuple[int, int]]],
    intermediate_format: str,
    parquet_compression: str,
    write_csv_compat: bool,
    n_jobs: int,
    maxiter: Optional[int] = None,
    optimizer_ftol: Optional[float] = None,
) -> Path:
    try:
        from noiseK_latent import fit_noiseK_latent_powerlaw  # noqa: F401
    except Exception as e:
        raise ImportError(
            "Unable to import noiseK_latent.fit_noiseK_latent_powerlaw. "
            "Put noiseK_latent.py in the same folder as this script or set PYTHONPATH."
        ) from e

    ensure_dir(results_dir)
    ensure_dir(perclone_dir)

    pattern = re.compile(pattern_str)
    all_pairs = index_pairs(data_dir, pattern)
    # Step 2 intentionally uses every complete replicate pair. Upstream QC is
    # not consulted and no pair is removed before latent inference.
    selected_pairs = all_pairs
    print(f"Found {len(all_pairs)} complete replicate pairs in {data_dir}")
    print("[all-pairs] Upstream QC filtering is disabled; keeping every complete pair.")

    if not all_pairs:
        examples = sorted([fp.name for fp in Path(data_dir).iterdir() if fp.is_file()])[:10]
        raise ValueError(
            "No complete replicate pairs found. Check --pattern.\n"
            f"Pattern: {pattern_str}\nFirst files: {examples}"
        )
    if not selected_pairs:
        raise ValueError("No complete replicate pairs available for latent inference.")

    n_jobs = max(1, min(int(n_jobs), len(selected_pairs)))
    configure_arrow_runtime(n_jobs)

    jobs = []
    for subject, time, fp1, fp2 in selected_pairs:
        jobs.append({
            "subject": int(subject),
            "time": int(time),
            "fp1": str(fp1),
            "fp2": str(fp2),
            "perclone_dir": str(perclone_dir),
            "file_sep": file_sep,
            "gamma_init": float(gamma_init),
            "k_init": float(k_init),
            "grid_size": int(grid_size),
            "chunk_size": int(chunk_size),
            "maxiter": maxiter,
            "optimizer_ftol": optimizer_ftol,
            "intermediate_format": intermediate_format,
            "parquet_compression": parquet_compression,
            "write_csv_compat": bool(write_csv_compat),
            "allowed_pairs_was_used": False,
        })

    print(f"[latent-fit] n_pairs={len(jobs)} n_jobs={n_jobs} grid_size={grid_size} chunk_size={chunk_size}")
    rows: List[dict] = []

    if n_jobs == 1:
        for i, job in enumerate(jobs, start=1):
            pair_id = f"{job['subject']}_{job['time']}"
            print(f"\n--- [{i}/{len(jobs)}] Processing {pair_id}: {Path(job['fp1']).name} + {Path(job['fp2']).name}")
            row = _fit_one_pair_worker(job)
            if row.get("success", False):
                print(f"Saved latent per-clone file: {row.get('perclone_file', '')}")
            else:
                eprint(f"FAILED {pair_id}: {row.get('message', '')}")
            rows.append(row)
            gc.collect()
    else:
        print("[latent-fit] Parallel mode enabled.")
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            future_to_pair = {ex.submit(_fit_one_pair_worker, job): f"{job['subject']}_{job['time']}" for job in jobs}
            for completed, fut in enumerate(as_completed(future_to_pair), start=1):
                pair_id = future_to_pair[fut]
                try:
                    row = fut.result()
                except Exception as e:
                    row = {"pair_id": pair_id, "success": False, "message": f"ERROR: {type(e).__name__}: {e}"}
                if row.get("success", False):
                    print(f"[latent-fit] [{completed}/{len(jobs)}] DONE {pair_id}")
                else:
                    eprint(f"[latent-fit] [{completed}/{len(jobs)}] FAILED {pair_id}: {row.get('message', '')}")
                rows.append(row)
                gc.collect()

    rows = sorted(rows, key=lambda r: (int(r.get("subject", 10**12)), int(r.get("time", 10**12))))
    out_table = Path(results_dir) / "latent_params_by_pair.csv"
    pd.DataFrame(rows).to_csv(out_table, index=False)
    n_ok = int(sum(bool(r.get("success", False)) for r in rows))
    print(f"\nSaved parameter table: {out_table}")
    print(f"Saved per-clone latent files in: {perclone_dir}")
    print(f"[latent-fit] completed: success={n_ok}, failed={len(rows) - n_ok}")
    if n_ok == 0:
        raise RuntimeError(f"Latent noise fitting failed for all pairs. See {out_table}")
    return out_table


# =============================================================================
# Aggregation, empirical p-values and downstream-compatible outputs
# =============================================================================


def infer_clone_and_logp_cols(df: pd.DataFrame) -> Tuple[str, str]:
    clone_candidates = ["aaSeqCDR3", "cloneId", "cloneID", "clonotype", "cdr3"]
    logp_candidates = ["logP", "logp", "log_prob", "log_probability", "logP_noise", "logPnull"]
    clone_col = next((c for c in clone_candidates if c in df.columns), None)
    logp_col = next((c for c in logp_candidates if c in df.columns), None)
    if clone_col is None:
        raise ValueError(f"Cannot infer clone column. Columns: {list(df.columns)}")
    if logp_col is None:
        raise ValueError(f"Cannot infer logP column. Columns: {list(df.columns)}")
    return clone_col, logp_col


def load_perclone_latent(perclone_dir: Path, file_sep: str = "\t") -> pd.DataFrame:
    files = list_perclone_files(perclone_dir)
    if not files:
        raise FileNotFoundError(f"No per-clone latent files found in {perclone_dir}")
    rows = []
    print(f"Reading {len(files)} per-clone latent files from {perclone_dir}")
    for i, fp in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] {fp.name}")
        df = read_table(fp, file_sep=file_sep)
        if {"subject", "time", "pair_id"}.issubset(df.columns):
            rows.append(df)
            continue
        # Fallback from filename latent_per_clone_logP_<subject>_<time>
        m = re.search(r"(?P<subject>\d+)_(?P<time>\d+)", fp.stem)
        if not m:
            eprint(f"[skip] could not parse subject/time from {fp.name}")
            continue
        df["subject"] = int(m.group("subject"))
        df["time"] = int(m.group("time"))
        df["pair_id"] = f"{int(m.group('subject'))}_{int(m.group('time'))}"
        rows.append(df)
    if not rows:
        raise ValueError(f"No usable per-clone latent files found in {perclone_dir}")
    return pd.concat(rows, ignore_index=True)


def assign_pools(df: pd.DataFrame, choice: str, min_pairs_per_subject: int) -> pd.DataFrame:
    df = df.copy()
    choice = str(choice).lower().strip()
    if choice == "clone":
        df["pool_id"] = df["pair_id"].astype(str)
    elif choice == "global":
        df["pool_id"] = "GLOBAL"
    elif choice == "subject":
        pair_counts = (
            df[["subject", "pair_id"]]
            .drop_duplicates()
            .groupby("subject")
            .size()
            .rename("n_pairs_subject")
            .to_dict()
        )
        df["n_pairs_subject"] = df["subject"].map(pair_counts).fillna(0).astype(int)
        df["pool_id"] = np.where(
            df["n_pairs_subject"] >= int(min_pairs_per_subject),
            df["subject"].astype(str),
            "GLOBAL",
        )
    else:
        raise ValueError("--null must be one of: clone, subject, global")
    return df


def add_empirical_logp_pvalues(df: pd.DataFrame, logp_col: str, alpha: float, tail: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    pools = {str(pid): np.sort(g[logp_col].dropna().to_numpy(dtype=float)) for pid, g in df.groupby("pool_id")}
    if "GLOBAL" not in pools:
        pools["GLOBAL"] = np.sort(df[logp_col].dropna().to_numpy(dtype=float))

    p_low = np.empty(len(df), dtype=float)
    pool_ids = df["pool_id"].astype(str).to_numpy()
    vals = df[logp_col].to_numpy(dtype=float)
    for i, (pid, x) in enumerate(zip(pool_ids, vals)):
        arr = pools.get(pid, pools["GLOBAL"])
        p_low[i] = empirical_cdf(arr, float(x))
    df["p_emp_low"] = p_low
    df["p_emp_high"] = 1.0 - p_low

    tail = str(tail).lower().strip()
    if tail == "low":
        df["p_value"] = df["p_emp_low"]
    elif tail == "high":
        df["p_value"] = df["p_emp_high"]
    else:
        raise ValueError("--tail must be low or high")
    df["observable"] = df["p_value"] < float(alpha)

    cutoff_rows = []
    for pid, arr in pools.items():
        if arr.size == 0:
            continue
        cutoff_rows.append({
            "pool_id": pid,
            "alpha": float(alpha),
            "tail": tail,
            "logP_cutoff_low_tail": float(np.quantile(arr, float(alpha))),
            "logP_cutoff_high_tail": float(np.quantile(arr, 1.0 - float(alpha))),
            "n_rows_pool": int(arr.size),
        })
    return df, pd.DataFrame(cutoff_rows)


def aggregate_params(params_csv: Path, results_dir: Path, choice: str) -> Optional[Path]:
    params = pd.read_csv(params_csv)
    if "success" in params.columns:
        params = params[params["success"] == True].copy()
    if params.empty:
        return None
    cols = [c for c in ["gamma", "k", "fmin", "N_total_mean", "N_total_min", "N_total_max", "grid_size", "chunk_size", "log_likelihood"] if c in params.columns]
    if not cols:
        return None
    out = Path(results_dir) / "latent_params_aggregated.csv"
    if choice == "subject" and "subject" in params.columns:
        agg = params.groupby("subject")[cols].median(numeric_only=True).reset_index()
    else:
        med = params[cols].median(numeric_only=True)
        agg = pd.DataFrame([med])
        agg.insert(0, "pool_id", "GLOBAL")
        agg.insert(1, "aggregation_level", choice)
        agg["n_pairs_used"] = int(len(params))
    agg.to_csv(out, index=False)
    return out


def write_downstream_outputs(df: pd.DataFrame, clone_col: str, null_choice: str, results_dir: Path) -> Tuple[Path, Path]:
    results_dir = Path(results_dir)
    # Per clone/time observability summary.
    obs = (
        df.groupby(["subject", clone_col, "time"], as_index=False)
        .agg(
            observable=("observable", "max"),
            p_emp_low=("p_emp_low", "min"),
            p_value=("p_value", "min"),
            n_rows=(clone_col, "size"),
        )
    )
    out_obs = results_dir / f"latent_trajectory_observability_{null_choice}.csv"
    obs.to_csv(out_obs, index=False)

    # State summary, one row per subject/clone/time. Since each time already
    # corresponds to one replicate pair, first/median are equivalent; median is
    # robust if duplicated rows occur. State-level detectability is propagated as
    # the primary continuous quality metric; replicate-specific values remain QC.
    agg_dict = {}
    for c in [
        "count_rep1", "count_rep2", "depth_rep1", "depth_rep2",
        "freq_rep1", "freq_rep2", "freq_geo",
        "f_latent_mean", "f_latent_median", "f_latent_mode", "f_latent_q025", "f_latent_q975",
        "x_latent_mean", "x_latent_median", "x_latent_mode", "x_latent_q025", "x_latent_q975", "x_latent_sd",
        "p_detect_state", "p_dropout_state",
        "p_detect_rep1", "p_detect_rep2", "p_dropout_rep1", "p_dropout_rep2",
        "logP", "p_value", "observable",
    ]:
        if c in df.columns:
            if c == "observable":
                agg_dict[c] = (c, "max")
            elif c.startswith("count_") or c.startswith("depth_"):
                agg_dict[c] = (c, "first")
            else:
                agg_dict[c] = (c, "median")
    freq = df.groupby(["subject", clone_col, "time"], as_index=False).agg(**agg_dict)
    out_freq = results_dir / f"latent_trajectory_frequency_{null_choice}.csv"
    freq.to_csv(out_freq, index=False)
    return out_obs, out_freq


def run_aggregate(
    results_dir: Path,
    perclone_dir: Path,
    null_choice: str,
    alpha: float,
    tail: str,
    min_pairs_per_subject: int,
    intermediate_format: str,
    parquet_compression: str,
    file_sep: str,
) -> Path:
    df = load_perclone_latent(perclone_dir, file_sep=file_sep)
    clone_col, logp_col = infer_clone_and_logp_cols(df)
    df = assign_pools(df, null_choice, min_pairs_per_subject)
    df, cutoffs = add_empirical_logp_pvalues(df, logp_col=logp_col, alpha=alpha, tail=tail)

    ensure_dir(results_dir)
    cutoffs_out = Path(results_dir) / f"latent_logP_cutoffs_{null_choice}.csv"
    cutoffs.to_csv(cutoffs_out, index=False)

    main_out = write_table(
        df,
        Path(results_dir) / f"per_clone_latent_{null_choice}",
        fmt=intermediate_format,
        compression=parquet_compression,
        index=False,
    )
    out_obs, out_freq = write_downstream_outputs(df, clone_col, null_choice, Path(results_dir))

    print("Saved aggregate latent outputs:")
    print(" -", main_out)
    print(" -", cutoffs_out)
    print(" -", out_obs)
    print(" -", out_freq)
    return main_out


# =============================================================================
# CLI
# =============================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Latent-frequency parallel Step 2 for ClonoDynamics."
    )
    p.add_argument("--data-dir", type=Path, help="Directory containing raw replicate files.")
    p.add_argument("--results-dir", type=Path, required=True, help="Output directory.")
    p.add_argument("--perclone-dir", type=Path, default=None, help="Optional existing per-clone latent directory for --aggregate-only.")
    p.add_argument("--aggregate-only", action="store_true", help="Skip fitting and aggregate existing per-clone latent files.")

    p.add_argument("--file-sep", default="\t", help="Separator for TSV-like input files.")
    p.add_argument("--pattern", default=PATTERN_DEFAULT, help="Filename regex with subject/time/replica groups.")

    p.add_argument("--gamma-init", type=float, default=1.6)
    p.add_argument("--k-init", type=float, default=50.0)
    p.add_argument("--grid-size", type=int, default=500)
    p.add_argument("--chunk-size", type=int, default=20000)
    p.add_argument("--maxiter", type=int, default=None)
    p.add_argument("--optimizer-ftol", type=float, default=None)

    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--null", choices=["clone", "subject", "global"], default="subject")
    p.add_argument("--tail", choices=["low", "high"], default="low")
    p.add_argument("--min-pairs-per-subject", type=int, default=2)

    p.add_argument("--use-qc", action="store_true", help="Deprecated compatibility flag; ignored. Step 2 always uses all complete pairs.")
    p.add_argument("--qc-csv", type=Path, default=None, help="Deprecated compatibility option; ignored.")

    p.add_argument("--intermediate-format", choices=["parquet", "feather", "csv"], default="parquet")
    p.add_argument("--parquet-compression", default="snappy")
    p.add_argument("--write-csv-compat", action="store_true")

    p.add_argument("--n-jobs", type=int, default=1, help="Workers. Use <=0 for automatic macOS-friendly choice.")
    p.add_argument("--n-jobs-profile", choices=["memory", "balanced", "speed"], default="balanced")
    p.add_argument("--fit-cache-dir", type=Path, default=None)
    p.add_argument("--aggregate-results-dir", type=Path, default=None)
    p.add_argument("--skip-fit", action="store_true")

    p.add_argument("--ignore-qc", action="store_true", default=True)
    p.add_argument("--skip-alpha-sensitivity", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # ============================================================
    # AGGREGATE-ONLY MODE
    # ============================================================
    # Used for alpha sensitivity, e.g. p_02 / p_01.
    # It reads fitted per-clone latent files from the main alpha cache
    # and writes re-aggregated outputs to the requested alpha directory.
    if args.aggregate_only:
        if args.fit_cache_dir is None:
            raise ValueError("--aggregate-only requires --fit-cache-dir")

        if args.aggregate_results_dir is None:
            raise ValueError("--aggregate-only requires --aggregate-results-dir")

        if args.perclone_dir is not None:
            perclone_dir = Path(args.perclone_dir)
        else:
            perclone_dir = Path(args.fit_cache_dir) / "per_clone_latent_logP"

        results_dir = ensure_dir(Path(args.aggregate_results_dir))

        if not perclone_dir.exists():
            raise FileNotFoundError(
                f"Per-clone latent cache directory not found: {perclone_dir}"
            )

        run_aggregate(
            results_dir=results_dir,
            perclone_dir=perclone_dir,
            null_choice=args.null,
            alpha=args.alpha,
            tail=args.tail,
            min_pairs_per_subject=args.min_pairs_per_subject,
            intermediate_format=args.intermediate_format,
            parquet_compression=args.parquet_compression,
            file_sep=args.file_sep,
        )
        return 0

    # ============================================================
    # FIT + AGGREGATE MODE
    # ============================================================
    # Used for the main alpha, e.g. p_05.
    # It fits noiseK latent from raw replicate files, saves per-pair
    # latent outputs, then aggregates them.
    if args.data_dir is None:
        raise ValueError("--data-dir is required unless --aggregate-only is used")

    results_dir = ensure_dir(Path(args.results_dir))

    if args.perclone_dir is not None:
        perclone_dir = ensure_dir(Path(args.perclone_dir))
    else:
        perclone_dir = ensure_dir(results_dir / "per_clone_latent_logP")

    allowed_pairs = None
    if args.use_qc or args.qc_csv is not None:
        print("[warning] --use-qc/--qc-csv are ignored: this Step 2 version always processes all complete pairs.")

    n_jobs = resolve_macos_n_jobs(args.n_jobs, args.n_jobs_profile)

    params_csv = run_latent_fit(
        data_dir=args.data_dir,
        results_dir=results_dir,
        perclone_dir=perclone_dir,
        file_sep=args.file_sep,
        pattern_str=args.pattern,
        gamma_init=args.gamma_init,
        k_init=args.k_init,
        grid_size=args.grid_size,
        chunk_size=args.chunk_size,
        allowed_pairs=allowed_pairs,
        intermediate_format=args.intermediate_format,
        parquet_compression=args.parquet_compression,
        write_csv_compat=args.write_csv_compat,
        n_jobs=n_jobs,
        maxiter=args.maxiter,
        optimizer_ftol=args.optimizer_ftol,
    )

    build_latent_pair_qc(params_csv, results_dir)

    agg_params = aggregate_params(params_csv, results_dir, args.null)
    if agg_params:
        print("Saved aggregated params:", agg_params)

    run_aggregate(
        results_dir=results_dir,
        perclone_dir=perclone_dir,
        null_choice=args.null,
        alpha=args.alpha,
        tail=args.tail,
        min_pairs_per_subject=args.min_pairs_per_subject,
        intermediate_format=args.intermediate_format,
        parquet_compression=args.parquet_compression,
        file_sep=args.file_sep,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
