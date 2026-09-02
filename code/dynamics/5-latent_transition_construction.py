#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5-latent_transition_construction.py

Memory-bounded construction and posterior propagation of finite-time latent
clonotype transitions for the ClonoDynamics pipeline.

Overview
--------
This script implements Step 5 of the ClonoDynamics pipeline.

It converts the standardized longitudinal latent-state trajectories generated
by Step 4 (`4-latent_trajectory_construction.py`) into explicit finite-time
transitions between pairs of sampling time points.

For each clonotype represented at two eligible time points t0 < t1, the script
constructs a transition containing the endpoint latent states, the finite-time
latent displacement, the temporal lag, alternative conditioning coordinates,
endpoint observability and detectability, optional count/observed-frequency
representations, and propagated posterior uncertainty.

When the complete discrete Step-2 posterior files are supplied, the script
also propagates endpoint uncertainty to the displacement distribution by
Monte-Carlo sampling from the two endpoint posteriors.

The script is designed for very large longitudinal repertoire datasets. It is
memory-bounded: one subject and one temporal pair are processed at a time,
intermediate temporal-pair outputs are resumable Parquet parts, posterior
propagation is performed in Arrow batches, and completed parts are concatenated
without loading the complete transition table into memory.

This step constructs transition-level quantities but does not estimate drift,
conditional variance, diffusion coefficients or temporal-scaling laws. Those
belong to subsequent ClonoDynamics analysis steps.


Position in the ClonoDynamics pipeline
--------------------------------------

    Step 1
    1-repertoire_characterization.py
        -> observed repertoire architecture and heavy-tail characterization

    Step 2
    2-noise_aware_latent_inference.py
        -> posterior latent clonotype state, uncertainty and detectability

    Step 3
    3-latent_posterior_quality_characterization.py
        -> posterior quality, replicate agreement and detectability validation

    Step 4
    4-latent_trajectory_construction.py
        -> standardized subject x clonotype x time latent trajectories

    Step 5
    5-latent_transition_construction.py
        -> finite-time transitions (x0, x1, Delta x, dt)
        -> conditioning coordinates
        -> observability classes
        -> endpoint uncertainty/detectability
        -> optional full-posterior displacement propagation


Input
-----
The primary input is the standardized long-format trajectory table produced by
Step 4 and supplied with:

    --trajectories

Canonical input:

    latent_trajectories_long.parquet

Supported trajectory formats:

    .parquet
    .pq
    .csv
    .tsv

Parquet is strongly recommended for large analyses.


Required trajectory identifiers
--------------------------------
The trajectory table must contain identifiers for:

    subject
        subject identifier

    aaSeqCDR3
        amino-acid CDR3 clonotype identifier

    time
        sampling-time identifier

The input column names can be changed using:

    --subject-col
    --clone-col
    --time-col

with defaults:

    subject
    aaSeqCDR3
    time


Latent abundance input
----------------------
The script searches for a latent-frequency coordinate in the following order:

    --freq-col
    freq_latent
    freq
    f_latent_median
    f_latent_mean

and for a latent log-frequency coordinate in:

    --log-freq-col
    x_latent
    log_freq
    x_latent_median
    x_latent_mean

At least one usable latent abundance representation must be available.

If latent log-frequency is unavailable but latent frequency is present, the
script reconstructs

    x = ln(max(f, f_min))

where `f_min` is controlled by:

    --min-freq

and defaults to:

    1e-15.

For the primary ClonoDynamics analysis, the intended coordinate is the latent
log-frequency posterior summary propagated from Step 2 through Step 4.


State-level detectability input
-------------------------------
By default the trajectory table must contain:

    p_detect_state
    p_dropout_state

These are the posterior-predictive detectability and dropout probabilities of
the latent state inferred in Step 2 and propagated through Step 4.

Legacy trajectory tables without these quantities can be accepted only by
explicitly using:

    --allow-missing-state-detectability

No detectability model is fitted in this step.


Observability input
-------------------
The script uses the trajectory-level field:

    observable

(or a custom column supplied with `--observable-col`).

For each transition it defines endpoint observability as:

    obs0 = T if the t0 state is observable, otherwise F
    obs1 = T if the t1 state is observable, otherwise F

and the transition class:

    obs_class = obs0 + obs1

with four possible classes:

    TT    observable at both endpoints
    TF    observable at t0 but not t1
    FT    not observable at t0 but observable at t1
    FF    not observable at either endpoint

A subset can optionally be retained with:

    --keep-obs-class

If no class restriction is specified, all classes are retained.


Optional complete posterior input
---------------------------------
Complete latent-frequency posterior propagation is enabled by supplying:

    --posterior-dir

The expected Step-2 posterior file for subject s and time t is:

    latent_per_clone_logP_<subject>_<time>.posterior.npz

Each NPZ file must contain:

    f_grid
        one-dimensional latent-frequency grid

    posterior_matrix
        discrete posterior probabilities for each clonotype across f_grid

    aaSeqCDR3
        clonotype identifiers corresponding to posterior_matrix rows

The script transforms the grid to natural-log space as:

    x_grid = ln(f_grid).

Complete posterior propagation is optional. If `--posterior-dir` is omitted,
transition point estimates and endpoint-summary uncertainty quantities are
still constructed.


Canonical state preparation
---------------------------
Before transitions are formed, each subject is processed independently.
Trajectory rows are grouped by:

    subject x aaSeqCDR3 x time

so that at most one canonical state remains for each clonotype-timepoint.

If duplicate rows exist:

    - numeric variables are averaged;
    - `observable` is aggregated using the Boolean maximum.

This grouping is ordinarily an identity operation for a canonical Step-4
trajectory table containing one row per subject-clonotype-time state.


Mathematical definitions
------------------------

1. Finite-time transition
~~~~~~~~~~~~~~~~~~~~~~~~~
For clonotype i from subject s represented at two sampling times t0 < t1,
let

    x0 = x_i,s(t0)
    x1 = x_i,s(t1)

be the latent log-frequency states.

The temporal lag is:

    dt = t1 - t0.

A finite-time transition is therefore represented by:

    T_i,s(t0 -> t1) = {x0, x1, dt}.


2. Latent displacement
~~~~~~~~~~~~~~~~~~~~~~
The primary finite-time latent displacement is:

    Delta x = x1 - x0.

Output fields:

    x0_latent
    x1_latent
    dx_latent

Positive `dx_latent` indicates an increase in latent log-frequency over the
interval; negative values indicate a decrease.

The corresponding linear-frequency endpoint fields are:

    freq0_latent
    freq1_latent.


3. Initial-state conditioning coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The conventional forward conditioning coordinate is:

    x0_condition = x0.

This permits downstream estimation of conditional quantities such as:

    E[Delta x | x0].

This script constructs the coordinate but does not estimate the conditional
drift function.


4. Midpoint conditioning coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The endpoint-symmetric midpoint is:

              x0 + x1
    xmid = ------------.
                 2

It is stored as:

    xmid_latent.


5. Uncertainty-weighted conditioning coordinate (xstar)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When posterior log-frequency standard deviations are available at both
endpoints, let:

    sigma0 = SD[x(t0) | data]
    sigma1 = SD[x(t1) | data].

The script defines inverse-variance weights:

         1
    w0 = --------
         sigma0^2

         1
    w1 = --------
         sigma1^2

(with a small numerical stabilizer in the implementation), and constructs:

             w0*x0 + w1*x1
    xstar = ----------------.
                 w0 + w1

Output:

    xstar_latent

The source is recorded as:

    inverse_variance_weighted

when valid endpoint uncertainties are available. Otherwise:

    xstar = xmid

and:

    xstar_source = midpoint_fallback.


6. Endpoint-summary displacement uncertainty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When endpoint posterior standard deviations are available, the script
constructs:

    dx_sd = sqrt(sigma0^2 + sigma1^2)

and:

    dx_var = sigma0^2 + sigma1^2.

Output fields:

    x0_sd
    x1_sd
    dx_sd
    dx_var

This calculation treats the endpoint uncertainty contributions as independent.


7. Approximate Gaussian directional probability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Using the point-estimate displacement and endpoint-summary uncertainty:

                 Delta x
    z = --------------------------
        sqrt(sigma0^2 + sigma1^2)

and, under a Gaussian approximation:

    p_positive_approx = Phi(z)
    p_negative_approx = 1 - Phi(z),

where Phi is the standard-normal cumulative distribution function.

These are approximation-based diagnostics and must be distinguished from the
full-posterior Monte-Carlo directional probabilities described below.


8. Endpoint credible intervals and conservative displacement bounds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When endpoint posterior quantiles are available, the script preserves:

    x0_q025
    x0_q975
    x1_q025
    x1_q975.

It also constructs:

    dx_q025 = x1_q025 - x0_q975
    dx_q975 = x1_q975 - x0_q025.

These are endpoint-derived conservative bounds. They are NOT the exact 2.5th
and 97.5th percentiles of the posterior distribution of Delta x.

When complete posterior propagation is enabled, the direct Monte-Carlo
posterior displacement quantiles are instead:

    dx_post_q025
    dx_post_q975.


9. Observed-frequency displacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When observed-frequency coordinates are available at both endpoints:

    dx_obs = x1_obs - x0_obs.

The backward-compatible alias:

    dx_freq_obs

contains the same quantity.

Endpoint outputs include:

    x0_obs
    x1_obs
    freq0_obs
    freq1_obs.


10. Count-based displacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When summed-count log-coordinates are available:

    x_count_sum(t) = ln(count_sum(t) + 1),

then:

    dx_count_sum = x_count_sum(t1) - x_count_sum(t0)

and:

    xmid_count_sum = [x_count_sum(t0) + x_count_sum(t1)] / 2.


11. Closure component
~~~~~~~~~~~~~~~~~~~~~
When both observed-frequency and count-based displacements are available:

    closure_component = dx_obs - dx_count_sum.

This quantity separates the difference between the observed relative-frequency
change and the corresponding summed-count change, including denominator/depth
contributions introduced by the relative-frequency representation.


Transition selection
--------------------

12. All-pairs mode
~~~~~~~~~~~~~~~~~~
With:

    --transition-mode all

the script considers every subject-level pair of sampling times satisfying:

    t0 < t1

and:

    min_dt <= t1 - t0 <= max_dt.

For each temporal pair, endpoint states are inner-joined by `aaSeqCDR3`, so a
transition is formed only for clonotypes represented at both endpoints.


13. Adjacent mode
~~~~~~~~~~~~~~~~~
With:

    --transition-mode adjacent

or:

    --transition-mode adj

transitions are constructed between consecutive available states of each
individual clonotype after sorting its trajectory by time.

Adjacency is therefore clonotype-specific. The resulting transition is kept
only if its temporal separation satisfies the selected dt range.


14. Both mode
~~~~~~~~~~~~~
With:

    --transition-mode both

both all-pairs and clone-wise adjacent transitions are generated.

The output field:

    pattern

records either:

    all
    adj.


Complete posterior propagation
------------------------------

15. Monte-Carlo displacement posterior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When `--posterior-dir` is supplied, the script loads the discrete endpoint
posteriors:

    P0(x) = p(x_i(t0) | data_t0)
    P1(x) = p(x_i(t1) | data_t1).

For each transition it independently samples:

    x0^(m) ~ P0
    x1^(m) ~ P1

for:

    m = 1, ..., M,

where M is controlled by:

    --n-posterior-samples

and defaults to:

    M = 1000.

For each Monte-Carlo draw:

    Delta x^(m) = x1^(m) - x0^(m).

This produces a direct empirical approximation to the posterior distribution
of the finite-time displacement under independent endpoint posteriors.


16. Posterior transition summaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Monte-Carlo displacement propagation generates:

    x0_post_mean
    x1_post_mean

    x0_post_median
    x1_post_median

    dx_post_mean
    dx_post_median

    dx_post_q025
    dx_post_q975

    dx_post_sd

    p_dx_gt0
    p_dx_lt0.

The directional probabilities are:

                 1
    p_dx_gt0 = ----- sum_m I(Delta x^(m) > 0)
                 M

and:

                 1
    p_dx_lt0 = ----- sum_m I(Delta x^(m) < 0).
                 M

These full-posterior quantities are distinct from the Gaussian summary-based
`p_positive_approx` and `p_negative_approx` fields.


17. Posterior propagation status
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For every transition, the script records:

    posterior_success
    posterior_message.

Posterior propagation may fail if, for example:

    - an endpoint posterior file is unavailable;
    - the clonotype is absent from one endpoint posterior file;
    - posterior probabilities are empty or non-finite.

Such failures are recorded explicitly; they are not silently replaced by a
full-posterior estimate.


Deterministic reproducibility
-----------------------------
Posterior Monte-Carlo sampling uses a deterministic pair-level seed derived
from:

    base posterior seed
    subject
    t0
    t1
    transition pattern.

For fixed inputs and `--posterior-seed`, posterior propagation is therefore
reproducible independently of the processing order of temporal pairs.


Structural transition weights
-----------------------------
Optional structural weights can be generated with:

    --add-weights.

These weights do not change the transition variables. They provide optional
weights for downstream summaries.


18. Clone weight
~~~~~~~~~~~~~~~~
For clonotype i within one transition pattern, let n_i be the number of
transition rows contributed by that clonotype. Then:

              1
    w_clone = ---.
              n_i

This prevents clonotypes contributing many temporal pairs from automatically
receiving proportionally greater total weight.


19. Temporal-lag weight
~~~~~~~~~~~~~~~~~~~~~~~
When:

    --balance-dt

is enabled, let n_dt be the number of transition rows for temporal lag dt.
Then:

           1
    w_dt = ----.
          n_dt

Otherwise:

    w_dt = 1.


20. Observability-class weight
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The class weight is:

    w_class = 1
        for TT;

    w_class = gamma_cross
        for TF and FT;

    w_class = gamma_ff
        for FF.

Default values:

    gamma_cross = 0.25
    gamma_ff    = 0.10.

These values are user-configurable structural analysis weights, not posterior
probabilities.


21. Combined structural weight
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The final structural weight is:

    w = w_clone * w_dt * w_class.

If:

    --normalize-weight-components

is enabled, each positive component is divided by its dataset mean before the
three components are multiplied.


Replicate-decoupled conditioning
--------------------------------
This Step 5 implementation DOES NOT construct the replicate-decoupled AB/BA
forward-dynamics estimand.

Step 4 preserves replicate-specific posterior summaries so that a downstream
analysis can condition on one technical replicate while calculating the
finite-time displacement from the other, for example:

    conditioning: x0,A
    displacement: Delta x_B = x1,B - x0,B

and reciprocally:

    conditioning: x0,B
    displacement: Delta x_A = x1,A - x0,A.

Those cross-fitted AB/BA quantities are not generated by this script and must
not be inferred from `x0_condition`, `xmid_latent` or `xstar_latent`.


Filtering and inclusion policy
------------------------------
Transition construction itself does not impose an abundance-based quality
threshold.

A transition is determined by:

    - availability of valid endpoint latent states;
    - the selected transition mode;
    - the requested min/max temporal lag;
    - an optional explicit observability-class restriction.

TT, TF, FT and FF transitions can all be retained.

This allows later analyses to define estimand-specific subsets without
reconstructing the transition table.


Output
------
The principal output is the Parquet transition table specified with:

    --out

Canonical name:

    latent_transitions.parquet

The memory-bounded implementation requires Parquet as its primary output
format. A CSV compatibility copy can optionally be requested with:

    --write-csv-compat.


Core transition columns
-----------------------
The principal transition fields include:

    subject
    aaSeqCDR3

    t0
    t1
    dt

    x0_latent
    x1_latent
    dx_latent

    freq0_latent
    freq1_latent

    x0_condition
    xmid_latent
    xstar_latent
    xstar_source

    obs0
    obs1
    obs_class

    pattern
    n_steps.


Observed-state outputs
----------------------
When available:

    x0_obs
    x1_obs
    dx_obs
    dx_freq_obs

    freq0_obs
    freq1_obs.


Count/depth outputs
-------------------
When available, endpoint quantities include:

    count_rep1_t0 / count_rep1_t1
    count_rep2_t0 / count_rep2_t1
    count_sum_t0  / count_sum_t1

    depth_rep1_t0 / depth_rep1_t1
    depth_rep2_t0 / depth_rep2_t1
    depth_sum_t0  / depth_sum_t1
    depth_mean_t0 / depth_mean_t1

    freq_count_sum_t0 / freq_count_sum_t1
    x_count_sum_t0    / x_count_sum_t1.

Derived transition fields include:

    dx_count_sum
    xmid_count_sum
    closure_component.


Endpoint uncertainty outputs
----------------------------
When available:

    x0_sd
    x1_sd
    dx_sd
    dx_var

    x0_q025
    x0_q975
    x1_q025
    x1_q975

    dx_q025
    dx_q975

    freq0_q025
    freq0_q975
    freq1_q025
    freq1_q975.


Detectability outputs
---------------------
When available upstream:

    p_detect_state_t0
    p_detect_state_t1

    p_dropout_state_t0
    p_dropout_state_t1

    p_detect_rep1_t0
    p_detect_rep1_t1

    p_detect_rep2_t0
    p_detect_rep2_t1

    p_dropout_rep1_t0
    p_dropout_rep1_t1

    p_dropout_rep2_t0
    p_dropout_rep2_t1.


Full-posterior outputs
----------------------
When `--posterior-dir` is supplied:

    x0_post_mean
    x1_post_mean
    x0_post_median
    x1_post_median

    dx_post_mean
    dx_post_median
    dx_post_q025
    dx_post_q975
    dx_post_sd

    p_dx_gt0
    p_dx_lt0

    posterior_success
    posterior_message.


Optional weight outputs
-----------------------
When `--add-weights` is enabled:

    w_clone
    w_dt
    w_class
    w.


Memory-bounded processing architecture
--------------------------------------
Processing is organized as follows:

    1. one subject is prepared at a time;

    2. one eligible temporal pair is constructed at a time;

    3. the raw temporal-pair transition block is written to temporary Parquet;

    4. optional complete-posterior propagation is performed in Arrow batches;

    5. the completed temporal pair is stored as a resumable Parquet part;

    6. completed parts are concatenated incrementally into the final Parquet
       output without loading them simultaneously.

This architecture allows transition datasets substantially larger than
available RAM to be generated.


Resume and restart
------------------
The work directory stores independently completed temporal-pair parts.

With resume enabled, completed parts are reused. Empty temporal-pair
combinations are represented by explicit JSON markers so that they are not
recomputed.

The work directory contains a manifest describing the input and key analysis
parameters. If a pre-existing work directory was created using incompatible
inputs or settings, execution stops instead of mixing results.

Use:

    --restart

to explicitly discard an existing work directory and rebuild the analysis.


Build report
------------
Unless an alternative path is supplied with `--report-md`, the script writes:

    latent_transition_build_report.md

The report records the trajectory input, number of prepared states and
subjects, latent-coordinate source, transition mode and lag limits,
observability restriction, posterior-propagation settings, completed/empty
temporal-pair parts, counts by dt and observability class, posterior success,
and availability of endpoint intervals, xstar, closure decomposition,
state-level detectability and full-posterior propagation.


Typical usage
-------------
Primary all-pairs ClonoDynamics analysis:

    python3 5-latent_transition_construction.py \
        --trajectories \
        ./results/4-latent_trajectory_construction/latent_trajectories_long.parquet \
        --posterior-dir \
        ./results/2-noise_aware_latent_inference/per_clone_latent_logP \
        --out \
        ./results/5-latent_transition_construction/latent_transitions.parquet \
        --transition-mode all \
        --min-dt 1 \
        --max-dt 5 \
        --n-posterior-samples 1000 \
        --posterior-seed 123 \
        --posterior-chunk-size 5000 \
        --parquet-compression zstd \
        --parquet-row-group-size 100000

Adjacent-transition analysis:

    python3 5-latent_transition_construction.py \
        --trajectories \
        ./results/4-latent_trajectory_construction/latent_trajectories_long.parquet \
        --posterior-dir \
        ./results/2-noise_aware_latent_inference/per_clone_latent_logP \
        --out \
        ./results/5-latent_transition_construction/latent_transitions_adj.parquet \
        --transition-mode adjacent \
        --min-dt 1 \
        --max-dt 5


Interpretation
--------------
Step 4 represents longitudinal clonotype states:

    x(t1), x(t2), ..., x(tT).

Step 5 converts pairs of those states into finite-time observations:

    (x0, x1, Delta x, dt).

It therefore marks the transition from latent-state reconstruction to
dynamical analysis. The resulting transition table is the direct input for
subsequent ClonoDynamics analyses of conditional drift, fluctuation magnitude,
temporal scaling and representation sensitivity.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

try:
    import polars as pl
except ImportError as exc:
    raise SystemExit(
        "polars is required. Activate the environment containing Polars 1.36.x."
    ) from exc

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError as exc:
    raise SystemExit("pyarrow is required for incremental Parquet writing.") from exc


PARQUET_SUFFIXES = {".parquet", ".pq"}
CSV_SUFFIXES = {".csv", ".tsv"}
OBS_CLASSES = {"TT", "TF", "FT", "FF"}
TRUE_STRINGS = {
    "true",
    "t",
    "1",
    "yes",
    "y",
    "observable",
    "observed",
    "detected",
    "detectable",
}

POSTERIOR_FLOAT_COLUMNS = [
    "x0_post_mean",
    "x1_post_mean",
    "x0_post_median",
    "x1_post_median",
    "dx_post_mean",
    "dx_post_median",
    "dx_post_q025",
    "dx_post_q975",
    "dx_post_sd",
    "p_dx_gt0",
    "p_dx_lt0",
]

COUNT_COLUMNS = [
    "count_rep1",
    "count_rep2",
    "depth_rep1",
    "depth_rep2",
    "count_sum",
    "depth_sum",
    "depth_mean",
    "freq_count_sum",
    "x_count_sum",
]

QUALITY_COLUMNS = [
    "p_detect_state",
    "p_dropout_state",
    "p_dropout_rep1",
    "p_dropout_rep2",
    "p_detect_rep1",
    "p_detect_rep2",
]

UNCERTAINTY_COLUMNS = [
    "freq_q025",
    "freq_q975",
    "x_q025",
    "x_q975",
    "x_sd",
]


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def eprint(*args: Any, **kwargs: Any) -> None:
    print(*args, file=sys.stderr, **kwargs)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def remove_if_exists(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def collect_streaming(lazy_frame: pl.LazyFrame) -> pl.DataFrame:
    """Collect using the Polars streaming engine with a compatibility fallback."""
    try:
        return lazy_frame.collect(engine="streaming")
    except TypeError:
        return lazy_frame.collect(streaming=True)


def sink_parquet_streaming(
    lazy_frame: pl.LazyFrame,
    path: Path,
    compression: str,
    row_group_size: int,
) -> None:
    """Write a LazyFrame to Parquet without a global collect."""
    ensure_parent(path)
    kwargs: Dict[str, Any] = {
        "compression": compression,
        "statistics": True,
        "maintain_order": True,
    }
    if row_group_size > 0:
        kwargs["row_group_size"] = int(row_group_size)

    try:
        lazy_frame.sink_parquet(path, engine="streaming", mkdir=True, **kwargs)
    except TypeError:
        try:
            lazy_frame.sink_parquet(path, **kwargs)
        except TypeError:
            kwargs.pop("row_group_size", None)
            lazy_frame.sink_parquet(path, **kwargs)


def scan_table(path: Path, separator: str) -> pl.LazyFrame:
    suffix = path.suffix.lower()
    if suffix in PARQUET_SUFFIXES:
        return pl.scan_parquet(path)
    if suffix == ".tsv":
        return pl.scan_csv(path, separator="\t", infer_schema_length=10000)
    if suffix == ".csv":
        return pl.scan_csv(path, separator=separator, infer_schema_length=10000)
    raise ValueError(
        f"Unsupported trajectory format: {path.suffix}. Use Parquet, CSV or TSV."
    )


def parse_keep_values(text: Optional[str]) -> Optional[Set[str]]:
    if text is None or not str(text).strip():
        return None
    values = {item.strip().upper() for item in str(text).split(",") if item.strip()}
    invalid = values.difference(OBS_CLASSES)
    if invalid:
        raise ValueError(f"Invalid observability classes: {sorted(invalid)}")
    return values


def first_existing(names: Set[str], candidates: Sequence[str]) -> Optional[str]:
    for candidate in candidates:
        if candidate in names:
            return candidate
    return None


def numeric_expr(source: str, alias: Optional[str] = None) -> pl.Expr:
    return pl.col(source).cast(pl.Float64, strict=False).alias(alias or source)


def is_numeric_dtype(dtype: Any) -> bool:
    numeric_types = {
        pl.Int8,
        pl.Int16,
        pl.Int32,
        pl.Int64,
        pl.UInt8,
        pl.UInt16,
        pl.UInt32,
        pl.UInt64,
        pl.Float32,
        pl.Float64,
    }
    return dtype in numeric_types


def observable_expr(
    source: Optional[str],
    schema: Any,
) -> pl.Expr:
    if source is None or source not in schema.names():
        return pl.lit(False).alias("observable")

    dtype = schema[source]
    if dtype == pl.Boolean:
        return pl.col(source).fill_null(False).cast(pl.Boolean).alias("observable")

    if is_numeric_dtype(dtype):
        return (
            pl.col(source)
            .cast(pl.Float64, strict=False)
            .fill_null(0.0)
            .ne(0.0)
            .alias("observable")
        )

    return (
        pl.col(source)
        .cast(pl.String, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
        .is_in(sorted(TRUE_STRINGS))
        .fill_null(False)
        .alias("observable")
    )


def input_signature(path: Path) -> Dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


# -----------------------------------------------------------------------------
# Canonical trajectory preparation
# -----------------------------------------------------------------------------


def build_canonical_expressions(
    schema: Any,
    args: argparse.Namespace,
) -> Tuple[List[pl.Expr], str, List[str]]:
    names = set(schema.names())

    required = [args.subject_col, args.clone_col, args.time_col]
    missing = [name for name in required if name not in names]
    if missing:
        raise ValueError(
            f"Trajectory table is missing required columns: {missing}. "
            f"Available columns: {sorted(names)}"
        )

    if not args.allow_missing_state_detectability:
        state_missing = [
            name
            for name in ["p_detect_state", "p_dropout_state"]
            if name not in names
        ]
        if state_missing:
            raise ValueError(
                "State-level detectability columns are missing: "
                f"{state_missing}. Use --allow-missing-state-detectability only "
                "for legacy Step-4 outputs."
            )

    freq_used = first_existing(
        names,
        [
            args.freq_col,
            "freq_latent",
            "freq",
            "f_latent_median",
            "f_latent_mean",
        ],
    )
    log_used = first_existing(
        names,
        [
            args.log_freq_col,
            "x_latent",
            "log_freq",
            "x_latent_median",
            "x_latent_mean",
        ],
    )

    if freq_used is None and log_used is None:
        raise ValueError(
            "No latent frequency or latent log-frequency column was found."
        )

    expressions: List[pl.Expr] = [
        pl.col(args.subject_col).cast(pl.Int64, strict=False).alias("subject"),
        pl.col(args.clone_col).cast(pl.String, strict=False).alias("aaSeqCDR3"),
        pl.col(args.time_col).cast(pl.Int64, strict=False).alias("time"),
    ]

    if freq_used is not None:
        freq_expr = pl.col(freq_used).cast(pl.Float64, strict=False)
        expressions.append(freq_expr.alias("freq_latent"))
    else:
        freq_expr = pl.lit(float("nan"), dtype=pl.Float64)
        expressions.append(freq_expr.alias("freq_latent"))

    if log_used is not None:
        expressions.append(
            pl.col(log_used).cast(pl.Float64, strict=False).alias("x_latent")
        )
        x_source = log_used
    else:
        expressions.append(
            freq_expr
            .clip(lower_bound=float(args.min_freq))
            .log()
            .alias("x_latent")
        )
        x_source = f"ln({freq_used})"

    if "freq_obs" in names:
        freq_obs_expr = pl.col("freq_obs").cast(pl.Float64, strict=False)
        expressions.append(freq_obs_expr.alias("freq_obs"))
    else:
        freq_obs_expr = pl.lit(float("nan"), dtype=pl.Float64)
        expressions.append(freq_obs_expr.alias("freq_obs"))

    if "x_obs" in names:
        expressions.append(
            pl.col("x_obs").cast(pl.Float64, strict=False).alias("x_obs")
        )
    elif "freq_obs" in names:
        expressions.append(
            freq_obs_expr
            .clip(lower_bound=float(args.min_freq))
            .log()
            .alias("x_obs")
        )
    else:
        expressions.append(pl.lit(float("nan"), dtype=pl.Float64).alias("x_obs"))

    available_numeric: List[str] = [
        "freq_latent",
        "x_latent",
        "freq_obs",
        "x_obs",
    ]

    # Original count/depth columns.
    for column in ["count_rep1", "count_rep2", "depth_rep1", "depth_rep2"]:
        if column in names:
            expressions.append(numeric_expr(column))
            available_numeric.append(column)

    if "count_sum" in names:
        count_sum_expr: Optional[pl.Expr] = pl.col("count_sum").cast(
            pl.Float64, strict=False
        )
        expressions.append(count_sum_expr.alias("count_sum"))
        available_numeric.append("count_sum")
    elif {"count_rep1", "count_rep2"}.issubset(names):
        count_sum_expr = (
            pl.col("count_rep1").cast(pl.Float64, strict=False).fill_null(0.0)
            + pl.col("count_rep2").cast(pl.Float64, strict=False).fill_null(0.0)
        )
        expressions.append(count_sum_expr.alias("count_sum"))
        available_numeric.append("count_sum")
    else:
        count_sum_expr = None

    if "depth_sum" in names:
        depth_sum_expr: Optional[pl.Expr] = pl.col("depth_sum").cast(
            pl.Float64, strict=False
        )
        expressions.append(depth_sum_expr.alias("depth_sum"))
        available_numeric.append("depth_sum")
    elif {"depth_rep1", "depth_rep2"}.issubset(names):
        depth_sum_expr = (
            pl.col("depth_rep1").cast(pl.Float64, strict=False).fill_null(0.0)
            + pl.col("depth_rep2").cast(pl.Float64, strict=False).fill_null(0.0)
        )
        expressions.append(depth_sum_expr.alias("depth_sum"))
        available_numeric.append("depth_sum")
    else:
        depth_sum_expr = None

    if "depth_mean" in names:
        expressions.append(numeric_expr("depth_mean"))
        available_numeric.append("depth_mean")
    elif {"depth_rep1", "depth_rep2"}.issubset(names):
        expressions.append(
            (
                0.5
                * (
                    pl.col("depth_rep1").cast(pl.Float64, strict=False)
                    + pl.col("depth_rep2").cast(pl.Float64, strict=False)
                )
            ).alias("depth_mean")
        )
        available_numeric.append("depth_mean")

    if "freq_count_sum" in names:
        expressions.append(numeric_expr("freq_count_sum"))
        available_numeric.append("freq_count_sum")
    elif count_sum_expr is not None and depth_sum_expr is not None:
        expressions.append(
            pl.when(depth_sum_expr > 0)
            .then(count_sum_expr / depth_sum_expr)
            .otherwise(pl.lit(float("nan"), dtype=pl.Float64))
            .alias("freq_count_sum")
        )
        available_numeric.append("freq_count_sum")

    if "x_count_sum" in names:
        expressions.append(numeric_expr("x_count_sum"))
        available_numeric.append("x_count_sum")
    elif count_sum_expr is not None:
        expressions.append(
            (count_sum_expr.clip(lower_bound=0.0) + 1.0)
            .log()
            .alias("x_count_sum")
        )
        available_numeric.append("x_count_sum")

    observable_source = (
        args.observable_col if args.observable_col in names else None
    )
    expressions.append(observable_expr(observable_source, schema))

    optional_sources = {
        "freq_q025": ["freq_q025"],
        "freq_q975": ["freq_q975"],
        "x_q025": ["x_q025", "log_freq_q025"],
        "x_q975": ["x_q975", "log_freq_q975"],
        "x_sd": ["x_sd", "log_freq_sd"],
        "p_detect_state": ["p_detect_state"],
        "p_dropout_state": ["p_dropout_state"],
        "p_dropout_rep1": ["p_dropout_rep1"],
        "p_dropout_rep2": ["p_dropout_rep2"],
        "p_detect_rep1": ["p_detect_rep1"],
        "p_detect_rep2": ["p_detect_rep2"],
    }

    for target, candidates in optional_sources.items():
        source = first_existing(names, candidates)
        if source is not None:
            expressions.append(numeric_expr(source, target))
            available_numeric.append(target)

    return expressions, x_source, available_numeric


def prepare_subject(
    trajectories: Path,
    separator: str,
    subject: int,
    args: argparse.Namespace,
    canonical_expressions: Sequence[pl.Expr],
    canonical_numeric: Sequence[str],
) -> pl.DataFrame:
    # Filter before the canonical projection so Parquet predicate pushdown can act
    # on the original subject column.
    lazy_frame = (
        scan_table(trajectories, separator)
        .filter(
            pl.col(args.subject_col).cast(pl.Int64, strict=False) == int(subject)
        )
        .select(list(canonical_expressions))
        .filter(
            pl.col("subject").is_not_null()
            & pl.col("time").is_not_null()
            & pl.col("aaSeqCDR3").is_not_null()
            & (pl.col("aaSeqCDR3").str.len_chars() > 0)
            & pl.col("x_latent").is_not_null()
        )
    )

    aggregations: List[pl.Expr] = [
        pl.col("observable").max().alias("observable")
    ]
    for column in canonical_numeric:
        aggregations.append(pl.col(column).mean().alias(column))

    prepared = (
        lazy_frame
        .group_by(["subject", "aaSeqCDR3", "time"])
        .agg(aggregations)
        .sort(["aaSeqCDR3", "time"])
    )
    return collect_streaming(prepared)


# -----------------------------------------------------------------------------
# Transition construction
# -----------------------------------------------------------------------------


def normal_cdf_approx_expr(z_value: pl.Expr) -> pl.Expr:
    """Vectorised approximation to the standard normal cumulative distribution."""
    a1 = 0.319381530
    a2 = -0.356563782
    a3 = 1.781477937
    a4 = -1.821255978
    a5 = 1.330274429
    p_value = 0.2316419

    absolute = z_value.abs()
    t_value = 1.0 / (1.0 + p_value * absolute)
    polynomial = (
        ((((a5 * t_value + a4) * t_value + a3) * t_value + a2) * t_value + a1)
        * t_value
    )
    density = (-0.5 * absolute.pow(2)).exp() / math.sqrt(2.0 * math.pi)
    positive_cdf = 1.0 - density * polynomial
    return (
        pl.when(z_value >= 0)
        .then(positive_cdf)
        .otherwise(1.0 - positive_cdf)
    )


def endpoint_frame(
    subject_frame: pl.DataFrame,
    time_value: int,
    suffix: str,
) -> pl.DataFrame:
    endpoint = subject_frame.filter(pl.col("time") == int(time_value)).drop(
        ["subject", "time"]
    )
    rename_map = {
        column: f"{column}_{suffix}"
        for column in endpoint.columns
        if column != "aaSeqCDR3"
    }
    return endpoint.rename(rename_map)


def transition_projection(
    joined: pl.DataFrame,
    subject: int,
    t0: int,
    t1: int,
    pattern: str,
    keep_classes: Optional[Set[str]],
) -> pl.DataFrame:
    dt_value = int(t1 - t0)

    observability0 = (
        pl.when(pl.col("observable_t0"))
        .then(pl.lit("T"))
        .otherwise(pl.lit("F"))
    )
    observability1 = (
        pl.when(pl.col("observable_t1"))
        .then(pl.lit("T"))
        .otherwise(pl.lit("F"))
    )

    expressions: List[pl.Expr] = [
        pl.lit(int(subject), dtype=pl.Int64).alias("subject"),
        pl.col("aaSeqCDR3"),
        pl.lit(int(t0), dtype=pl.Int64).alias("t0"),
        pl.lit(int(t1), dtype=pl.Int64).alias("t1"),
        pl.lit(dt_value, dtype=pl.Int64).alias("dt"),
        pl.col("x_latent_t0").alias("x0_latent"),
        pl.col("x_latent_t1").alias("x1_latent"),
        (pl.col("x_latent_t1") - pl.col("x_latent_t0")).alias("dx_latent"),
        pl.col("freq_latent_t0").alias("freq0_latent"),
        pl.col("freq_latent_t1").alias("freq1_latent"),
        pl.col("x_obs_t0").alias("x0_obs"),
        pl.col("x_obs_t1").alias("x1_obs"),
        (pl.col("x_obs_t1") - pl.col("x_obs_t0")).alias("dx_obs"),
        (pl.col("x_obs_t1") - pl.col("x_obs_t0")).alias("dx_freq_obs"),
        pl.col("freq_obs_t0").alias("freq0_obs"),
        pl.col("freq_obs_t1").alias("freq1_obs"),
        (0.5 * (pl.col("x_latent_t0") + pl.col("x_latent_t1"))).alias(
            "xmid_latent"
        ),
        pl.col("x_latent_t0").alias("x0_condition"),
        observability0.alias("obs0"),
        observability1.alias("obs1"),
        pl.concat_str([observability0, observability1]).alias("obs_class"),
        pl.lit(pattern).alias("pattern"),
        pl.lit(dt_value, dtype=pl.Int64).alias("n_steps"),
    ]

    for column in COUNT_COLUMNS:
        endpoint0 = f"{column}_t0"
        endpoint1 = f"{column}_t1"
        if endpoint0 in joined.columns:
            expressions.append(pl.col(endpoint0))
        if endpoint1 in joined.columns:
            expressions.append(pl.col(endpoint1))

    if {"x_count_sum_t0", "x_count_sum_t1"}.issubset(joined.columns):
        dx_count = pl.col("x_count_sum_t1") - pl.col("x_count_sum_t0")
        expressions.extend(
            [
                dx_count.alias("dx_count_sum"),
                (
                    0.5
                    * (pl.col("x_count_sum_t0") + pl.col("x_count_sum_t1"))
                ).alias("xmid_count_sum"),
                (
                    (pl.col("x_obs_t1") - pl.col("x_obs_t0")) - dx_count
                ).alias("closure_component"),
            ]
        )
    else:
        expressions.extend(
            [
                pl.lit(float("nan"), dtype=pl.Float64).alias("dx_count_sum"),
                pl.lit(float("nan"), dtype=pl.Float64).alias(
                    "xmid_count_sum"
                ),
                pl.lit(float("nan"), dtype=pl.Float64).alias(
                    "closure_component"
                ),
            ]
        )

    midpoint = 0.5 * (pl.col("x_latent_t0") + pl.col("x_latent_t1"))
    if {"x_sd_t0", "x_sd_t1"}.issubset(joined.columns):
        sd0 = pl.col("x_sd_t0")
        sd1 = pl.col("x_sd_t1")
        dx_sd = (sd0.pow(2) + sd1.pow(2)).sqrt()
        valid_sd = (
            sd0.is_finite()
            & sd1.is_finite()
            & (sd0 > 0)
            & (sd1 > 0)
        )
        weight0 = 1.0 / (sd0.pow(2) + 1e-12)
        weight1 = 1.0 / (sd1.pow(2) + 1e-12)
        xstar = (
            pl.when(valid_sd)
            .then(
                (
                    weight0 * pl.col("x_latent_t0")
                    + weight1 * pl.col("x_latent_t1")
                )
                / (weight0 + weight1)
            )
            .otherwise(midpoint)
        )
        z_value = (
            pl.col("x_latent_t1") - pl.col("x_latent_t0")
        ) / dx_sd
        positive_probability = (
            pl.when(dx_sd.is_finite() & (dx_sd > 0))
            .then(normal_cdf_approx_expr(z_value))
            .otherwise(pl.lit(float("nan"), dtype=pl.Float64))
        )

        expressions.extend(
            [
                xstar.alias("xstar_latent"),
                pl.when(valid_sd)
                .then(pl.lit("inverse_variance_weighted"))
                .otherwise(pl.lit("midpoint_fallback"))
                .alias("xstar_source"),
                sd0.alias("x0_sd"),
                sd1.alias("x1_sd"),
                dx_sd.alias("dx_sd"),
                dx_sd.pow(2).alias("dx_var"),
                positive_probability.alias("p_positive_approx"),
                (1.0 - positive_probability).alias("p_negative_approx"),
            ]
        )
    else:
        expressions.extend(
            [
                midpoint.alias("xstar_latent"),
                pl.lit("midpoint_fallback").alias("xstar_source"),
                pl.lit(float("nan"), dtype=pl.Float64).alias("x0_sd"),
                pl.lit(float("nan"), dtype=pl.Float64).alias("x1_sd"),
                pl.lit(float("nan"), dtype=pl.Float64).alias("dx_sd"),
                pl.lit(float("nan"), dtype=pl.Float64).alias("dx_var"),
                pl.lit(float("nan"), dtype=pl.Float64).alias(
                    "p_positive_approx"
                ),
                pl.lit(float("nan"), dtype=pl.Float64).alias(
                    "p_negative_approx"
                ),
            ]
        )

    if {"x_q025_t0", "x_q975_t0", "x_q025_t1", "x_q975_t1"}.issubset(
        joined.columns
    ):
        expressions.extend(
            [
                pl.col("x_q025_t0").alias("x0_q025"),
                pl.col("x_q975_t0").alias("x0_q975"),
                pl.col("x_q025_t1").alias("x1_q025"),
                pl.col("x_q975_t1").alias("x1_q975"),
                (pl.col("x_q025_t1") - pl.col("x_q975_t0")).alias(
                    "dx_q025"
                ),
                (pl.col("x_q975_t1") - pl.col("x_q025_t0")).alias(
                    "dx_q975"
                ),
            ]
        )
    else:
        for column in [
            "x0_q025",
            "x0_q975",
            "x1_q025",
            "x1_q975",
            "dx_q025",
            "dx_q975",
        ]:
            expressions.append(
                pl.lit(float("nan"), dtype=pl.Float64).alias(column)
            )

    for source in ["freq_q025", "freq_q975"]:
        endpoint0 = f"{source}_t0"
        endpoint1 = f"{source}_t1"
        output0 = source.replace("freq_", "freq0_")
        output1 = source.replace("freq_", "freq1_")
        if endpoint0 in joined.columns:
            expressions.append(pl.col(endpoint0).alias(output0))
        if endpoint1 in joined.columns:
            expressions.append(pl.col(endpoint1).alias(output1))

    for column in QUALITY_COLUMNS:
        endpoint0 = f"{column}_t0"
        endpoint1 = f"{column}_t1"
        if endpoint0 in joined.columns:
            expressions.append(pl.col(endpoint0))
        if endpoint1 in joined.columns:
            expressions.append(pl.col(endpoint1))

    output = joined.select(expressions)
    if keep_classes is not None:
        output = output.filter(pl.col("obs_class").is_in(sorted(keep_classes)))

    return output.sort(["subject", "aaSeqCDR3", "t0", "t1"])


def build_all_pair(
    subject_frame: pl.DataFrame,
    subject: int,
    t0: int,
    t1: int,
    keep_classes: Optional[Set[str]],
) -> pl.DataFrame:
    endpoint0 = endpoint_frame(subject_frame, t0, "t0")
    endpoint1 = endpoint_frame(subject_frame, t1, "t1")
    joined = endpoint0.join(endpoint1, on="aaSeqCDR3", how="inner")
    del endpoint0, endpoint1
    return transition_projection(
        joined=joined,
        subject=subject,
        t0=t0,
        t1=t1,
        pattern="all",
        keep_classes=keep_classes,
    )


def build_adjacent_blocks(
    subject_frame: pl.DataFrame,
    subject: int,
    min_dt: int,
    max_dt: Optional[int],
    keep_classes: Optional[Set[str]],
) -> Iterable[Tuple[int, int, pl.DataFrame]]:
    """Yield clone-wise adjacent transitions, grouped by temporal pair."""
    ordered = subject_frame.sort(["aaSeqCDR3", "time"])
    endpoint_columns = [
        column
        for column in ordered.columns
        if column not in {"subject", "aaSeqCDR3", "time"}
    ]

    shifted_expressions: List[pl.Expr] = [
        pl.col("time").shift(-1).over("aaSeqCDR3").alias("time_next")
    ]
    for column in endpoint_columns:
        shifted_expressions.append(
            pl.col(column)
            .shift(-1)
            .over("aaSeqCDR3")
            .alias(f"{column}_next")
        )

    shifted = (
        ordered
        .with_columns(shifted_expressions)
        .filter(pl.col("time_next").is_not_null())
        .with_columns((pl.col("time_next") - pl.col("time")).alias("_dt"))
        .filter(pl.col("_dt") >= int(min_dt))
    )
    if max_dt is not None:
        shifted = shifted.filter(pl.col("_dt") <= int(max_dt))

    pair_rows = (
        shifted
        .select(["time", "time_next"])
        .unique()
        .sort(["time", "time_next"])
        .iter_rows()
    )

    for raw_t0, raw_t1 in pair_rows:
        t0 = int(raw_t0)
        t1 = int(raw_t1)
        pair_block = shifted.filter(
            (pl.col("time") == t0) & (pl.col("time_next") == t1)
        )

        select_expressions: List[pl.Expr] = [pl.col("aaSeqCDR3")]
        for column in endpoint_columns:
            select_expressions.append(pl.col(column).alias(f"{column}_t0"))
            select_expressions.append(
                pl.col(f"{column}_next").alias(f"{column}_t1")
            )

        joined_like = pair_block.select(select_expressions)
        yield (
            t0,
            t1,
            transition_projection(
                joined=joined_like,
                subject=subject,
                t0=t0,
                t1=t1,
                pattern="adj",
                keep_classes=keep_classes,
            ),
        )


# -----------------------------------------------------------------------------
# Posterior propagation
# -----------------------------------------------------------------------------


def posterior_path(
    posterior_dir: Path,
    subject: int,
    time_value: int,
) -> Path:
    return posterior_dir / (
        f"latent_per_clone_logP_{int(subject)}_{int(time_value)}.posterior.npz"
    )


def load_posterior_file(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(str(path))

    with np.load(path, allow_pickle=True) as archive:
        required = {"f_grid", "posterior_matrix", "aaSeqCDR3"}
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(
                f"Posterior file {path.name} is missing arrays: {sorted(missing)}"
            )

        f_grid = np.asarray(archive["f_grid"], dtype=np.float64)
        posterior_matrix = np.asarray(
            archive["posterior_matrix"], dtype=np.float32
        )
        clones = np.asarray(archive["aaSeqCDR3"]).astype(str)

    if f_grid.ndim != 1:
        raise ValueError(f"f_grid must be one-dimensional in {path.name}")
    if posterior_matrix.ndim != 2:
        raise ValueError(
            f"posterior_matrix must be two-dimensional in {path.name}"
        )

    if posterior_matrix.shape[0] == f_grid.size and posterior_matrix.shape[1] == clones.size:
        posterior_matrix = posterior_matrix.T

    if posterior_matrix.shape[0] != clones.size:
        raise ValueError(
            f"Posterior rows ({posterior_matrix.shape[0]:,}) do not match "
            f"clones ({clones.size:,}) in {path.name}"
        )
    if posterior_matrix.shape[1] != f_grid.size:
        raise ValueError(
            f"Posterior columns ({posterior_matrix.shape[1]:,}) do not match "
            f"f_grid ({f_grid.size:,}) in {path.name}"
        )

    safe_grid = np.maximum(f_grid, np.finfo(np.float64).tiny)
    x_grid = np.log(safe_grid)
    clone_to_index = {clone: index for index, clone in enumerate(clones)}

    return {
        "x_grid": x_grid,
        "posterior_matrix": posterior_matrix,
        "clone_to_index": clone_to_index,
    }


def sample_posterior(
    x_grid: np.ndarray,
    probabilities: np.ndarray,
    n_samples: int,
    rng: np.random.Generator,
) -> np.ndarray:
    probabilities64 = np.asarray(probabilities, dtype=np.float64)
    probabilities64 = np.where(
        np.isfinite(probabilities64) & (probabilities64 > 0),
        probabilities64,
        0.0,
    )
    total = float(probabilities64.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Posterior probabilities are empty or non-finite")
    probabilities64 /= total
    return rng.choice(
        x_grid,
        size=int(n_samples),
        replace=True,
        p=probabilities64,
    )


def deterministic_pair_seed(
    base_seed: int,
    subject: int,
    t0: int,
    t1: int,
    pattern: str,
) -> int:
    payload = f"{base_seed}|{subject}|{t0}|{t1}|{pattern}"
    digest = hashlib.blake2b(payload.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="little", signed=False) % (2 ** 32)


def append_posterior_columns(
    record_batch: pa.RecordBatch,
    posterior0: Optional[Dict[str, Any]],
    posterior1: Optional[Dict[str, Any]],
    file_error: Optional[str],
    n_samples: int,
    rng: np.random.Generator,
) -> pa.Table:
    table = pa.Table.from_batches([record_batch])
    clones = table.column("aaSeqCDR3").to_pylist()
    row_count = len(clones)

    output_arrays: Dict[str, np.ndarray] = {
        column: np.full(row_count, np.nan, dtype=np.float64)
        for column in POSTERIOR_FLOAT_COLUMNS
    }
    success = np.zeros(row_count, dtype=np.bool_)
    messages: List[str] = [""] * row_count

    if file_error is not None or posterior0 is None or posterior1 is None:
        message = file_error or "posterior_file_error"
        messages = [message] * row_count
    else:
        map0 = posterior0["clone_to_index"]
        map1 = posterior1["clone_to_index"]
        grid0 = posterior0["x_grid"]
        grid1 = posterior1["x_grid"]
        matrix0 = posterior0["posterior_matrix"]
        matrix1 = posterior1["posterior_matrix"]

        posterior_ddof = 1 if int(n_samples) > 1 else 0

        for row_index, raw_clone in enumerate(clones):
            clone = str(raw_clone)
            try:
                index0 = map0.get(clone)
                index1 = map1.get(clone)
                if index0 is None:
                    raise KeyError(f"clone missing at t0: {clone}")
                if index1 is None:
                    raise KeyError(f"clone missing at t1: {clone}")

                samples0 = sample_posterior(
                    grid0,
                    matrix0[index0],
                    n_samples,
                    rng,
                )
                samples1 = sample_posterior(
                    grid1,
                    matrix1[index1],
                    n_samples,
                    rng,
                )
                displacement = samples1 - samples0

                output_arrays["x0_post_mean"][row_index] = float(
                    np.mean(samples0)
                )
                output_arrays["x1_post_mean"][row_index] = float(
                    np.mean(samples1)
                )
                output_arrays["x0_post_median"][row_index] = float(
                    np.median(samples0)
                )
                output_arrays["x1_post_median"][row_index] = float(
                    np.median(samples1)
                )
                output_arrays["dx_post_mean"][row_index] = float(
                    np.mean(displacement)
                )
                output_arrays["dx_post_median"][row_index] = float(
                    np.median(displacement)
                )
                output_arrays["dx_post_q025"][row_index] = float(
                    np.quantile(displacement, 0.025)
                )
                output_arrays["dx_post_q975"][row_index] = float(
                    np.quantile(displacement, 0.975)
                )
                output_arrays["dx_post_sd"][row_index] = float(
                    np.std(displacement, ddof=posterior_ddof)
                )
                output_arrays["p_dx_gt0"][row_index] = float(
                    np.mean(displacement > 0)
                )
                output_arrays["p_dx_lt0"][row_index] = float(
                    np.mean(displacement < 0)
                )
                success[row_index] = True

            except Exception as exc:
                messages[row_index] = (
                    f"transition_error: {type(exc).__name__}: {exc}"
                )

    for column in POSTERIOR_FLOAT_COLUMNS:
        table = table.append_column(
            column,
            pa.array(output_arrays[column], type=pa.float64()),
        )
    table = table.append_column(
        "posterior_success",
        pa.array(success, type=pa.bool_()),
    )
    table = table.append_column(
        "posterior_message",
        pa.array(messages, type=pa.string()),
    )
    return table


def enrich_pair_with_posterior(
    raw_path: Path,
    output_path: Path,
    posterior_dir: Path,
    subject: int,
    t0: int,
    t1: int,
    pattern: str,
    n_samples: int,
    seed: int,
    batch_size: int,
    compression: str,
    row_group_size: int,
) -> None:
    posterior0: Optional[Dict[str, Any]] = None
    posterior1: Optional[Dict[str, Any]] = None
    file_error: Optional[str] = None

    try:
        posterior0 = load_posterior_file(
            posterior_path(posterior_dir, subject, t0)
        )
        if t1 == t0:
            posterior1 = posterior0
        else:
            posterior1 = load_posterior_file(
                posterior_path(posterior_dir, subject, t1)
            )
    except Exception as exc:
        file_error = f"posterior_file_error: {type(exc).__name__}: {exc}"

    random_generator = np.random.default_rng(
        deterministic_pair_seed(seed, subject, t0, t1, pattern)
    )

    parquet_file = pq.ParquetFile(raw_path)
    writer: Optional[pq.ParquetWriter] = None
    processed_rows = 0

    try:
        for record_batch in parquet_file.iter_batches(
            batch_size=int(batch_size),
            use_threads=True,
        ):
            enriched = append_posterior_columns(
                record_batch=record_batch,
                posterior0=posterior0,
                posterior1=posterior1,
                file_error=file_error,
                n_samples=int(n_samples),
                rng=random_generator,
            )

            if writer is None:
                ensure_parent(output_path)
                writer = pq.ParquetWriter(
                    output_path,
                    enriched.schema,
                    compression=compression,
                    use_dictionary=True,
                    write_statistics=True,
                )

            writer.write_table(
                enriched,
                row_group_size=int(row_group_size),
            )
            processed_rows += int(enriched.num_rows)
            print(
                f"[INFO] Posterior block subject={subject} "
                f"{t0}->{t1} {pattern}: {processed_rows:,} rows"
            )

            del enriched, record_batch

    finally:
        if writer is not None:
            writer.close()
        del posterior0, posterior1
        gc.collect()


def copy_pair_without_posterior(
    raw_path: Path,
    output_path: Path,
    batch_size: int,
    compression: str,
    row_group_size: int,
) -> None:
    parquet_file = pq.ParquetFile(raw_path)
    writer: Optional[pq.ParquetWriter] = None

    try:
        for record_batch in parquet_file.iter_batches(
            batch_size=int(batch_size),
            use_threads=True,
        ):
            table = pa.Table.from_batches([record_batch])
            if writer is None:
                ensure_parent(output_path)
                writer = pq.ParquetWriter(
                    output_path,
                    table.schema,
                    compression=compression,
                    use_dictionary=True,
                    write_statistics=True,
                )
            writer.write_table(table, row_group_size=int(row_group_size))
            del table, record_batch
    finally:
        if writer is not None:
            writer.close()


# -----------------------------------------------------------------------------
# Pair files, resume and final concatenation
# -----------------------------------------------------------------------------


def write_raw_pair(
    transitions: pl.DataFrame,
    raw_path: Path,
    compression: str,
    row_group_size: int,
) -> int:
    if transitions.is_empty():
        return 0
    ensure_parent(raw_path)
    transitions.write_parquet(
        raw_path,
        compression=compression,
        statistics=True,
        row_group_size=int(row_group_size),
    )
    return int(transitions.height)


def finalize_pair_part(
    raw_path: Path,
    part_path: Path,
    args: argparse.Namespace,
    subject: int,
    t0: int,
    t1: int,
    pattern: str,
) -> int:
    if not raw_path.exists():
        return 0

    partial_path = part_path.with_suffix(part_path.suffix + ".partial")
    remove_if_exists(partial_path)

    if args.posterior_dir is not None:
        enrich_pair_with_posterior(
            raw_path=raw_path,
            output_path=partial_path,
            posterior_dir=args.posterior_dir,
            subject=subject,
            t0=t0,
            t1=t1,
            pattern=pattern,
            n_samples=int(args.n_posterior_samples),
            seed=int(args.posterior_seed),
            batch_size=int(args.posterior_chunk_size),
            compression=args.parquet_compression,
            row_group_size=int(args.parquet_row_group_size),
        )
    else:
        copy_pair_without_posterior(
            raw_path=raw_path,
            output_path=partial_path,
            batch_size=int(args.posterior_chunk_size),
            compression=args.parquet_compression,
            row_group_size=int(args.parquet_row_group_size),
        )

    if not partial_path.exists():
        raise RuntimeError(
            f"Pair output was not created for subject={subject}, {t0}->{t1}."
        )

    partial_path.replace(part_path)
    remove_if_exists(raw_path)
    return int(pq.ParquetFile(part_path).metadata.num_rows)


def empty_marker_path(part_path: Path) -> Path:
    return part_path.with_suffix(part_path.suffix + ".empty.json")


def write_empty_marker(
    marker_path: Path,
    subject: int,
    t0: int,
    t1: int,
    pattern: str,
) -> None:
    ensure_parent(marker_path)
    marker_path.write_text(
        json.dumps(
            {
                "subject": int(subject),
                "t0": int(t0),
                "t1": int(t1),
                "pattern": pattern,
                "rows": 0,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def pair_is_complete(part_path: Path, resume: bool) -> bool:
    if not resume:
        return False
    return part_path.exists() or empty_marker_path(part_path).exists()


def combine_parts(
    parts: Sequence[Path],
    output: Path,
    compression: str,
    batch_size: int,
    row_group_size: int,
) -> None:
    if not parts:
        raise RuntimeError("No non-empty transition parts were generated.")

    ensure_parent(output)
    partial_output = output.with_suffix(output.suffix + ".partial")
    remove_if_exists(partial_output)

    writer: Optional[pq.ParquetWriter] = None
    target_schema: Optional[pa.Schema] = None

    try:
        for part in parts:
            parquet_file = pq.ParquetFile(part)
            for record_batch in parquet_file.iter_batches(
                batch_size=int(batch_size),
                use_threads=True,
            ):
                table = pa.Table.from_batches([record_batch])

                if writer is None:
                    target_schema = table.schema
                    writer = pq.ParquetWriter(
                        partial_output,
                        target_schema,
                        compression=compression,
                        use_dictionary=True,
                        write_statistics=True,
                    )
                elif target_schema is not None and table.schema != target_schema:
                    table = table.cast(target_schema, safe=False)

                writer.write_table(
                    table,
                    row_group_size=int(row_group_size),
                )
                del table, record_batch

    finally:
        if writer is not None:
            writer.close()

    if not partial_output.exists():
        raise RuntimeError("Final Parquet output was not created.")
    partial_output.replace(output)


# -----------------------------------------------------------------------------
# Optional structural weights
# -----------------------------------------------------------------------------


def add_structural_weights_streaming(
    output: Path,
    args: argparse.Namespace,
) -> None:
    lazy_frame = pl.scan_parquet(output)

    clone_counts = lazy_frame.group_by(
        ["pattern", "subject", "aaSeqCDR3"]
    ).agg(pl.len().alias("_n_clone"))

    weighted = lazy_frame.join(
        clone_counts,
        on=["pattern", "subject", "aaSeqCDR3"],
        how="left",
    ).with_columns(
        (1.0 / pl.col("_n_clone").cast(pl.Float64)).alias("w_clone")
    )

    if args.balance_dt:
        dt_counts = lazy_frame.group_by(["pattern", "dt"]).agg(
            pl.len().alias("_n_dt")
        )
        weighted = weighted.join(
            dt_counts,
            on=["pattern", "dt"],
            how="left",
        ).with_columns(
            (1.0 / pl.col("_n_dt").cast(pl.Float64)).alias("w_dt")
        )
    else:
        weighted = weighted.with_columns(
            pl.lit(1.0, dtype=pl.Float64).alias("w_dt")
        )

    weighted = weighted.with_columns(
        pl.when(pl.col("obs_class") == "TT")
        .then(pl.lit(1.0))
        .when(pl.col("obs_class").is_in(["TF", "FT"]))
        .then(pl.lit(float(args.gamma_cross)))
        .when(pl.col("obs_class") == "FF")
        .then(pl.lit(float(args.gamma_ff)))
        .otherwise(pl.lit(1.0))
        .alias("w_class")
    )

    first_pass = output.with_suffix(".weights_first_pass.parquet")
    sink_parquet_streaming(
        weighted,
        first_pass,
        args.parquet_compression,
        int(args.parquet_row_group_size),
    )

    second_pass = pl.scan_parquet(first_pass)
    if args.normalize_weight_components:
        means_frame = collect_streaming(
            second_pass.select(
                [
                    pl.col("w_clone")
                    .filter(pl.col("w_clone") > 0)
                    .mean()
                    .alias("w_clone"),
                    pl.col("w_dt")
                    .filter(pl.col("w_dt") > 0)
                    .mean()
                    .alias("w_dt"),
                    pl.col("w_class")
                    .filter(pl.col("w_class") > 0)
                    .mean()
                    .alias("w_class"),
                ]
            )
        )
        means = means_frame.row(0, named=True)
        second_pass = second_pass.with_columns(
            [
                (pl.col("w_clone") / float(means["w_clone"])).alias(
                    "w_clone"
                ),
                (pl.col("w_dt") / float(means["w_dt"])).alias("w_dt"),
                (pl.col("w_class") / float(means["w_class"])).alias(
                    "w_class"
                ),
            ]
        )

    second_pass = second_pass.with_columns(
        (pl.col("w_clone") * pl.col("w_dt") * pl.col("w_class")).alias("w")
    )

    schema_names = set(second_pass.collect_schema().names())
    drop_columns = [
        column for column in ["_n_clone", "_n_dt"] if column in schema_names
    ]
    if drop_columns:
        second_pass = second_pass.drop(drop_columns)

    weighted_output = output.with_suffix(".weighted.partial.parquet")
    sink_parquet_streaming(
        second_pass,
        weighted_output,
        args.parquet_compression,
        int(args.parquet_row_group_size),
    )
    weighted_output.replace(output)
    remove_if_exists(first_pass)


# -----------------------------------------------------------------------------
# Reporting and work-directory manifest
# -----------------------------------------------------------------------------


def summarize_output(output: Path) -> Dict[str, Any]:
    lazy_frame = pl.scan_parquet(output)
    total_rows = int(pq.ParquetFile(output).metadata.num_rows)

    counts_by_dt_frame = collect_streaming(
        lazy_frame.group_by("dt").agg(pl.len().alias("n")).sort("dt")
    )
    counts_by_class_frame = collect_streaming(
        lazy_frame
        .group_by("obs_class")
        .agg(pl.len().alias("n"))
        .sort("obs_class")
    )

    summary: Dict[str, Any] = {
        "transition_rows": total_rows,
        "counts_by_dt": {
            str(row[0]): int(row[1])
            for row in counts_by_dt_frame.iter_rows()
        },
        "counts_by_obs_class": {
            str(row[0]): int(row[1])
            for row in counts_by_class_frame.iter_rows()
        },
    }

    schema_names = set(lazy_frame.collect_schema().names())
    if "posterior_success" in schema_names:
        posterior_summary = collect_streaming(
            lazy_frame.select(
                [
                    pl.col("posterior_success").sum().alias("n_success"),
                    pl.col("posterior_success").mean().alias("fraction"),
                ]
            )
        ).row(0, named=True)
        summary["posterior_success_n"] = int(
            posterior_summary["n_success"] or 0
        )
        summary["posterior_success_fraction"] = float(
            posterior_summary["fraction"] or 0.0
        )

    return summary


def write_report(path: Path, payload: Dict[str, Any]) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Latent transition build report\n\n")
        for key, raw_value in payload.items():
            value: Any = raw_value
            if isinstance(value, (dict, list, set, tuple)):
                value = json.dumps(value, ensure_ascii=False)
            handle.write(f"- {key}: {value}\n")


def build_manifest(
    args: argparse.Namespace,
    x_source: str,
) -> Dict[str, Any]:
    return {
        "input": input_signature(args.trajectories),
        "posterior_dir": (
            str(args.posterior_dir.resolve())
            if args.posterior_dir is not None
            else None
        ),
        "transition_mode": args.transition_mode,
        "min_dt": int(args.min_dt),
        "max_dt": int(args.max_dt) if args.max_dt is not None else None,
        "keep_obs_class": args.keep_obs_class,
        "n_posterior_samples": int(args.n_posterior_samples),
        "posterior_seed": int(args.posterior_seed),
        "posterior_chunk_size": int(args.posterior_chunk_size),
        "parquet_compression": args.parquet_compression,
        "parquet_row_group_size": int(args.parquet_row_group_size),
        "x_source": x_source,
        "script_architecture": "subject_pair_streaming_v3",
    }


def prepare_work_directory(
    work_directory: Path,
    manifest: Dict[str, Any],
    restart: bool,
) -> None:
    if restart and work_directory.exists():
        shutil.rmtree(work_directory)

    work_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = work_directory / "manifest.json"

    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous != manifest:
            raise RuntimeError(
                "The existing work directory was created with different input "
                "or parameters. Use --restart or choose another --work-dir."
            )
    else:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Memory-bounded construction and posterior propagation of latent clonotype transitions (ClonoDynamics Step 5)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--trajectories", type=Path, required=True)
    parser.add_argument("--posterior-dir", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--output-format",
        choices=["parquet", "csv", "feather"],
        default="parquet",
    )
    parser.add_argument("--parquet-compression", default="zstd")
    parser.add_argument("--parquet-row-group-size", type=int, default=100000)
    parser.add_argument("--file-sep", default=",")
    parser.add_argument("--write-csv-compat", action="store_true")

    parser.add_argument("--subject-col", default="subject")
    parser.add_argument("--clone-col", default="aaSeqCDR3")
    parser.add_argument("--time-col", default="time")
    parser.add_argument("--freq-col", default="freq")
    parser.add_argument("--log-freq-col", default="log_freq")
    parser.add_argument("--observable-col", default="observable")
    parser.add_argument("--min-freq", type=float, default=1e-15)
    parser.add_argument(
        "--allow-missing-state-detectability",
        action="store_true",
    )

    parser.add_argument(
        "--transition-mode",
        choices=["all", "adjacent", "adj", "both"],
        default="all",
    )
    parser.add_argument("--min-dt", type=int, default=1)
    parser.add_argument("--max-dt", type=int, default=6)
    parser.add_argument("--keep-obs-class", default=None)

    parser.add_argument("--n-posterior-samples", type=int, default=1000)
    parser.add_argument("--posterior-seed", type=int, default=123)
    parser.add_argument(
        "--posterior-chunk-size",
        type=int,
        default=5000,
        help="Arrow batch size used during posterior propagation",
    )

    parser.add_argument("--add-weights", action="store_true")
    parser.add_argument("--balance-dt", action="store_true")
    parser.add_argument("--gamma-cross", type=float, default=0.25)
    parser.add_argument("--gamma-ff", type=float, default=0.10)
    parser.add_argument(
        "--normalize-weight-components",
        action="store_true",
    )

    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no-resume", dest="resume", action="store_false")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--keep-work-dir", action="store_true")
    parser.add_argument("--report-md", type=Path, default=None)

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.output_format != "parquet":
        raise ValueError(
            "The memory-bounded builder writes Parquet. Use "
            "--output-format parquet. A CSV compatibility copy can be requested "
            "with --write-csv-compat."
        )
    if not args.trajectories.exists():
        raise FileNotFoundError(
            f"Trajectory file not found: {args.trajectories}"
        )
    if args.posterior_dir is not None and not args.posterior_dir.exists():
        raise FileNotFoundError(
            f"Posterior directory not found: {args.posterior_dir}"
        )
    if int(args.min_dt) < 1:
        raise ValueError("--min-dt must be at least 1")
    if args.max_dt is not None and int(args.max_dt) < int(args.min_dt):
        raise ValueError("--max-dt must be greater than or equal to --min-dt")
    if int(args.posterior_chunk_size) <= 0:
        raise ValueError("--posterior-chunk-size must be positive")
    if int(args.parquet_row_group_size) <= 0:
        raise ValueError("--parquet-row-group-size must be positive")
    if args.posterior_dir is not None and int(args.n_posterior_samples) <= 0:
        raise ValueError("--n-posterior-samples must be positive")


def main() -> None:
    args = build_arg_parser().parse_args()
    validate_arguments(args)

    print("[INFO] Step 5 latent transition construction — memory-bounded v3")

    input_scan = scan_table(args.trajectories, args.file_sep)
    input_schema = input_scan.collect_schema()
    canonical_expressions, x_source, canonical_numeric = (
        build_canonical_expressions(input_schema, args)
    )
    keep_classes = parse_keep_values(args.keep_obs_class)

    subjects_frame = collect_streaming(
        input_scan
        .select(
            pl.col(args.subject_col)
            .cast(pl.Int64, strict=False)
            .alias("subject")
        )
        .drop_nulls()
        .unique()
        .sort("subject")
    )
    subjects = [
        int(value)
        for value in subjects_frame.get_column("subject").to_list()
    ]

    input_rows = int(
        collect_streaming(input_scan.select(pl.len().alias("n"))).item()
    )
    print(f"[INFO] Input trajectory rows: {input_rows:,}")
    print(f"[INFO] Subjects: {subjects}")
    print(f"[INFO] x source: {x_source}")

    work_directory = args.work_dir or (
        args.out.parent / f".{args.out.stem}_step5_work"
    )
    manifest = build_manifest(args, x_source)
    prepare_work_directory(
        work_directory=work_directory,
        manifest=manifest,
        restart=bool(args.restart),
    )

    raw_directory = work_directory / "raw"
    parts_directory = work_directory / "parts"
    raw_directory.mkdir(parents=True, exist_ok=True)
    parts_directory.mkdir(parents=True, exist_ok=True)

    prepared_rows = 0
    expected_parts: List[Path] = []
    empty_pairs = 0

    for subject in subjects:
        print(f"[INFO] Preparing subject {subject}")
        subject_frame = prepare_subject(
            trajectories=args.trajectories,
            separator=args.file_sep,
            subject=subject,
            args=args,
            canonical_expressions=canonical_expressions,
            canonical_numeric=canonical_numeric,
        )

        prepared_rows += int(subject_frame.height)
        print(
            f"[INFO] Subject {subject} prepared rows: "
            f"{subject_frame.height:,}"
        )

        if subject_frame.is_empty():
            del subject_frame
            gc.collect()
            continue

        times = sorted(
            int(value)
            for value in subject_frame.get_column("time").unique().to_list()
        )

        if args.transition_mode in {"all", "both"}:
            for index, t0 in enumerate(times[:-1]):
                for t1 in times[index + 1 :]:
                    dt_value = int(t1 - t0)
                    if dt_value < int(args.min_dt):
                        continue
                    if args.max_dt is not None and dt_value > int(args.max_dt):
                        continue

                    stem = f"s{subject:04d}_t{t0:03d}_t{t1:03d}_all"
                    raw_path = raw_directory / f"{stem}.raw.parquet"
                    part_path = parts_directory / f"{stem}.parquet"

                    if pair_is_complete(part_path, bool(args.resume)):
                        if part_path.exists():
                            expected_parts.append(part_path)
                            print(
                                f"[INFO] Reusing completed part: {part_path.name}"
                            )
                        else:
                            empty_pairs += 1
                            print(
                                f"[INFO] Reusing empty-pair marker: "
                                f"{empty_marker_path(part_path).name}"
                            )
                        continue

                    print(
                        f"[INFO] Building subject={subject} pair={t0}->{t1} "
                        "pattern=all"
                    )
                    transitions = build_all_pair(
                        subject_frame=subject_frame,
                        subject=subject,
                        t0=t0,
                        t1=t1,
                        keep_classes=keep_classes,
                    )
                    marker = empty_marker_path(part_path)
                    remove_if_exists(marker)
                    remove_if_exists(raw_path)
                    raw_rows = write_raw_pair(
                        transitions=transitions,
                        raw_path=raw_path,
                        compression=args.parquet_compression,
                        row_group_size=int(args.parquet_row_group_size),
                    )
                    del transitions
                    gc.collect()
                    if raw_rows == 0:
                        write_empty_marker(
                            marker, subject, t0, t1, "all"
                        )
                        rows = 0
                    else:
                        rows = finalize_pair_part(
                            raw_path=raw_path,
                            part_path=part_path,
                            args=args,
                            subject=subject,
                            t0=t0,
                            t1=t1,
                            pattern="all",
                        )
                    if rows > 0:
                        expected_parts.append(part_path)
                    else:
                        empty_pairs += 1
                    print(f"[INFO] Completed part rows: {rows:,}")
                    gc.collect()

        if args.transition_mode in {"adjacent", "adj", "both"}:
            for t0, t1, transitions in build_adjacent_blocks(
                subject_frame=subject_frame,
                subject=subject,
                min_dt=int(args.min_dt),
                max_dt=args.max_dt,
                keep_classes=keep_classes,
            ):
                stem = f"s{subject:04d}_t{t0:03d}_t{t1:03d}_adj"
                raw_path = raw_directory / f"{stem}.raw.parquet"
                part_path = parts_directory / f"{stem}.parquet"

                if pair_is_complete(part_path, bool(args.resume)):
                    del transitions
                    if part_path.exists():
                        expected_parts.append(part_path)
                        print(
                            f"[INFO] Reusing completed part: {part_path.name}"
                        )
                    else:
                        empty_pairs += 1
                        print(
                            f"[INFO] Reusing empty-pair marker: "
                            f"{empty_marker_path(part_path).name}"
                        )
                    continue

                print(
                    f"[INFO] Building subject={subject} pair={t0}->{t1} "
                    "pattern=adj"
                )
                marker = empty_marker_path(part_path)
                remove_if_exists(marker)
                remove_if_exists(raw_path)
                raw_rows = write_raw_pair(
                    transitions=transitions,
                    raw_path=raw_path,
                    compression=args.parquet_compression,
                    row_group_size=int(args.parquet_row_group_size),
                )
                del transitions
                gc.collect()
                if raw_rows == 0:
                    write_empty_marker(
                        marker, subject, t0, t1, "adj"
                    )
                    rows = 0
                else:
                    rows = finalize_pair_part(
                        raw_path=raw_path,
                        part_path=part_path,
                        args=args,
                        subject=subject,
                        t0=t0,
                        t1=t1,
                        pattern="adj",
                    )
                if rows > 0:
                    expected_parts.append(part_path)
                else:
                    empty_pairs += 1
                print(f"[INFO] Completed part rows: {rows:,}")
                gc.collect()

        del subject_frame
        gc.collect()

    completed_parts = sorted(
        {part for part in expected_parts if part.exists()}
    )
    print(f"[INFO] Combining {len(completed_parts):,} completed pair parts")

    combine_parts(
        parts=completed_parts,
        output=args.out,
        compression=args.parquet_compression,
        batch_size=int(args.posterior_chunk_size),
        row_group_size=int(args.parquet_row_group_size),
    )

    if args.add_weights:
        print("[INFO] Adding structural weights in a streaming second pass")
        add_structural_weights_streaming(args.out, args)

    if args.write_csv_compat:
        csv_path = args.out.with_suffix(".csv")
        print(
            f"[WARN] CSV compatibility output may be extremely large: "
            f"{csv_path}"
        )
        csv_lazy = pl.scan_parquet(args.out)
        try:
            csv_lazy.sink_csv(csv_path)
        except (AttributeError, TypeError):
            collect_streaming(csv_lazy).write_csv(csv_path)

    output_summary = summarize_output(args.out)
    output_schema_names = set(
        pl.scan_parquet(args.out).collect_schema().names()
    )
    report_path = args.report_md or (
        args.out.parent / "latent_transition_build_report.md"
    )

    report_payload: Dict[str, Any] = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "architecture": "subject_pair_streaming_with_resumable_parts_v3",
        "trajectories": str(args.trajectories),
        "output": str(args.out),
        "trajectory_input_rows": int(input_rows),
        "prepared_trajectory_rows": int(prepared_rows),
        "subjects": len(subjects),
        "x_source": x_source,
        "transition_mode": args.transition_mode,
        "min_dt": int(args.min_dt),
        "max_dt": int(args.max_dt) if args.max_dt is not None else None,
        "keep_obs_class": args.keep_obs_class,
        "posterior_dir": (
            str(args.posterior_dir)
            if args.posterior_dir is not None
            else None
        ),
        "n_posterior_samples": int(args.n_posterior_samples),
        "posterior_chunk_size": int(args.posterior_chunk_size),
        "completed_parts": len(completed_parts),
        "empty_pairs": int(empty_pairs),
        "has_endpoint_intervals": {
            "x0_q025",
            "x0_q975",
            "x1_q025",
            "x1_q975",
        }.issubset(output_schema_names),
        "has_dx_sd": "dx_sd" in output_schema_names,
        "has_xstar_latent": "xstar_latent" in output_schema_names,
        "has_closure_component": "closure_component" in output_schema_names,
        "has_state_detectability": {
            "p_detect_state_t0",
            "p_detect_state_t1",
        }.issubset(output_schema_names),
        "has_full_posterior_propagation": (
            "dx_post_mean" in output_schema_names
        ),
    }
    report_payload.update(output_summary)
    write_report(report_path, report_payload)

    print(f"[INFO] Wrote transitions: {args.out}")
    print(f"[INFO] Wrote report: {report_path}")

    if not args.keep_work_dir:
        shutil.rmtree(work_directory, ignore_errors=True)

    print("[DONE]")


if __name__ == "__main__":
    main()