#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2-multirepresentation_clonotype_state_inference.py
=================================================

ClonoDynamics Step 2: pair-specific latent-state inference, replicate-resolved
observed measurements, and observation-state annotation from paired technical
replicate TCR repertoire sequencing data.

Purpose
-------
Step 2 converts each complete subject-timepoint technical-replicate pair into a
single clonotype-level state table while deliberately preserving the different
kinds of information carried by the two replicates.

The production architecture is replicate-resolved. The paired Negative-
Binomial model provides a noise-aware latent consensus representation and its
uncertainty, but the technical replicates are NOT collapsed into a single
measurement for every downstream estimand. Step 2 therefore exports three
conceptually distinct layers:

    1. model-based latent consensus representation;
    2. replicate-resolved observed measurements;
    3. observation-state annotations.

An auxiliary conditional-positive summary is retained for compatibility and
representation diagnostics, but it is not the primary fluctuation estimator in
the final ClonoDynamics pipeline.

Final downstream contract
-------------------------
The intended production roles are:

    latent consensus representation
        x_latent_median and the paired posterior
        -> noise-aware state representation
        -> conditioning geometry (including midpoint/xstar constructions)
        -> posterior uncertainty characterization
        -> latent-representation sensitivity analyses

    replicate-resolved observed measurements
        count_rep1/count_rep2
        freq_rep1/freq_rep2
        x_obs_rep1/x_obs_rep2
        -> primary replicate-decoupled forward estimator downstream:
           AB: condition on replicate A and measure displacement in B
           BA: condition on replicate B and measure displacement in A
        -> primary replicate-consistent fluctuation estimator downstream:
           Cov(delta x_A, delta x_B) on common4 eligibility

    observation-state annotations
        realized replicate detection
        model-based posterior-predictive detectability
        p_value-based operational observability
        -> operational TT/TF/FT/FF transition annotation
        -> observation-domain and alpha-threshold sensitivity analyses

The primary fluctuation conditioning coordinate is constructed downstream from
the latent consensus representation (xmid_latent). Step 2 itself does not
estimate forward drift, cross-replicate covariance, temporal scaling, or TT-
restricted dynamics.

Important interpretive boundary
--------------------------------
Step 2 does not define a single numerical abundance estimate as universally
optimal. In particular:

    - x_latent_median is a model-based consensus representation, not ground
      truth abundance;
    - x_obs_rep1/x_obs_rep2 preserve the independent observed measurements
      needed by replicate-decoupled estimands;
    - a zero count is a valid observation under the latent NB likelihood but is
      not interpreted as biological absence;
    - p_detect_state is posterior-predictive model-based detectability;
    - observable is a separate empirical operational annotation derived from
      p_value;
    - the Step-1 heavy-tail fit is descriptive and does not supply gamma or
      kappa to Step 2.

Position in the pipeline
------------------------

    Step 1  1-repertoire_characterization.py
        -> observed repertoire architecture and descriptive heavy-tail fitting
        -> no filtering and no parameter transfer to Step 2

    Step 2  2-multirepresentation_clonotype_state_inference.py
        -> pair-specific NB latent consensus + posterior uncertainty
        -> replicate-resolved observed measurements
        -> realized detection + model-based detectability
        -> p_value-based operational observability

    Step 3  3-latent_state_and_observation_model_diagnostics.py
        -> diagnostic characterization only; no filtering

    Step 4  replicate-resolved trajectory construction
    Step 5  generic transition construction

Later steps select estimand-specific representations. The primary forward
analysis uses replicate-decoupled observed measurements; the primary
fluctuation analysis uses cross-replicate displacement covariance with latent
conditioning.

Input
-----
The script expects repertoire files named by default as:

    <subject>_<time>-<replica>.<extension>

with replica 1 or 2, for example:

    1_1-1.tsv
    1_1-2.tsv
    1_2-1.tsv
    1_2-2.tsv

Supported input formats are TSV/text, CSV, Parquet and Feather/Arrow.
Only complete replicate pairs are fitted. No upstream PASS/FAIL filter is
applied in the primary workflow.

Each file must contain at least:

    aaSeqCDR3
    readCount

Within each technical replicate:

    - readCount is coerced to numeric;
    - non-finite and non-positive rows are removed;
    - duplicate aaSeqCDR3 rows are collapsed by summing readCount.

Thus the analysis unit is a unique amino-acid CDR3 clonotype within a
repertoire, and readCount denotes assigned sequencing reads, not cells, UMIs
or templates.

Technical-replicate pairing
---------------------------
For one subject-timepoint, the two repertoires are outer-joined by aaSeqCDR3.
Missing clonotypes in one replicate receive an observed count of zero:

    (c_i1 > 0, c_i2 = 0)
    (c_i1 = 0, c_i2 > 0).

A (0,0) row is not generated because the union is formed from clonotypes
observed in at least one replicate.

Replicate sequencing depths are:

    N_1 = sum_i c_i1
    N_2 = sum_i c_i2.

Replicate-specific observed relative frequencies are:

    f_i1^obs = c_i1 / N_1
    f_i2^obs = c_i2 / N_2.

For positive counts, the log-frequency measurements are:

    x_i1^obs = ln(f_i1^obs)
    x_i2^obs = ln(f_i2^obs).

They are exported as x_obs_rep1/x_obs_rep2, with explicit aliases
x_observed_rep1/x_observed_rep2 and f_observed_rep1/f_observed_rep2.
For a zero count the corresponding observed log-frequency is NaN, because a
non-detection is not assigned an arbitrary finite log abundance.

Layer 1: pair-specific latent consensus representation
------------------------------------------------------
For clonotype i and replicate r, the observation model is:

    c_ir | f_i, N_r, kappa
        ~ NegativeBinomial(mean = N_r f_i, size = kappa),

with:

    E[C | f]   = N_r f
    Var(C | f) = N_r f + (N_r f)^2 / kappa.

Both positive and zero counts enter the likelihood. Therefore a zero technical
count contributes evidence about the shared latent state under the model, but
does not imply f_i = 0.

The latent frequency is evaluated on a log-spaced grid of G points. With the
production defaults:

    G = 500
    f_max = 1
    f_min = 1 / max(mean(N_1,N_2), 1).

The prior implemented by noiseK_latent.py is normalized DISCRETE probability
mass on this grid:

    pi_j proportional to f_j^(-gamma),

without grid-cell-width quadrature weights. It should therefore be described
as a power-law-shaped discrete prior on a log-spaced frequency grid, not as an
exact continuous Pareto density integration.

For each complete technical-replicate pair, gamma and kappa are estimated
jointly by marginal likelihood. They are pair-specific and are not inherited
from Step 1.

For clonotype i the paired posterior is:

    P(f_j | c_i1,c_i2)
        proportional to
        P(c_i1 | f_j,N_1,kappa)
        P(c_i2 | f_j,N_2,kappa)
        pi_j.

Principal joint-posterior fields include:

    f_latent_mean, f_latent_median, f_latent_mode
    f_latent_q025, f_latent_q975
    x_latent_mean, x_latent_median, x_latent_mode
    x_latent_q025, x_latent_q975
    x_latent_sd
    posterior_entropy

where x = ln(f). The posterior median x_latent_median is the canonical
model-based consensus coordinate.

Single-replicate latent posteriors
----------------------------------
Using the same fitted pair-specific gamma, kappa and latent grid, the engine
also constructs:

    P(f | c_rep1)
    P(f | c_rep2).

Their summaries are exported as rep1_*/rep2_* and stable aliases including:

    x_rep1_latent_median
    x_rep2_latent_median.

These are retained for model diagnostics and latent-representation sensitivity
analyses. They are not the primary observed AB/BA forward measurements and do
not replace the joint consensus state.

The exact discrete-grid Bhattacharyya coefficient between the two
single-replicate posteriors is also exported for Step-3 diagnostics.

Layer 2: replicate-resolved observed measurements
-------------------------------------------------
The independent technical measurements are preserved explicitly through:

    count_rep1, count_rep2
    depth_rep1, depth_rep2
    freq_rep1, freq_rep2
    x_obs_rep1, x_obs_rep2
    present_rep1, present_rep2.

These fields are the source measurements used by the final replicate-resolved
dynamical estimands. Downstream, realized-positive eligibility is imposed by
the relevant estimand. In particular, common4 denotes clonotypes positive in
both replicates at both transition endpoints; common4 is constructed later and
is not a Step-2 filter.

Auxiliary conditional-positive summary
--------------------------------------
For backward compatibility and representation diagnostics, Step 2 also retains
an auxiliary positive-only point summary:

    f_positive_current
    x_positive_current.

If both technical replicates are positive:

    f_positive_current = sqrt(f_rep1 * f_rep2)
    x_positive_current = (x_obs_rep1 + x_obs_rep2)/2.

If only one replicate is positive, the positive replicate alone defines the
summary. The script records:

    n_positive_replicates
    positive_state_class       ONE_POSITIVE / TWO_POSITIVE
    positive_replicate_source  REP1_ONLY / REP2_ONLY / BOTH
    positive_abundance_method.

Historical aliases are retained:

    freq_geo == f_positive_current
    x_obs_geo == x_positive_current.

This positive-only summary is NOT the primary fluctuation displacement in the
final pipeline. The primary fluctuation estimand uses the covariance between
replicate-specific observed displacements; x_positive_current remains useful
for diagnostics, backward compatibility and representation comparisons.

Layer 3A: realized technical detection
--------------------------------------
The script records:

    present_rep1 = (count_rep1 > 0)
    present_rep2 = (count_rep2 > 0)
    present_both
    present_only_one.

These are realized sequencing outcomes. They must not be described as
biological presence/absence.

Layer 3B: model-based posterior-predictive detectability
--------------------------------------------------------
noiseK_latent.py returns replicate-specific and paired-state posterior-
predictive quantities:

    p_detect_rep1, p_detect_rep2
    p_dropout_rep1, p_dropout_rep2
    p_detect_state, p_dropout_state.

For the default detection threshold of one read, conditional on latent f:

    P(detect in replicate r | f)
        = 1 - P(C_r = 0 | f,N_r,kappa).

The paired state is considered technically detected if at least one replicate
reaches the count threshold. Its posterior-predictive dropout is:

    p_dropout_state = E_post[q_1(f) q_2(f)]
    p_detect_state  = 1 - p_dropout_state,

where q_r(f) is the conditional replicate dropout probability and the
expectation is under the shared paired posterior.

These are model-based quantities and are not independently calibrated
probabilities of biological presence.

Layer 3C: operational observability
-----------------------------------
For every clonotype, noiseK_latent.py supplies the fitted-model marginal log
probability logP of the observed count pair.

After successful pair fits are combined, empirical logP pools are constructed
according to --null:

    clone    pair-specific pool (historical option name)
    subject  subject-level pool, with global fallback for subjects with fewer
             than --min-pairs-per-subject fitted pairs
    global   one global pool.

With the production low-tail convention:

    p_value_i = empirical_CDF_pool(logP_i)
    observable_i = (p_value_i < alpha).

Operational observability is therefore distinct from both realized detection
and posterior-predictive detectability. Downstream TT/TF/FT/FF classes are
formed from observable unless a later sensitivity analysis explicitly
redefines alpha.

Re-aggregation across alpha
---------------------------
Changing alpha does not require refitting the latent NB model. In
--aggregate-only mode the script reuses the successful per-pair fit cache and
recomputes the empirical operational annotation for the requested alpha.

Pair-level parameter estimation
-------------------------------
The latent model is fitted independently for each complete subject-timepoint
technical-replicate pair using noiseK_latent.fit_noiseK_latent_powerlaw.
Production defaults include:

    gamma_init = 1.6
    k_init     = 50
    grid_size  = 500
    chunk_size = 20,000.

The engine uses L-BFGS-B with positive lower bounds of 1e-6 on gamma and
kappa.

Fit-success policy
------------------
Pair-level post-fit diagnostics are descriptive and never used as cohort-
relative filtering criteria. However, numerical optimizer failure is not a
valid inferred state. Therefore:

    - every complete technical-replicate pair is attempted;
    - only pairs with fit.success == True are written to the canonical
      per-clonotype fit cache and admitted to aggregate Step-2 outputs;
    - failed fits remain documented in latent_params_by_pair.csv and
      latent_pair_qc.csv;
    - aggregation is explicitly restricted to successful pair IDs from the
      current fit manifest/table, preventing stale files from a previous run
      from entering the current dataset.

No abundance, posterior-width, overlap, ONE_POSITIVE/TWO_POSITIVE or
posterior-predictive detectability threshold is used to exclude a successfully
fitted state.

Full posterior cache
--------------------
For each successful pair, the engine can store a compressed NPZ containing:

    f_grid
    posterior_matrix
    aaSeqCDR3.

These files support later full-posterior latent sensitivity calculations. They
do not define the primary replicate-resolved observed displacement.

Output structure
----------------
For <results-dir>/, principal outputs are:

1. per_clone_latent_logP/
   Successful pair-specific clonotype tables and full-posterior NPZ files.
   Historical directory/file naming is retained for compatibility.

2. latent_params_by_pair.csv
   One row per attempted complete pair, including optimizer success/failure,
   fitted parameters and pair-level diagnostics when available.

3. latent_pair_qc.csv
   Non-filtering post-fit diagnostic table. Numerical fit failure is documented
   as failure; cohort-relative warning flags do not exclude successful fits.

4. latent_params_aggregated.csv
   Descriptive aggregation of fitted pair-specific parameters.

5. latent_logP_cutoffs_<null>.csv
   Empirical operational-observability cutoffs.

6. per_clone_latent_<null>.<format>
   Main combined clonotype-state table containing latent, observed replicate,
   auxiliary positive-only and observation-state fields.

7. latent_trajectory_observability_<null>.csv
   Compact operational-observability table.

8. latent_trajectory_frequency_<null>.csv
   Compact downstream state table preserving the replicate-resolved observed
   measurements, latent consensus representation and observation annotations.

9. state_representation_manifest.json
   Machine-readable contract describing the final roles of the state layers,
   operational-observability definition and fit-success aggregation policy.

Latent count proxy
------------------
For count-scale descriptive comparisons only, the script retains:

    count_latent_median
        = f_latent_median * mean(N_1,N_2)

    x_count_latent_median
        = ln(count_latent_median + 1).

These are derived read-count-scale proxies, not independently fitted latent
cell or template counts.

No upstream QC filtering
------------------------
The deprecated command-line options --use-qc and --qc-csv are accepted only
for compatibility and are ignored. Step 2 always attempts every complete
technical-replicate pair found by the filename pattern.

Recommended terminology
-----------------------
Use:

    model-based latent consensus representation
        for f_latent_* / x_latent_*;

    replicate-resolved observed measurement
        for count_rep*, freq_rep*, x_obs_rep*;

    auxiliary conditional-positive summary
        for f_positive_current / x_positive_current;

    realized technical detection
        for present_rep1 / present_rep2;

    model-based posterior-predictive detectability
        for p_detect_* / p_dropout_*;

    operational observability
        for p_value / observable and downstream TT/TF/FT/FF annotation.

Avoid:

    true abundance;
    biological presence/absence for a zero technical count;
    independently calibrated detection probability for p_detect_state;
    describing x_positive_current as the primary fluctuation estimator;
    implying that Step-1 heavy-tail parameters feed the Step-2 prior.

Typical usage
-------------
Fit and aggregate the reference alpha=0.05 state table:

    python3 2-multirepresentation_clonotype_state_inference.py \
        --data-dir ./dataset_longitudinal \
        --results-dir ./results/2-clonotype_state_inference/p_05 \
        --alpha 0.05 \
        --null subject \
        --tail low \
        --n-jobs 4 \
        --intermediate-format parquet \
        --parquet-compression snappy

Re-aggregate the same successful fit cache at another operational alpha without
refitting:

    python3 2-multirepresentation_clonotype_state_inference.py \
        --aggregate-only \
        --fit-cache-dir ./results/2-clonotype_state_inference/p_05 \
        --aggregate-results-dir ./results/2-clonotype_state_inference/p_025 \
        --alpha 0.025 \
        --null subject \
        --tail low

Conceptual summary
------------------
Step 2 is a representation layer, not a dynamical inference step. The paired
NB model supplies a noise-aware consensus coordinate, while the original
replicate measurements are deliberately preserved so their independence can
be used downstream to identify forward and fluctuation estimands. Observation
state is carried separately from both abundance representations.
"""

from __future__ import annotations

import argparse
import gc
import json
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
    r"(?:\.(?:tsv|txt|csv|parquet|pq|feather|arrow))?$"
)

PARQUET_EXTS = {".parquet", ".pq"}
FEATHER_EXTS = {".feather", ".arrow"}
CSV_EXTS = {".csv", ".tsv", ".txt"}

SCRIPT_VERSION = "v9-replicate-resolved-state-contract-2026-09-09"


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

    # Remove same-pair cache artifacts from any earlier run before fitting.
    # This prevents a stale file in a different intermediate format from being
    # preferred later by cache discovery if the current run changes format or
    # if the current optimizer fit fails. Other pair IDs are protected at
    # aggregation time by successful-pair filtering.
    cache_stem = Path(job["perclone_dir"]) / f"latent_per_clone_logP_{pair_id}"
    for suffix in (".parquet", ".feather", ".csv", ".posterior.npz"):
        candidate = cache_stem.with_suffix(suffix) if suffix != ".posterior.npz" else Path(str(cache_stem) + suffix)
        try:
            if candidate.exists():
                candidate.unlink()
        except OSError:
            pass

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
        # -----------------------------------------------------------------
        # Auxiliary conditional-positive abundance summary.
        #
        # Historical aliases `freq_geo` and `x_obs_geo` are retained exactly
        # for backward compatibility and representation diagnostics. In the
        # final pipeline this summary is not the primary fluctuation
        # displacement; replicate-resolved observed measurements are preserved
        # separately for the primary downstream estimands.
        # -----------------------------------------------------------------
        merged["freq_geo"] = freq_geo
        merged["x_obs_rep1"] = np.where(c1 > 0, np.log(np.maximum(f1, 1e-300)), np.nan)
        merged["x_obs_rep2"] = np.where(c2 > 0, np.log(np.maximum(f2, 1e-300)), np.nan)
        merged["x_obs_geo"] = np.log(np.maximum(freq_geo, 1e-300))

        # Stable explicit aliases for the replicate-resolved observed channel.
        # These are numerically identical to freq_rep*/x_obs_rep* and make the
        # final replicate-resolved contract self-documenting downstream.
        merged["f_observed_rep1"] = merged["freq_rep1"].astype(float)
        merged["f_observed_rep2"] = merged["freq_rep2"].astype(float)
        merged["x_observed_rep1"] = merged["x_obs_rep1"].astype(float)
        merged["x_observed_rep2"] = merged["x_obs_rep2"].astype(float)

        merged["present_rep1"] = c1 > 0
        merged["present_rep2"] = c2 > 0
        merged["present_both"] = both
        merged["present_only_one"] = only1 | only2

        n_positive = (c1 > 0).astype(np.int8) + (c2 > 0).astype(np.int8)
        merged["n_positive_replicates"] = n_positive
        merged["positive_state_class"] = np.where(
            n_positive == 2, "TWO_POSITIVE", "ONE_POSITIVE"
        )
        merged["positive_replicate_source"] = np.select(
            [both, only1, only2],
            ["BOTH", "REP1_ONLY", "REP2_ONLY"],
            default="UNDEFINED",
        )
        merged["positive_abundance_method"] = np.select(
            [both, only1, only2],
            [
                "GEOMETRIC_MEAN_TWO_POSITIVE",
                "SINGLE_POSITIVE_REP1",
                "SINGLE_POSITIVE_REP2",
            ],
            default="UNDEFINED",
        )
        merged["f_positive_current"] = merged["freq_geo"].astype(float)
        merged["x_positive_current"] = merged["x_obs_geo"].astype(float)

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
        fit_success = bool(getattr(fit, "success", False))

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

        # Canonical per-clone files are written only for numerically successful
        # optimizer fits. noiseK_latent may still return posterior summaries
        # when optimizer.success is False; those are retained only for failure
        # diagnostics in the pair-level table and are not admitted downstream.
        stem = f"latent_per_clone_logP_{pair_id}"
        out_fp = None
        if fit_success:
            out_fp = write_table(
                per_clone,
                Path(job["perclone_dir"]) / stem,
                fmt=job["intermediate_format"],
                compression=job["parquet_compression"],
                index=False,
            )
            if bool(job.get("write_csv_compat", False)) and job["intermediate_format"] != "csv":
                per_clone.to_csv((Path(job["perclone_dir"]) / stem).with_suffix(".csv"), index=False)
        else:
            # The engine may have written the requested full-posterior cache
            # before exposing optimizer.success. Remove it so a failed fit can
            # never be mistaken for a canonical successful cache entry.
            try:
                if posterior_npz_path.exists():
                    posterior_npz_path.unlink()
            except OSError:
                pass

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
        n_one_positive = int((pd.to_numeric(per_clone.get("n_positive_replicates"), errors="coerce") == 1).sum()) if len(per_clone) else 0
        n_two_positive = int((pd.to_numeric(per_clone.get("n_positive_replicates"), errors="coerce") == 2).sum()) if len(per_clone) else 0
        fraction_one_positive = float(n_one_positive / len(per_clone)) if len(per_clone) else float("nan")
        fraction_two_positive = float(n_two_positive / len(per_clone)) if len(per_clone) else float("nan")
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
            "n_one_positive": n_one_positive,
            "n_two_positive": n_two_positive,
            "fraction_one_positive": fraction_one_positive,
            "fraction_two_positive": fraction_two_positive,
            "fraction_joint_posterior_wide": fraction_joint_posterior_wide,
            "success": fit_success,
            "message": str(getattr(fit, "message", "")),
            "log_likelihood": float(getattr(fit, "log_likelihood", np.nan)),
            "QC_pass": True if allowed_pairs_was_used else np.nan,
            "QC_used_for_filtering": allowed_pairs_was_used,
            "perclone_file": str(out_fp) if out_fp is not None else "",
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
        "fraction_present_only_one", "n_one_positive", "n_two_positive",
        "fraction_one_positive", "fraction_two_positive",
        "fraction_joint_posterior_wide",
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


def load_perclone_latent(
    perclone_dir: Path,
    file_sep: str = "\t",
    successful_pair_ids: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """Load canonical per-clone fit files, optionally restricted to successful pairs.

    Restricting by successful_pair_ids protects aggregation from stale cache
    files left by an earlier run and guarantees that numerical fit failures do
    not enter the current combined state table.
    """
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

    out = pd.concat(rows, ignore_index=True)
    if successful_pair_ids is not None:
        allowed = {str(x) for x in successful_pair_ids}
        if not allowed:
            raise ValueError("No successful pair IDs were supplied for aggregation.")
        if "pair_id" not in out.columns:
            raise ValueError("Per-clone cache lacks pair_id required for successful-fit filtering.")
        out = out[out["pair_id"].astype(str).isin(allowed)].copy()
        found = set(out["pair_id"].astype(str).unique())
        missing = sorted(allowed - found)
        if missing:
            raise FileNotFoundError(
                "Successful fit pair(s) missing from per-clone cache: " + ", ".join(missing)
            )
        if out.empty:
            raise ValueError("Successful-fit filtering removed all per-clone rows.")
    return out


def successful_pair_ids_from_params(params_csv: Path) -> Set[str]:
    """Return pair IDs whose numerical optimizer fit succeeded."""
    params = pd.read_csv(params_csv)
    if params.empty or "pair_id" not in params.columns or "success" not in params.columns:
        raise ValueError(
            f"Cannot derive successful pair IDs from {params_csv}; expected pair_id and success columns."
        )
    success = parse_bool_series(params["success"])
    ids = set(params.loc[success, "pair_id"].astype(str))
    if not ids:
        raise ValueError(f"No successful latent fits recorded in {params_csv}")
    return ids


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
        params = params[parse_bool_series(params["success"])].copy()
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
    """Write compact downstream tables while preserving the final state layers.

    Historical filenames are retained for compatibility. The compact state
    table preserves the latent consensus representation, replicate-resolved
    observed measurements, auxiliary positive-only summary, realized detection,
    model-based detectability and operational observability.
    """
    results_dir = Path(results_dir)

    # Per clone/time operational-observability summary.
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

    # Canonically there is one row per subject/clone/time, so most aggregation
    # operations are identities. Rules are nevertheless explicit so duplicated
    # rows cannot silently coerce string/bool state descriptors through numeric
    # aggregation.
    agg_dict = {}

    numeric_median_fields = [
        # observed/positive abundance
        "freq_rep1", "freq_rep2", "freq_geo", "x_obs_rep1", "x_obs_rep2", "x_obs_geo",
        "f_observed_rep1", "f_observed_rep2", "x_observed_rep1", "x_observed_rep2",
        "f_positive_current", "x_positive_current", "n_positive_replicates",
        # joint latent state
        "f_latent_mean", "f_latent_median", "f_latent_mode", "f_latent_q025", "f_latent_q975",
        "x_latent_mean", "x_latent_median", "x_latent_mode", "x_latent_q025", "x_latent_q975", "x_latent_sd",
        "posterior_entropy",
        # model-based detectability
        "p_detect_state", "p_dropout_state",
        "p_detect_rep1", "p_detect_rep2", "p_dropout_rep1", "p_dropout_rep2",
        # model compatibility / operational threshold statistic
        "logP", "p_value", "p_emp_low", "p_emp_high",
    ]
    first_fields = [
        "count_rep1", "count_rep2", "depth_rep1", "depth_rep2", "depth_mean",
        "present_rep1", "present_rep2", "present_both", "present_only_one",
        "positive_state_class", "positive_replicate_source", "positive_abundance_method",
    ]

    # Preserve single-replicate latent summaries and stable aliases when present.
    for c in df.columns:
        if c.startswith("rep1_") or c.startswith("rep2_") or c.startswith("f_rep1_") or c.startswith("f_rep2_") or c.startswith("x_rep1_") or c.startswith("x_rep2_"):
            if c not in numeric_median_fields:
                numeric_median_fields.append(c)

    for c in numeric_median_fields:
        if c in df.columns:
            agg_dict[c] = (c, "median")
    for c in first_fields:
        if c in df.columns:
            agg_dict[c] = (c, "first")
    if "observable" in df.columns:
        agg_dict["observable"] = ("observable", "max")

    freq = df.groupby(["subject", clone_col, "time"], as_index=False).agg(**agg_dict)
    out_freq = results_dir / f"latent_trajectory_frequency_{null_choice}.csv"
    freq.to_csv(out_freq, index=False)
    return out_obs, out_freq



def write_state_representation_manifest(
    results_dir: Path,
    null_choice: str,
    alpha: float,
    tail: str,
    n_successful_pairs: Optional[int] = None,
) -> Path:
    """Write the production contract for Step-2 state representations."""
    payload = {
        "script": "2-multirepresentation_clonotype_state_inference.py",
        "script_version": SCRIPT_VERSION,
        "design": "replicate_resolved_state_representation_layer",
        "step_role": (
            "Pair-specific latent consensus inference plus preservation of independent observed replicate "
            "measurements and separate observation-state annotations. Dynamic estimands are defined downstream."
        ),
        "layers": {
            "latent_consensus": {
                "point_state": "x_latent_median",
                "frequency_state": "f_latent_median",
                "zero_count_enters_likelihood": True,
                "model": "paired Negative-Binomial likelihood + discrete power-law-shaped prior",
                "prior_parameter_source": "pair-specific Step-2 marginal-likelihood fit; not Step 1",
                "intended_roles": [
                    "noise_aware_state_representation",
                    "conditioning_geometry",
                    "xmid_xstar_construction_downstream",
                    "posterior_uncertainty_characterization",
                    "latent_representation_sensitivity",
                ],
                "is_ground_truth": False,
            },
            "replicate_resolved_observed": {
                "counts": ["count_rep1", "count_rep2"],
                "frequencies": ["freq_rep1", "freq_rep2"],
                "log_frequencies": ["x_obs_rep1", "x_obs_rep2"],
                "explicit_aliases": {
                    "f_observed_rep1": "freq_rep1",
                    "f_observed_rep2": "freq_rep2",
                    "x_observed_rep1": "x_obs_rep1",
                    "x_observed_rep2": "x_obs_rep2",
                },
                "zero_observed_log_frequency": "NaN/non-detection; no arbitrary finite log abundance assigned",
                "primary_downstream_roles": [
                    "replicate_decoupled_observed_AB_BA_forward",
                    "cross_replicate_displacement_covariance_on_common4",
                ],
            },
            "auxiliary_conditional_positive": {
                "point_state": "x_positive_current",
                "frequency_state": "f_positive_current",
                "legacy_aliases": {
                    "x_positive_current": "x_obs_geo",
                    "f_positive_current": "freq_geo",
                },
                "two_positive_rule": "geometric mean of positive replicate relative frequencies",
                "one_positive_rule": "relative frequency of the single positive replicate",
                "final_pipeline_status": "auxiliary/backward-compatible representation; not primary fluctuation displacement",
                "probabilistic_uncertainty_calibrated": False,
            },
            "observation": {
                "realized_detection": [
                    "present_rep1", "present_rep2", "present_both", "present_only_one",
                    "n_positive_replicates", "positive_state_class",
                ],
                "model_based_posterior_predictive_detectability": [
                    "p_detect_state", "p_dropout_state",
                    "p_detect_rep1", "p_detect_rep2",
                    "p_dropout_rep1", "p_dropout_rep2",
                ],
                "operational_observability": ["logP", "p_value", "observable"],
                "detectability_independently_calibrated": False,
                "operational_state_is_biological_presence_absence": False,
            },
        },
        "primary_estimands_defined_downstream": {
            "forward": "observed replicate-decoupled AB/BA",
            "fluctuations": "Cov(delta_x_rep1, delta_x_rep2) on common4",
            "fluctuation_conditioning": "xmid_latent",
        },
        "operational_observability": {
            "null_pool": str(null_choice),
            "alpha": float(alpha),
            "tail": str(tail),
            "definition": "observable = (empirical p_value from fitted-model logP pool < alpha)",
        },
        "fit_success_policy": {
            "complete_pairs_attempted_without_upstream_qc_filter": True,
            "aggregate_successful_optimizer_fits_only": True,
            "n_successful_pairs": int(n_successful_pairs) if n_successful_pairs is not None else None,
            "stale_cache_protection": "aggregation restricted by successful pair IDs when a fit-parameter table is available",
            "postfit_cohort_relative_qc_filters_data": False,
        },
        "step1_relationship": {
            "step1_is_descriptive": True,
            "step1_gamma_passed_to_step2": False,
            "step1_parameterizes_step2_prior": False,
            "step1_validates_step2_prior": False,
        },
        "compatibility": {
            "historical_output_names_retained": True,
            "noiseK_latent_engine_changed_in_this_step_revision": False,
        },
    }
    path = Path(results_dir) / "state_representation_manifest.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path

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
    successful_pair_ids: Optional[Set[str]] = None,
) -> Path:
    df = load_perclone_latent(
        perclone_dir,
        file_sep=file_sep,
        successful_pair_ids=successful_pair_ids,
    )
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
    state_manifest = write_state_representation_manifest(
        Path(results_dir),
        null_choice=null_choice,
        alpha=alpha,
        tail=tail,
        n_successful_pairs=(len(successful_pair_ids) if successful_pair_ids is not None else None),
    )

    print("Saved aggregate multirepresentation state outputs:")
    print(" -", main_out)
    print(" -", cutoffs_out)
    print(" -", out_obs)
    print(" -", out_freq)
    print(" -", state_manifest)
    return main_out


# =============================================================================
# CLI
# =============================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pair-specific latent consensus inference, replicate-resolved observed states and observation annotation (ClonoDynamics Step 2)."
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

    if args.skip_fit:
        raise ValueError(
            "--skip-fit was a legacy unused option and is no longer silently accepted. "
            "Use --aggregate-only with --fit-cache-dir and --aggregate-results-dir."
        )

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

        params_csv = Path(args.fit_cache_dir) / "latent_params_by_pair.csv"
        successful_pair_ids = None
        if params_csv.exists():
            successful_pair_ids = successful_pair_ids_from_params(params_csv)
            print(
                f"[aggregate-only] Restricting cache to {len(successful_pair_ids)} successful pair fits "
                f"from {params_csv}"
            )
        else:
            print(
                "[warning] latent_params_by_pair.csv not found in fit cache; "
                "aggregate-only cannot verify successful-fit pair membership."
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
            successful_pair_ids=successful_pair_ids,
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

    successful_pair_ids = successful_pair_ids_from_params(params_csv)
    print(f"[aggregate] Restricting to {len(successful_pair_ids)} successful pair fits from current run.")

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
        successful_pair_ids=successful_pair_ids,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
