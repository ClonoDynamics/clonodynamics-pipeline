#!/usr/bin/env python3
"""
noiseK_latent.py
================

Latent-frequency inference engine for the paired technical-replicate
observation model used by ClonoDynamics.

Overview
--------
This module is a reusable statistical dependency of:

    2-multirepresentation_clonotype_state_inference.py

and is not a numbered pipeline step by itself.

Its role is deliberately limited to inference under the technical-replicate
observation model. For each biological sampling time, Step 2 supplies paired
technical-replicate clonotype counts and replicate sequencing depths. This
module then fits a pair-specific Negative-Binomial observation model with a
discrete power-law-shaped prior on latent clonotype relative frequency and
reconstructs clonotype-specific posterior distributions.

In the final ClonoDynamics architecture, the joint paired-replicate posterior
is a model-based latent consensus representation. It is useful for
noise-aware state representation, conditioning geometry, posterior uncertainty
and latent sensitivity analyses, but it is not treated as ground-truth
clonotype abundance and is not the sole input to downstream dynamics.

The final primary dynamical estimands are defined outside this module:

    - primary forward drift uses replicate-decoupled OBSERVED single-measure
      AB/BA estimators;

    - primary fluctuations use the cross-replicate covariance of OBSERVED
      log-frequency displacements on common4, with xmid_latent used as the
      primary conditioning coordinate.

Accordingly, the single-replicate latent posteriors produced here are retained
for model diagnostics, latent-state sensitivities and optional posterior
propagation. They should not be described as the primary observed AB/BA
forward estimator or as the primary fluctuation estimator.

The module provides:

    - pair-specific estimation of the prior exponent gamma and the
      Negative-Binomial size/overdispersion parameter kappa;

    - joint paired-replicate latent-frequency posteriors;

    - single-replicate latent posteriors evaluated with the same fitted
      pair-specific gamma, kappa and frequency grid;

    - exact discrete-grid Bhattacharyya overlap between the two
      single-replicate posteriors;

    - posterior-predictive replicate-specific and paired-state detectability;

    - optional storage of the full joint posterior matrix;

    - utility functions for Monte-Carlo propagation of independent latent
      endpoint posteriors into log-frequency displacement summaries.

The module does NOT impose a longitudinal dynamical model, does NOT define the
operational observation classes TT/TF/FT/FF, and does NOT perform biological
presence/absence inference.

Position in the final ClonoDynamics architecture
------------------------------------------------
The relevant early-pipeline relationship is:

    Step 1
    1-repertoire_characterization.py
        -> descriptive observed repertoire architecture
        -> descriptive heavy-tail characterization
        -> does NOT pass its fitted tail exponent to Step 2

    Step 2
    2-multirepresentation_clonotype_state_inference.py
        -> reads and pairs technical-replicate repertoires
        -> outer-joins clonotypes and represents replicate absence as count 0
        -> calls this module independently for each complete replicate pair
        -> retains the joint latent consensus representation
        -> retains replicate-resolved observed measurements
        -> constructs realized-detection and operational-observability fields

    noiseK_latent.py
        -> fits pair-specific gamma and kappa
        -> constructs joint and single-replicate latent posteriors
        -> returns model-based posterior uncertainty and detectability
        -> returns fitted-model marginal logP

    downstream dynamics
        -> observed replicate-decoupled AB/BA forward estimator
        -> observed cross-replicate displacement covariance on common4
        -> latent coordinates retained for conditioning and sensitivity

The heavy-tail analysis of Step 1 motivates a noise-aware representation but
does not numerically parameterize this module. gamma and kappa are estimated
anew for each complete technical-replicate pair.

Observation model
-----------------
For clonotype i measured in two technical replicates:

    c_i1 | f_i, N_1, kappa
        ~ NegativeBinomial(mean = f_i * N_1, size = kappa)

    c_i2 | f_i, N_2, kappa
        ~ NegativeBinomial(mean = f_i * N_2, size = kappa).

Here:

    c_ir
        observed read count in replicate r;

    N_r
        total sequencing depth of replicate r;

    f_i
        latent clonotype relative frequency shared by the two technical
        replicates under the fitted observation model;

    kappa
        pair-specific Negative-Binomial size/overdispersion parameter.

The implementation uses the mean/size parameterization:

    mu = f * N

    Var(C | f)
        = mu + mu^2 / kappa.

For observed count c, mean mu and size r = kappa:

    P(C=c | mu,r)
        =
        Gamma(c+r)
        -------------------------
        Gamma(r) Gamma(c+1)

        * [r/(r+mu)]^r

        * [mu/(r+mu)]^c.

Zero counts
-----------
Observed zero counts are valid Negative-Binomial observations.

Thus:

    c = 0

does NOT imply:

    f = 0.

A clonotype observed in one replicate and absent in the other remains a valid
paired observation. The zero contributes information to the joint latent
likelihood according to replicate depth and the fitted Negative-Binomial
model.

This is a model statement, not a biological absence statement.

Latent-frequency grid
---------------------
The posterior is evaluated on a log-spaced grid:

    f_1, ..., f_G

with default:

    G = 500.

The upper bound is:

    f_max = 1.

Unless supplied explicitly, the lower bound is:

                       1
    f_min = -------------------------
             mean(N_1, N_2).

More exactly, the implementation uses:

    f_min = 1 / max(mean(depths), 1).

The grid is then:

    logspace(log10(f_min), log10(1), G).

Power-law-shaped prior: exact implementation
--------------------------------------------
For grid point f_j, the implementation assigns normalized DISCRETE prior mass:

    pi_j
        proportional to f_j^(-gamma)

with:

    sum_j pi_j = 1.

Equivalently:

    log pi_j
        = -gamma * log(f_j) - log Z.

This distinction is important. The code treats f_j^(-gamma) as probability
MASS on the log-spaced grid. It does not multiply by grid-cell widths as a
quadrature approximation to a continuous density p(f) df.

Manuscript wording should therefore describe this as a discrete
power-law-shaped prior on a log-spaced latent-frequency grid unless the
numerical integration scheme is changed explicitly.

Joint posterior
---------------
For one clonotype with observed count pair:

    (c_1, c_2),

the joint latent posterior is:

    P(f_j | c_1,c_2)
        proportional to
        P(c_1 | f_j,N_1,kappa)
        P(c_2 | f_j,N_2,kappa)
        pi_j.

Normalization is performed with log-sum-exp.

Pair-level parameter estimation
-------------------------------
The public fitting function:

    fit_noiseK_latent_powerlaw(...)

estimates:

    gamma
    kappa

for one paired technical-replicate sample.

For clonotype i, its marginal count-pair probability on the discrete grid is:

    P(c_i1,c_i2 | gamma,kappa)
        =
        sum_j
            P(c_i1 | f_j,N_1,kappa)
            P(c_i2 | f_j,N_2,kappa)
            pi_j(gamma).

The pair-level objective is:

    log L(gamma,kappa)
        =
        sum_i
            log P(c_i1,c_i2 | gamma,kappa).

The implementation minimizes:

    -log L

using:

    scipy.optimize.minimize
    method = L-BFGS-B

with numerical lower bounds:

    gamma >= 1e-6
    kappa >= 1e-6.

Default starting values are:

    gamma_init = 1.6
    k_init     = 50.

Chunked likelihood evaluation
-----------------------------
The marginal likelihood is evaluated in clonotype chunks to limit peak memory.

Default:

    chunk_size = 20,000.

Chunking changes only the computation strategy, not the fitted objective.

Per-clonotype marginal log probability
--------------------------------------
After fitting gamma and kappa, the output field:

    logP

is:

    logP_i
        =
        log sum_j
            P(c_i1 | f_j,N_1,kappa_hat)
            P(c_i2 | f_j,N_2,kappa_hat)
            pi_j(gamma_hat).

Thus logP is the fitted-model marginal log probability of the observed
technical-replicate count pair under the discrete latent-frequency model.

This module does not convert logP into an operational observation state.
Step 2 constructs empirical logP pools, p_value and observable(alpha)
separately. Operational observability is therefore distinct from the
posterior-predictive detectability quantities calculated here.

Joint posterior summaries
-------------------------
For every clonotype the paired-replicate posterior provides:

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
    posterior_entropy

where:

    x = log(f).

Posterior means are probability-weighted discrete-grid means.

Posterior mode is the maximum-posterior grid point.

Median and quantiles are obtained by interpolation against the cumulative
discrete posterior on the f or log(f) grid.

Posterior entropy is:

    H
        = -sum_j p_j log(p_j).

Consensus latent representation
-------------------------------
The JOINT posterior:

    P(f | c_rep1,c_rep2)

provides the model-based latent consensus representation of one biological
timepoint.

The posterior median x_latent_median is the principal point summary retained by
Step 2 for latent-state representation and by downstream analyses for
conditioning geometry where specified.

This representation is not interpreted as an experimentally observed true
abundance. In particular, ONE_POSITIVE states can be shifted relative to the
positive replicate because the zero count in the other replicate contributes
to the joint likelihood.

Single-replicate latent posteriors
----------------------------------
Using the SAME fitted:

    gamma_hat
    kappa_hat
    f_grid,

the module also calculates:

    P(f | c_rep1)

and:

    P(f | c_rep2).

These yield:

    rep1_f_latent_mean
    rep1_f_latent_median
    rep1_f_latent_mode
    rep1_f_latent_q025
    rep1_f_latent_q975

    rep1_x_latent_mean
    rep1_x_latent_median
    rep1_x_latent_mode
    rep1_x_latent_q025
    rep1_x_latent_q975
    rep1_x_latent_sd
    rep1_posterior_entropy

and the corresponding rep2_* fields.

In the final pipeline these are retained primarily for:

    - technical diagnostics of replicate-specific latent behavior;
    - latent-state replicate agreement;
    - sensitivity analyses that use latent replicate-specific quantities;
    - optional posterior-based displacement propagation.

They do not replace the replicate-resolved observed measurements used by the
primary AB/BA forward estimator or by the primary cross-replicate fluctuation
estimator.

Stable replicate aliases
------------------------
For compatibility with downstream analysis code, the module also exposes:

    f_rep1_latent_median
    f_rep2_latent_median

    x_rep1_latent_median
    x_rep2_latent_median

    x_rep1_latent_q025
    x_rep1_latent_q975

    x_rep2_latent_q025
    x_rep2_latent_q975.

Replicate-posterior Bhattacharyya coefficient
---------------------------------------------
Because the two single-replicate posteriors are evaluated on the same discrete
frequency grid, their exact grid-level Bhattacharyya coefficient is:

    BC
        = sum_j sqrt[p_1(j) p_2(j)].

The output is:

    replicate_posterior_bhattacharyya.

The result is numerically clipped to [0,1]. This is the exact coefficient for
the two normalized discrete posteriors constructed by this module.

Replicate latent-difference diagnostics
---------------------------------------
The module provides convenience fields:

    dx_rep_latent_mean
    dx_rep_latent_median
    dx_rep_latent_mode
    dx_rep_latent_q025
    dx_rep_latent_q975
    dx_rep_latent_sd_approx.

Important distinction:

    dx_rep_latent_mean

is the difference of the two single-replicate posterior means and is exact for
the mean difference under an independent-posterior interpretation.

However:

    dx_rep_latent_median
    dx_rep_latent_mode

are differences of marginal summaries, not the median or mode of a formally
propagated difference posterior.

Likewise:

    dx_rep_latent_q025
        = rep2_x_q025 - rep1_x_q975

    dx_rep_latent_q975
        = rep2_x_q975 - rep1_x_q025

are conservative endpoint-derived bounds, not calibrated posterior quantiles
of the replicate difference.

The approximate standard deviation is:

    dx_rep_latent_sd_approx
        =
        sqrt(
            rep1_x_latent_sd^2
            +
            rep2_x_latent_sd^2
        ).

These are diagnostic or sensitivity quantities. They are not the observed
replicate-decoupled forward-drift estimator used as the final primary
forward estimand.

Posterior-predictive replicate detectability
--------------------------------------------
For latent frequency f and replicate depth N, detection is defined by:

    count >= detection_count_threshold.

Default:

    detection_count_threshold = 1.

For threshold 1:

    P(detect | f,N,kappa)
        = 1 - P(C=0 | f,N,kappa).

The reported:

    p_detect_rep1
    p_detect_rep2

are obtained by integrating the corresponding replicate-detection probability
over the JOINT paired-replicate posterior:

    P(f | c_rep1,c_rep2).

They are not calculated from the single-replicate posterior summaries.

Dropout fields are:

    p_dropout_rep1 = 1 - p_detect_rep1
    p_dropout_rep2 = 1 - p_detect_rep2.

Paired-state posterior-predictive detectability
-----------------------------------------------
The latent timepoint state is defined as technically detected when at least one
of the two technical replicates reaches the configured count threshold.

Conditional on f, replicate observations are independent under the fitted NB
model.

Let:

    q_1(f)
        = P(C_1 < threshold | f)

    q_2(f)
        = P(C_2 < threshold | f).

The implemented paired-state dropout probability is:

    p_dropout_state
        =
        E_post[
            q_1(f) q_2(f)
        ],

where the expectation is under the SHARED joint posterior:

    P(f | c_1,c_2).

Then:

    p_detect_state
        = 1 - p_dropout_state.

This is generally NOT identical to:

    1 - (1-p_detect_rep1)(1-p_detect_rep2),

because the replicate predictive probabilities are integrated over the same
uncertain shared latent frequency. Their difference is therefore a consequence
of posterior coupling and should not be interpreted as failure of an algebraic
identity.

These quantities are model-based posterior-predictive detectability measures.
They are not independently calibrated detection probabilities and are distinct
from the p_value-based operational observability annotation created by Step 2.

Full posterior storage
----------------------
When:

    posterior_npz_path

is supplied to:

    fit_noiseK_latent_powerlaw(...),

the module writes one compressed NPZ file containing exactly:

    f_grid
        float64 latent-frequency grid;

    posterior_matrix
        float32 matrix with shape:

            n_clonotypes x grid_size;

    aaSeqCDR3
        clonotype identifiers as strings.

No x-grid is stored because downstream code can reconstruct:

    x_grid = log(f_grid).

The cache supports latent-state diagnostics, sensitivity analyses and optional
posterior propagation. It does not define the primary observed displacement
estimators.

Endpoint-posterior displacement propagation
-------------------------------------------
The public function:

    posterior_logfold_from_latent_posteriors(...)

propagates two latent-frequency posteriors into a log-frequency displacement:

    dx = log(f_1) - log(f_0).

It samples independently from the two supplied marginal posterior grids:

    f_0^(m) ~ P_0(f)
    f_1^(m) ~ P_1(f)

and computes:

    dx^(m)
        = log f_1^(m) - log f_0^(m).

The returned summaries are:

    dx_latent_mean
    dx_latent_median
    dx_latent_q025
    dx_latent_q975
    p_dx_gt0
    dx_latent_sd
    n_samples.

The default requested Monte-Carlo sample count is:

    2000,

with a hard minimum of:

    100.

The function assumes INDEPENDENT endpoint posterior draws. It does not
construct a longitudinal joint posterior linking the two timepoints.

In the final pipeline these propagated latent displacements are optional latent
quantities and should not be conflated with the primary replicate-decoupled
observed forward estimator or primary observed cross-replicate covariance.

Single-count displacement compatibility function
------------------------------------------------
The convenience function:

    posterior_logfold_from_counts(...)

first constructs two single-count posteriors and then calls:

    posterior_logfold_from_latent_posteriors(...).

Separate:

    gamma0, kappa0
    gamma1, kappa1

can be supplied through the existing compatibility signature.

Public fit result
-----------------
The module returns:

    FitResult(
        success,
        message,
        log_likelihood,
        params,
        per_clone_logP
    ).

success and message come directly from the numerical optimizer.

The module constructs the fitted parameter vector and posterior outputs from
the optimizer result even if the optimizer success flag is False. Callers must
therefore inspect success explicitly. The production Step-2 workflow is
responsible for deciding whether an unsuccessful fit is eligible for canonical
downstream aggregation.

Returned fitted parameters
--------------------------
The params dictionary contains:

    gamma
    k
    fmin

    N_total_mean
    N_total_min
    N_total_max

    grid_size
    chunk_size

    detection_count_threshold.

Here:

    k

is the fitted Negative-Binomial size parameter kappa.

Backward compatibility
----------------------
The alias:

    fit_noiseK_nb_powerlaw

points to:

    fit_noiseK_latent_powerlaw.

This preserves compatibility with older callers expecting the earlier function
name while returning the expanded latent-posterior output structure.

Relationship to Step 2
----------------------
The intended production caller is:

    2-multirepresentation_clonotype_state_inference.py.

Step 2 is responsible for:

    parsing raw repertoire files;
    collapsing duplicate aaSeqCDR3 rows by summing readCount;
    pairing technical replicates;
    outer-merging clonotypes;
    assigning zero counts to replicate absences;
    supplying replicate depths;
    attaching subject/time metadata;
    constructing replicate-resolved observed-frequency fields;
    constructing auxiliary conditional-positive summaries;
    writing pair-level and subject-level output tables;
    creating empirical p_value-based operational observability annotations;
    deciding which successful fit outputs enter the canonical production
    aggregation.

This module is responsible only for the statistical latent-frequency inference
and posterior-predictive quantities described above.

Interpretive boundary
---------------------
The module reconstructs latent-frequency uncertainty under a specified
technical-replicate observation model.

It should not be described as:

    a longitudinal dynamical model;

    an estimator of biological birth, death, extinction or persistence;

    proof that observed count differences are purely technical;

    an independently calibrated detectability model;

    the primary forward-drift estimator;

    the primary fluctuation estimator;

    proof that the joint latent posterior is the true clonotype abundance;

    exact implementation of another software package unless separately
    demonstrated.

The principal inferential object of this module is the posterior distribution
over a model-based shared latent clonotype frequency for one biological
timepoint observed through paired technical replicates.

Dependencies
------------
Required:

    numpy
    pandas
    scipy.

Minimal smoke test
------------------
When executed directly as:

    python3 noiseK_latent.py

the module runs a small synthetic paired-replicate example and prints:

    optimizer success/message;
    fitted parameter dictionary;
    first rows of the per-clonotype output.
"""


from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp


@dataclass
class FitResult:
    success: bool
    message: str
    log_likelihood: float
    params: Dict[str, float]
    per_clone_logP: pd.DataFrame


# =============================================================================
# Numerical utilities
# =============================================================================


def _as_positive_float(x: float, fallback: float) -> float:
    try:
        x = float(x)
    except Exception:
        return float(fallback)
    if not np.isfinite(x) or x <= 0:
        return float(fallback)
    return x


def _make_frequency_grid(
    depths: List[int],
    grid_size: int = 500,
    fmin: Optional[float] = None,
    fmax: float = 1.0,
) -> Tuple[np.ndarray, float]:
    """Create the log-spaced latent-frequency grid."""
    if len(depths) == 0:
        raise ValueError("depths must contain at least one positive depth")

    depths = [int(d) for d in depths]
    if any(d <= 0 for d in depths):
        raise ValueError(f"All depths must be positive; got {depths}")

    if fmin is None:
        # Paired-replicate default: approximately one read at the mean replicate depth.
        n_total = float(np.mean(depths))
        fmin = 1.0 / max(n_total, 1.0)

    fmin = float(fmin)
    fmax = float(fmax)
    grid_size = int(max(8, grid_size))

    if not np.isfinite(fmin) or fmin <= 0:
        raise ValueError(f"fmin must be positive; got {fmin}")
    if not np.isfinite(fmax) or fmax <= fmin:
        raise ValueError(f"fmax must be finite and > fmin; got fmax={fmax}, fmin={fmin}")

    return np.logspace(np.log10(fmin), np.log10(fmax), grid_size).astype(np.float64), fmin


def _log_nb_pmf_matrix(
    counts: np.ndarray,
    mu_grid: np.ndarray,
    kappa: float,
) -> np.ndarray:
    """
    Negative-binomial log PMF for many counts against one mean grid.

    Parameterization:
        mean = mu
        size = kappa

    PMF:
        Gamma(c+r)/(Gamma(r)c!) * (r/(r+mu))^r * (mu/(r+mu))^c
    """
    c = np.asarray(counts, dtype=np.int64)
    mu = np.asarray(mu_grid, dtype=np.float64)
    r = float(kappa)

    if r <= 0 or not np.isfinite(r):
        raise ValueError(f"kappa must be positive; got {kappa}")
    if np.any(mu <= 0) or np.any(~np.isfinite(mu)):
        raise ValueError("mu_grid must be finite and positive")

    log_mu = np.log(mu)
    log_r = np.log(r)
    log_r_plus_mu = np.log(r + mu)

    return (
        gammaln(c[:, None] + r)
        - gammaln(r)
        - gammaln(c[:, None] + 1)
        + r * (log_r - log_r_plus_mu[None, :])
        + c[:, None] * (log_mu[None, :] - log_r_plus_mu[None, :])
    )


def _make_log_prior(f_grid: np.ndarray, gamma: float) -> np.ndarray:
    gamma = float(gamma)
    if gamma <= 0 or not np.isfinite(gamma):
        raise ValueError(f"gamma must be positive; got {gamma}")
    log_prior = -gamma * np.log(f_grid)
    log_prior -= logsumexp(log_prior)
    return log_prior


def _posterior_summaries_from_logpost(
    f_grid: np.ndarray,
    log_post: np.ndarray,
    q_low: float = 0.025,
    q_high: float = 0.975,
) -> Dict[str, float]:
    """Return frequency and log-frequency posterior summaries."""
    log_post = np.asarray(log_post, dtype=np.float64)
    log_post = log_post - logsumexp(log_post)
    p = np.exp(log_post)
    cdf = np.cumsum(p)
    cdf[-1] = 1.0

    mean_f = float(np.sum(p * f_grid))
    median_f = float(np.interp(0.5, cdf, f_grid))
    q025_f = float(np.interp(float(q_low), cdf, f_grid))
    q975_f = float(np.interp(float(q_high), cdf, f_grid))

    x_grid = np.log(f_grid)
    mean_x = float(np.sum(p * x_grid))
    median_x = float(np.interp(0.5, cdf, x_grid))
    q025_x = float(np.interp(float(q_low), cdf, x_grid))
    q975_x = float(np.interp(float(q_high), cdf, x_grid))

    # Posterior mode on the discrete grid.
    imax = int(np.argmax(p))
    mode_f = float(f_grid[imax])
    mode_x = float(x_grid[imax])

    return {
        "f_latent_mean": mean_f,
        "f_latent_median": median_f,
        "f_latent_mode": mode_f,
        "f_latent_q025": q025_f,
        "f_latent_q975": q975_f,
        "x_latent_mean": mean_x,
        "x_latent_median": median_x,
        "x_latent_mode": mode_x,
        "x_latent_q025": q025_x,
        "x_latent_q975": q975_x,
        "x_latent_sd": float(np.sqrt(np.sum(p * (x_grid - mean_x) ** 2))),
        "posterior_entropy": float(-np.sum(p * np.log(np.maximum(p, 1e-300)))),
    }


def _posterior_summaries_from_post(
    f_grid: np.ndarray,
    post: np.ndarray,
    q_low: float = 0.025,
    q_high: float = 0.975,
    prefix: str = "",
) -> Dict[str, float]:
    """
    Return posterior summaries from an already-normalized posterior probability
    vector. If prefix is provided, output keys are prefixed, e.g.
    prefix="rep1_" -> rep1_f_latent_median.
    """
    post = np.asarray(post, dtype=np.float64)
    s = float(np.sum(post))
    if not np.isfinite(s) or s <= 0:
        out = _posterior_summaries_from_logpost(f_grid, np.full_like(f_grid, -np.log(len(f_grid))))
    else:
        log_post = np.log(np.maximum(post / s, 1e-300))
        out = _posterior_summaries_from_logpost(f_grid, log_post, q_low=q_low, q_high=q_high)

    if prefix:
        return {f"{prefix}{k}": v for k, v in out.items()}
    return out


def _detection_prob_from_posterior(
    f_grid: np.ndarray,
    post: np.ndarray,
    depth: int,
    kappa: float,
    count_threshold: int = 1,
) -> float:
    """
    Posterior predictive probability that a clone would be detected in a replicate.

    Detection is defined as count >= count_threshold. Default threshold is 1 read.
    For threshold=1, P(detect | f) = 1 - P(c=0 | f).
    """
    depth = int(depth)
    count_threshold = int(count_threshold)
    if depth <= 0:
        return float("nan")
    if count_threshold <= 0:
        return 1.0

    mu = f_grid * float(depth)
    r = float(kappa)

    # Sum NB probabilities for c=0,...,threshold-1.
    probs_below = np.zeros_like(f_grid, dtype=np.float64)
    for c in range(count_threshold):
        logp = (
            gammaln(c + r)
            - gammaln(r)
            - gammaln(c + 1)
            + r * (np.log(r) - np.log(r + mu))
            + c * (np.log(mu) - np.log(r + mu))
        )
        probs_below += np.exp(logp)

    p_detect_grid = 1.0 - probs_below
    return float(np.sum(post * p_detect_grid))


def _paired_state_detection_prob_from_posterior(
    f_grid: np.ndarray,
    post: np.ndarray,
    depth1: int,
    depth2: int,
    kappa: float,
    count_threshold: int = 1,
) -> Tuple[float, float]:
    """
    Posterior-predictive detectability of one latent time-point state under the
    paired technical-replicate design.

    The latent state is considered detected when at least one of the two
    replicates yields a count >= count_threshold. Conditional on latent
    frequency f, the two technical replicates are independent under the fitted
    Negative-Binomial observation model. The joint dropout probability is
    therefore integrated over the shared posterior P(f | c1, c2):

        p_dropout_state = E_post[P(C1 < threshold | f)
                                   * P(C2 < threshold | f)]
        p_detect_state  = 1 - p_dropout_state

    Returns
    -------
    (p_detect_state, p_dropout_state)
    """
    depth1 = int(depth1)
    depth2 = int(depth2)
    count_threshold = int(count_threshold)

    if depth1 <= 0 or depth2 <= 0:
        return float("nan"), float("nan")
    if count_threshold <= 0:
        return 1.0, 0.0

    f_grid = np.asarray(f_grid, dtype=np.float64)
    post = np.asarray(post, dtype=np.float64)
    post_sum = float(np.sum(post))
    if not np.isfinite(post_sum) or post_sum <= 0:
        return float("nan"), float("nan")
    post = post / post_sum

    r = float(kappa)
    if not np.isfinite(r) or r <= 0:
        return float("nan"), float("nan")

    mu1 = f_grid * float(depth1)
    mu2 = f_grid * float(depth2)

    p_below1 = np.zeros_like(f_grid, dtype=np.float64)
    p_below2 = np.zeros_like(f_grid, dtype=np.float64)

    for c in range(count_threshold):
        logp1 = (
            gammaln(c + r)
            - gammaln(r)
            - gammaln(c + 1)
            + r * (np.log(r) - np.log(r + mu1))
            + c * (np.log(mu1) - np.log(r + mu1))
        )
        logp2 = (
            gammaln(c + r)
            - gammaln(r)
            - gammaln(c + 1)
            + r * (np.log(r) - np.log(r + mu2))
            + c * (np.log(mu2) - np.log(r + mu2))
        )
        p_below1 += np.exp(logp1)
        p_below2 += np.exp(logp2)

    p_dropout_state = float(np.sum(post * p_below1 * p_below2))
    p_dropout_state = float(np.clip(p_dropout_state, 0.0, 1.0))
    p_detect_state = float(1.0 - p_dropout_state)
    return p_detect_state, p_dropout_state


# =============================================================================
# Public posterior functions
# =============================================================================


def posterior_latent_frequency_from_replicates(
    count1: int,
    count2: int,
    depth1: int,
    depth2: int,
    gamma: float,
    kappa: float,
    fmin: Optional[float] = None,
    grid_size: int = 500,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute posterior P(f | count1, count2) on a frequency grid.

    Absent counts are allowed. For example, count1=0 and count2>0 is treated as
    a valid Negative-Binomial zero-count observation rather than as true zero frequency.
    """
    count1 = int(count1)
    count2 = int(count2)
    if count1 < 0 or count2 < 0:
        raise ValueError("counts must be non-negative")

    f_grid, _fmin = _make_frequency_grid([int(depth1), int(depth2)], grid_size=grid_size, fmin=fmin)
    log_prior = _make_log_prior(f_grid, gamma)

    mu1 = f_grid * float(depth1)
    mu2 = f_grid * float(depth2)

    log_like1 = _log_nb_pmf_matrix(np.asarray([count1]), mu1, kappa)[0]
    log_like2 = _log_nb_pmf_matrix(np.asarray([count2]), mu2, kappa)[0]

    log_post = log_like1 + log_like2 + log_prior
    log_post -= logsumexp(log_post)
    return f_grid, np.exp(log_post)


def posterior_frequency_from_count(
    count: int,
    depth: int,
    gamma: float,
    kappa: float,
    fmin: Optional[float] = None,
    grid_size: int = 500,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Single-replicate posterior P(f | count), retained for compatibility with
    earlier downstream usage.
    """
    count = int(count)
    if count < 0:
        raise ValueError("count must be non-negative")

    f_grid, _fmin = _make_frequency_grid([int(depth)], grid_size=grid_size, fmin=fmin)
    log_prior = _make_log_prior(f_grid, gamma)
    mu = f_grid * float(depth)
    log_like = _log_nb_pmf_matrix(np.asarray([count]), mu, kappa)[0]
    log_post = log_like + log_prior
    log_post -= logsumexp(log_post)
    return f_grid, np.exp(log_post)


def posterior_logfold_from_latent_posteriors(
    f_grid0: np.ndarray,
    post0: np.ndarray,
    f_grid1: np.ndarray,
    post1: np.ndarray,
    n_samples: int = 2000,
    random_state: Optional[int] = None,
) -> Dict[str, float]:
    """
    Monte-Carlo posterior summary of dx = log(f1) - log(f0).

    This function will be used downstream by latent transition reconstruction.
    It intentionally samples from posterior distributions instead of subtracting
    point estimates, thereby propagating endpoint posterior uncertainty.
    """
    rng = np.random.default_rng(random_state)
    n_samples = int(max(100, n_samples))

    f_grid0 = np.asarray(f_grid0, dtype=np.float64)
    f_grid1 = np.asarray(f_grid1, dtype=np.float64)
    post0 = np.asarray(post0, dtype=np.float64)
    post1 = np.asarray(post1, dtype=np.float64)

    post0 = post0 / np.sum(post0)
    post1 = post1 / np.sum(post1)

    x0 = np.log(rng.choice(f_grid0, size=n_samples, replace=True, p=post0))
    x1 = np.log(rng.choice(f_grid1, size=n_samples, replace=True, p=post1))
    dx = x1 - x0

    return {
        "dx_latent_mean": float(np.mean(dx)),
        "dx_latent_median": float(np.median(dx)),
        "dx_latent_q025": float(np.quantile(dx, 0.025)),
        "dx_latent_q975": float(np.quantile(dx, 0.975)),
        "p_dx_gt0": float(np.mean(dx > 0)),
        "dx_latent_sd": float(np.std(dx, ddof=1)),
        "n_samples": int(n_samples),
    }


# Compatibility alias expected by some older code names.
def posterior_logfold_from_counts(
    count0: int,
    count1: int,
    depth0: int,
    depth1: int,
    gamma0: float,
    k0: float,
    gamma1: Optional[float] = None,
    k1: Optional[float] = None,
    fmin0: Optional[float] = None,
    fmin1: Optional[float] = None,
    grid_size: int = 500,
    n_samples: int = 2000,
    random_state: Optional[int] = None,
) -> Dict[str, float]:
    """Single-count version of posterior log-fold change."""
    if gamma1 is None:
        gamma1 = gamma0
    if k1 is None:
        k1 = k0

    fg0, p0 = posterior_frequency_from_count(count0, depth0, gamma0, k0, fmin=fmin0, grid_size=grid_size)
    fg1, p1 = posterior_frequency_from_count(count1, depth1, gamma1, k1, fmin=fmin1, grid_size=grid_size)
    return posterior_logfold_from_latent_posteriors(
        fg0,
        p0,
        fg1,
        p1,
        n_samples=n_samples,
        random_state=random_state,
    )


# =============================================================================
# Fitting with latent posterior output
# =============================================================================


def fit_noiseK_latent_powerlaw(
    merged: pd.DataFrame,
    depths: List[int],
    count_cols: List[str],
    gamma_init: float = 1.6,
    k_init: float = 50.0,
    grid_size: int = 500,
    chunk_size: int = 20_000,
    maxiter: Optional[int] = None,
    optimizer_ftol: Optional[float] = None,
    detection_count_threshold: int = 1,
    posterior_npz_path: Optional[str] = None,
) -> FitResult:
    """
    Fit the NB + power-law prior model and return latent-frequency summaries.

    Public arguments retain compatibility with the earlier noise-model fitting API.
    This function additionally returns latent-frequency posterior summaries.

    Output columns added to per_clone_logP
    -------------------------------------
    logP

    Joint posterior from the paired technical replicates:
        f_latent_mean, f_latent_median, f_latent_mode, f_latent_q025, f_latent_q975
        x_latent_mean, x_latent_median, x_latent_mode, x_latent_q025, x_latent_q975
        x_latent_sd, posterior_entropy

    Single-replicate posterior summaries:
        rep1_f_latent_*, rep1_x_latent_*
        rep2_f_latent_*, rep2_x_latent_*

    Convenience aliases for downstream diagnostics:
        f_rep1_latent_median, f_rep2_latent_median
        x_rep1_latent_median, x_rep2_latent_median
        x_rep1_latent_q025/q975, x_rep2_latent_q025/q975
        dx_rep_latent_median, dx_rep_latent_mean, dx_rep_latent_mode
        dx_rep_latent_q025, dx_rep_latent_q975, dx_rep_latent_sd_approx
        replicate_posterior_bhattacharyya

    Detection/dropout posterior predictive summaries:
        p_detect_state, p_dropout_state
        p_detect_rep1, p_detect_rep2, p_dropout_rep1, p_dropout_rep2

    p_detect_state is the posterior-predictive probability that the clonotype
    would be detected in at least one of the two technical replicates.
    """
    if len(depths) != 2:
        raise ValueError("depths must contain exactly two values: [depth1, depth2]")
    if len(count_cols) != 2:
        raise ValueError("count_cols must contain exactly two column names")

    c1 = pd.to_numeric(merged[count_cols[0]], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    c2 = pd.to_numeric(merged[count_cols[1]], errors="coerce").fillna(0).to_numpy(dtype=np.int64)

    if c1.size == 0:
        raise ValueError("merged contains zero clonotypes")
    if np.any(c1 < 0) or np.any(c2 < 0):
        raise ValueError("counts must be non-negative")

    depth1, depth2 = int(depths[0]), int(depths[1])
    if depth1 <= 0 or depth2 <= 0:
        raise ValueError(f"depths must be positive; got depth1={depth1}, depth2={depth2}")

    gamma_init = _as_positive_float(gamma_init, 1.6)
    k_init = _as_positive_float(k_init, 50.0)
    grid_size = int(max(8, grid_size))
    chunk_size = int(max(1, chunk_size))

    f_grid, fmin = _make_frequency_grid([depth1, depth2], grid_size=grid_size)
    log_f_grid = np.log(f_grid)
    mu1 = f_grid * float(depth1)
    mu2 = f_grid * float(depth2)

    def neg_log_likelihood(theta: np.ndarray) -> float:
        gamma = float(theta[0])
        kappa = float(theta[1])
        if not np.isfinite(gamma) or not np.isfinite(kappa) or gamma <= 0.0 or kappa <= 0.0:
            return np.inf

        log_prior = _make_log_prior(f_grid, gamma)
        total = 0.0
        for start in range(0, c1.size, chunk_size):
            stop = min(start + chunk_size, c1.size)
            log_like1 = _log_nb_pmf_matrix(c1[start:stop], mu1, kappa)
            log_like2 = _log_nb_pmf_matrix(c2[start:stop], mu2, kappa)
            total += float(np.sum(logsumexp(log_like1 + log_like2 + log_prior[None, :], axis=1)))
        return -total

    options = {}
    if maxiter is not None:
        options["maxiter"] = int(maxiter)
    if optimizer_ftol is not None:
        options["ftol"] = float(optimizer_ftol)

    res = minimize(
        neg_log_likelihood,
        x0=np.array([gamma_init, k_init], dtype=np.float64),
        method="L-BFGS-B",
        bounds=[(1e-6, None), (1e-6, None)],
        options=options if options else None,
    )

    gamma_hat, k_hat = map(float, res.x)
    ll = -float(res.fun)
    log_prior_hat = _make_log_prior(f_grid, gamma_hat)

    out = merged.copy()

    # Preallocate outputs.
    n = c1.size
    posterior_matrix = None
    if posterior_npz_path is not None:
        posterior_matrix = np.empty((n, grid_size), dtype=np.float32)
    logP = np.empty(n, dtype=np.float64)
    summary_cols = {
        # Joint posterior from both technical replicates.
        "f_latent_mean": np.empty(n, dtype=np.float64),
        "f_latent_median": np.empty(n, dtype=np.float64),
        "f_latent_mode": np.empty(n, dtype=np.float64),
        "f_latent_q025": np.empty(n, dtype=np.float64),
        "f_latent_q975": np.empty(n, dtype=np.float64),
        "x_latent_mean": np.empty(n, dtype=np.float64),
        "x_latent_median": np.empty(n, dtype=np.float64),
        "x_latent_mode": np.empty(n, dtype=np.float64),
        "x_latent_q025": np.empty(n, dtype=np.float64),
        "x_latent_q975": np.empty(n, dtype=np.float64),
        "x_latent_sd": np.empty(n, dtype=np.float64),
        "posterior_entropy": np.empty(n, dtype=np.float64),

        # Single-replicate posterior summaries.
        "rep1_f_latent_mean": np.empty(n, dtype=np.float64),
        "rep1_f_latent_median": np.empty(n, dtype=np.float64),
        "rep1_f_latent_mode": np.empty(n, dtype=np.float64),
        "rep1_f_latent_q025": np.empty(n, dtype=np.float64),
        "rep1_f_latent_q975": np.empty(n, dtype=np.float64),
        "rep1_x_latent_mean": np.empty(n, dtype=np.float64),
        "rep1_x_latent_median": np.empty(n, dtype=np.float64),
        "rep1_x_latent_mode": np.empty(n, dtype=np.float64),
        "rep1_x_latent_q025": np.empty(n, dtype=np.float64),
        "rep1_x_latent_q975": np.empty(n, dtype=np.float64),
        "rep1_x_latent_sd": np.empty(n, dtype=np.float64),
        "rep1_posterior_entropy": np.empty(n, dtype=np.float64),

        "rep2_f_latent_mean": np.empty(n, dtype=np.float64),
        "rep2_f_latent_median": np.empty(n, dtype=np.float64),
        "rep2_f_latent_mode": np.empty(n, dtype=np.float64),
        "rep2_f_latent_q025": np.empty(n, dtype=np.float64),
        "rep2_f_latent_q975": np.empty(n, dtype=np.float64),
        "rep2_x_latent_mean": np.empty(n, dtype=np.float64),
        "rep2_x_latent_median": np.empty(n, dtype=np.float64),
        "rep2_x_latent_mode": np.empty(n, dtype=np.float64),
        "rep2_x_latent_q025": np.empty(n, dtype=np.float64),
        "rep2_x_latent_q975": np.empty(n, dtype=np.float64),
        "rep2_x_latent_sd": np.empty(n, dtype=np.float64),
        "rep2_posterior_entropy": np.empty(n, dtype=np.float64),

        # Convenience aliases for replicate-specific diagnostics.
        "f_rep1_latent_median": np.empty(n, dtype=np.float64),
        "f_rep2_latent_median": np.empty(n, dtype=np.float64),
        "x_rep1_latent_median": np.empty(n, dtype=np.float64),
        "x_rep2_latent_median": np.empty(n, dtype=np.float64),
        "x_rep1_latent_q025": np.empty(n, dtype=np.float64),
        "x_rep1_latent_q975": np.empty(n, dtype=np.float64),
        "x_rep2_latent_q025": np.empty(n, dtype=np.float64),
        "x_rep2_latent_q975": np.empty(n, dtype=np.float64),
        "dx_rep_latent_mean": np.empty(n, dtype=np.float64),
        "dx_rep_latent_median": np.empty(n, dtype=np.float64),
        "dx_rep_latent_mode": np.empty(n, dtype=np.float64),
        "dx_rep_latent_q025": np.empty(n, dtype=np.float64),
        "dx_rep_latent_q975": np.empty(n, dtype=np.float64),
        "dx_rep_latent_sd_approx": np.empty(n, dtype=np.float64),
        "replicate_posterior_bhattacharyya": np.empty(n, dtype=np.float64),

        # Detection/dropout from the joint posterior.
        # State-level values summarize the paired-replicate observation design;
        # replicate-specific values are retained for technical QC.
        "p_detect_state": np.empty(n, dtype=np.float64),
        "p_dropout_state": np.empty(n, dtype=np.float64),
        "p_detect_rep1": np.empty(n, dtype=np.float64),
        "p_detect_rep2": np.empty(n, dtype=np.float64),
        "p_dropout_rep1": np.empty(n, dtype=np.float64),
        "p_dropout_rep2": np.empty(n, dtype=np.float64),
    }

    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        log_like1 = _log_nb_pmf_matrix(c1[start:stop], mu1, k_hat)
        log_like2 = _log_nb_pmf_matrix(c2[start:stop], mu2, k_hat)
        log_joint = log_like1 + log_like2 + log_prior_hat[None, :]
        row_logP = logsumexp(log_joint, axis=1)
        logP[start:stop] = row_logP
        log_post = log_joint - row_logP[:, None]
        post = np.exp(log_post)
        if posterior_matrix is not None:
            posterior_matrix[start:stop, :] = post.astype(np.float32)

        # Vectorized posterior summaries where easy.
        cdf = np.cumsum(post, axis=1)
        cdf[:, -1] = 1.0
        x_grid = log_f_grid

        for j in range(stop - start):
            idx = start + j

            # -------------------------------------------------------------
            # Joint latent posterior from the paired technical replicates:
            #     p(f | c_rep1, c_rep2)
            # This is the consensus latent frequency used for trajectories.
            # -------------------------------------------------------------
            summaries = _posterior_summaries_from_logpost(f_grid, log_post[j])
            for key in [
                "f_latent_mean", "f_latent_median", "f_latent_mode", "f_latent_q025", "f_latent_q975",
                "x_latent_mean", "x_latent_median", "x_latent_mode", "x_latent_q025", "x_latent_q975",
                "x_latent_sd", "posterior_entropy",
            ]:
                summary_cols[key][idx] = summaries[key]

            # -------------------------------------------------------------
            # Single-replicate latent posteriors:
            #     p(f | c_rep1) and p(f | c_rep2)
            # These are NOT used as the consensus timepoint estimate. They
            # are saved to diagnose whether the latent model reduces the
            # replicate-level posterior disagreement and boundary effects.
            # -------------------------------------------------------------
            log_post_rep1 = log_like1[j] + log_prior_hat
            log_post_rep1 -= logsumexp(log_post_rep1)
            post_rep1 = np.exp(log_post_rep1)

            log_post_rep2 = log_like2[j] + log_prior_hat
            log_post_rep2 -= logsumexp(log_post_rep2)
            post_rep2 = np.exp(log_post_rep2)

            rep1_summaries = _posterior_summaries_from_post(f_grid, post_rep1, prefix="rep1_")
            rep2_summaries = _posterior_summaries_from_post(f_grid, post_rep2, prefix="rep2_")

            # Exact Bhattacharyya coefficient between the two single-replicate
            # posterior distributions, evaluated on the shared discrete f-grid.
            summary_cols["replicate_posterior_bhattacharyya"][idx] = float(
                np.clip(np.sum(np.sqrt(post_rep1 * post_rep2)), 0.0, 1.0)
            )

            for key, value in rep1_summaries.items():
                summary_cols[key][idx] = value
            for key, value in rep2_summaries.items():
                summary_cols[key][idx] = value

            # Stable downstream aliases.
            summary_cols["f_rep1_latent_median"][idx] = rep1_summaries["rep1_f_latent_median"]
            summary_cols["f_rep2_latent_median"][idx] = rep2_summaries["rep2_f_latent_median"]
            summary_cols["x_rep1_latent_median"][idx] = rep1_summaries["rep1_x_latent_median"]
            summary_cols["x_rep2_latent_median"][idx] = rep2_summaries["rep2_x_latent_median"]
            summary_cols["x_rep1_latent_q025"][idx] = rep1_summaries["rep1_x_latent_q025"]
            summary_cols["x_rep1_latent_q975"][idx] = rep1_summaries["rep1_x_latent_q975"]
            summary_cols["x_rep2_latent_q025"][idx] = rep2_summaries["rep2_x_latent_q025"]
            summary_cols["x_rep2_latent_q975"][idx] = rep2_summaries["rep2_x_latent_q975"]

            # Approximate replicate latent displacement summaries.
            # For independent single-replicate posteriors, the mean difference
            # is exact. Median/quantile approximations are deterministic and
            # sufficient for diagnostics; downstream formal transition posteriors
            # should sample full posterior grids when needed.
            summary_cols["dx_rep_latent_mean"][idx] = (
                rep2_summaries["rep2_x_latent_mean"] - rep1_summaries["rep1_x_latent_mean"]
            )
            summary_cols["dx_rep_latent_median"][idx] = (
                rep2_summaries["rep2_x_latent_median"] - rep1_summaries["rep1_x_latent_median"]
            )
            summary_cols["dx_rep_latent_mode"][idx] = (
                rep2_summaries["rep2_x_latent_mode"] - rep1_summaries["rep1_x_latent_mode"]
            )
            summary_cols["dx_rep_latent_q025"][idx] = (
                rep2_summaries["rep2_x_latent_q025"] - rep1_summaries["rep1_x_latent_q975"]
            )
            summary_cols["dx_rep_latent_q975"][idx] = (
                rep2_summaries["rep2_x_latent_q975"] - rep1_summaries["rep1_x_latent_q025"]
            )
            summary_cols["dx_rep_latent_sd_approx"][idx] = float(
                np.sqrt(rep1_summaries["rep1_x_latent_sd"] ** 2 + rep2_summaries["rep2_x_latent_sd"] ** 2)
            )

            pdet1 = _detection_prob_from_posterior(
                f_grid,
                post[j],
                depth1,
                k_hat,
                count_threshold=detection_count_threshold,
            )
            pdet2 = _detection_prob_from_posterior(
                f_grid,
                post[j],
                depth2,
                k_hat,
                count_threshold=detection_count_threshold,
            )
            pdet_state, pdrop_state = _paired_state_detection_prob_from_posterior(
                f_grid,
                post[j],
                depth1,
                depth2,
                k_hat,
                count_threshold=detection_count_threshold,
            )

            summary_cols["p_detect_state"][idx] = pdet_state
            summary_cols["p_dropout_state"][idx] = pdrop_state
            summary_cols["p_detect_rep1"][idx] = pdet1
            summary_cols["p_detect_rep2"][idx] = pdet2
            summary_cols["p_dropout_rep1"][idx] = 1.0 - pdet1
            summary_cols["p_dropout_rep2"][idx] = 1.0 - pdet2

    out["logP"] = logP
    for key, values in summary_cols.items():
        out[key] = values
    
    if posterior_matrix is not None:
        np.savez_compressed(
            posterior_npz_path,
            f_grid=f_grid.astype(np.float64),
            posterior_matrix=posterior_matrix,
            aaSeqCDR3=out["aaSeqCDR3"].astype(str).to_numpy(),
    )

    params = {
        "gamma": float(gamma_hat),
        "k": float(k_hat),
        "fmin": float(fmin),
        "N_total_mean": float(0.5 * (depth1 + depth2)),
        "N_total_min": float(min(depth1, depth2)),
        "N_total_max": float(max(depth1, depth2)),
        "grid_size": int(grid_size),
        "chunk_size": int(chunk_size),
        "detection_count_threshold": int(detection_count_threshold),
    }

    return FitResult(
        success=bool(res.success),
        message=str(res.message),
        log_likelihood=float(ll),
        params=params,
        per_clone_logP=out,
    )


# Backward-compatible alias: if a pipeline imports fit_noiseK_nb_powerlaw from
# this module, it will receive latent summaries while preserving the FitResult API.
fit_noiseK_nb_powerlaw = fit_noiseK_latent_powerlaw


if __name__ == "__main__":
    # Minimal smoke test on a tiny synthetic replicate pair.
    demo = pd.DataFrame(
        {
            "aaSeqCDR3": ["CASSA", "CASSB", "CASSC", "CASSD"],
            "count_rep1": [100, 10, 1, 0],
            "count_rep2": [90, 0, 2, 5],
        }
    )
    result = fit_noiseK_latent_powerlaw(
        demo,
        depths=[1_000_000, 950_000],
        count_cols=["count_rep1", "count_rep2"],
        grid_size=200,
        chunk_size=1000,
        maxiter=20,
    )
    print(result.success, result.message)
    print(result.params)
    print(result.per_clone_logP.head())