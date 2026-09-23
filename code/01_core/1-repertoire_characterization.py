#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1-repertoire_characterization.py
================================

Observed repertoire characterization and descriptive heavy-tail analysis for
ClonoDynamics Step 1.

Overview
--------
This script implements Step 1 of the ClonoDynamics pipeline. Each technical
replicate is analysed as an individual observed TCRbeta repertoire. The step
summarizes sequencing depth and observed clonotype number and characterizes the
high-frequency tail of the observed clonotype-frequency distribution.

Step 1 is deliberately DESCRIPTIVE. It does not infer clonotype latent states,
does not estimate the pair-specific Negative-Binomial parameters used in Step 2,
does not define operational observability, and does not filter repertoires or
clonotypes. In particular, the power-law exponent estimated here is a
descriptor of the observed repertoire tail; it is NOT passed to Step 2 and does
NOT numerically parameterize or validate the discrete power-law-shaped prior
used by the Step-2 latent-state model.

Position in the pipeline
------------------------

    Step 1
    1-repertoire_characterization.py
        -> observed repertoire architecture
        -> descriptive heavy-tail characterization
        -> Pareto-versus-truncated-log-normal tail comparison
        -> no filtering
                    |
                    v

    Step 2
    2-multirepresentation_clonotype_state_inference.py
        -> independent pair-specific gamma and kappa estimation
        -> latent consensus representation
        -> replicate-resolved observed measurements
        -> detection / operational-observation annotations

The two steps therefore address different statistical questions. Step 1 asks
whether the observed repertoire architecture is heavy-tailed and how well a
parsimonious Pareto description captures its upper tail. Step 2 independently
fits the technical-replicate observation model.

Input and clonotype definition
------------------------------
The script expects one file per subject x biological time x technical replicate.
By default filenames are parsed as:

    <subject>_<time>-<replica>.<extension>

where replica is 1 or 2. Supported formats are TSV/text, CSV, Parquet and
Feather/Arrow.

Each input file must contain at least:

    aaSeqCDR3
    readCount

For each repertoire:

    1. readCount is converted to numeric;
    2. non-finite and non-positive read-count rows are removed;
    3. rows with the same aaSeqCDR3 are collapsed by summing readCount;
    4. sequencing depth and relative clonotype frequencies are recalculated
       from the collapsed read counts.

Thus the analysis unit is the amino-acid CDR3 clonotype represented by
`aaSeqCDR3`. The script does not use an input `readFraction` column; frequencies
are recomputed as readCount divided by the total read count of the repertoire.

If more than one input file maps to the same subject x time x replica key, the
script stops rather than silently analysing duplicate representations of the
same technical replicate.

Observed repertoire summaries
------------------------------
For a repertoire containing N observed clonotypes with read counts c_i, total
sequencing depth is:

    D = sum_i c_i

and the observed relative frequency of clonotype i is:

    f_i = c_i / D.

The script reports D and N for every technical-replicate repertoire.

Primary descriptive power-law fit
---------------------------------
The high-frequency tail is described by a continuous Pareto model:

    p(f | f_min, gamma)
        = (gamma - 1) * f_min^(gamma - 1) * f^(-gamma),
          f >= f_min,

with gamma > 1.

For a fixed f_min, the maximum-likelihood estimator is:

    gamma_hat
        = 1 + n_tail / sum_i log(f_i / f_min),

where the sum runs over the n_tail clonotypes satisfying f_i >= f_min.

The lower tail boundary f_min is selected among observed frequency values by
minimizing the Kolmogorov-Smirnov distance:

    D_KS = sup_f |F_emp(f) - F_PL(f)|

between the empirical tail CDF and the fitted Pareto CDF. Candidate f_min values
can be thinned with `--candidate_step`; the default evaluates every distinct
observed candidate. A candidate tail must contain at least `--min_tail_size`
clonotypes (default 50).

The primary descriptive fit reports:

    pl_gamma_hat
    pl_xmin_hat
    pl_ks_distance
    pl_n_tail
    pl_tail_fraction
    pl_fit_success
    pl_fit_message

with:

    pl_tail_fraction = pl_n_tail / n_clonotypes.

This KS-selected fit is a repertoire-level structural descriptor. It is not an
upstream inclusion criterion and its fitted gamma is not inherited by Step 2.

Threshold-dependent tail-model comparison
-----------------------------------------
Because conclusions about a heavy-frequency tail can depend on the lower
boundary, the script additionally evaluates a fixed set of empirical threshold
quantiles. The default scan is:

    0.70, 0.75, 0.80, 0.85, 0.90, 0.95.

For each quantile q:

    f_min(q) = Q_q({f_i}),

and the same observed tail f_i >= f_min(q) is evaluated under:

    1. a continuous Pareto model;
    2. a lower-truncated log-normal (TLN) model.

For the Pareto model:

    log L_PL
        = n_tail * log(gamma - 1)
          + (gamma - 1) * n_tail * log(f_min)
          - gamma * sum_i log(f_i).

For the lower-truncated log-normal:

    p_TLN(f)
        = [1 / (f sigma sqrt(2*pi))]
          * exp[-(log(f) - mu)^2 / (2*sigma^2)]
          / [1 - Phi((log(f_min) - mu) / sigma)],

for f >= f_min.

The TLN parameters are estimated jointly by maximum likelihood including the
lower-truncation normalizing term. For numerical stability the implementation
uses the equivalent log-ratio coordinate:

    y = log(f / f_min) >= 0

with:

    p(y | beta, psi) proportional to exp(-beta*y - psi*y^2).

For psi > 0 this represents an interior lower-truncated log-normal model. The
boundary psi = 0 corresponds to the continuous Pareto family. The optimizer
therefore explicitly checks the TLN interior against the Pareto boundary.

Model-comparison interpretation
-------------------------------
The Pareto has one fitted parameter at fixed f_min; the TLN family has two.
Accordingly the primary descriptive model-comparison criterion is AICc, with
BIC retained as a sensitivity criterion:

    delta_aicc_tln_minus_pl
        = AICc_TLN - AICc_PL

    delta_bic_tln_minus_pl
        = BIC_TLN - BIC_PL.

Positive values favor the more parsimonious Pareto description.

Raw likelihood differences are also exported:

    delta_loglik_pl_minus_tln
        = logL_PL - logL_TLN.

They are DIAGNOSTIC ONLY. Because the implemented TLN family contains the
Pareto as a boundary case, a correctly optimized interior/boundary comparison
cannot use a larger Pareto raw likelihood as evidence for Pareto superiority;
model preference must account for the extra TLN parameter. Small positive raw
differences beyond numerical tolerance are therefore flagged for inspection of
optimization/boundary matching rather than interpreted biologically.

No likelihood-ratio significance test is performed by this step.

Filtering policy
----------------
No repertoire, replicate pair or clonotype is excluded on the basis of:

    sequencing depth;
    observed clonotype number;
    Pareto fit success;
    fitted gamma;
    KS distance;
    tail size;
    AICc/BIC preference;
    TLN convergence.

The script records descriptive fit failures and continues with the remaining
repertoires. Step 1 produces no PASS/FAIL QC classification.

Outputs
-------
The output directory contains:

1. `repertoire_characteristics.csv`

   One row per technical-replicate repertoire, including:

       subject, time, replica, pair_id, file
       depth, n_clonotypes
       pl_fit_success, pl_fit_message
       pl_gamma_hat, pl_xmin_hat, pl_ks_distance
       pl_n_tail, pl_tail_fraction
       threshold-scan summary statistics.

2. `repertoire_tail_model_scan.csv`

   One row per repertoire x threshold quantile, including:

       threshold_quantile, xmin_threshold, n_tail, tail_fraction
       Pareto gamma, log-likelihood and KS distance
       TLN convergence/boundary diagnostics and fitted parameters
       AIC, AICc and BIC for both models
       raw and per-tail-clonotype likelihood differences
       AICc/BIC model-preference indicators.

3. `repertoire_characterization_manifest.json`

   Machine-readable record of script version, input preprocessing,
   threshold-scan settings, model-comparison conventions and the explicit
   no-filtering / no-parameter-transfer policy.

Interpretive boundary
---------------------
This step supports statements about observed repertoire architecture and the
relative adequacy/parsimony of upper-tail descriptions. It does NOT establish
that the repertoire follows a power law over its full support, does NOT prove
that a power-law prior is the unique latent-state prior, and does NOT provide a
biological ground-truth clonotype abundance.

Recommended manuscript terminology is therefore:

    "observed repertoire architecture"
    "descriptive heavy-tail characterization"
    "continuous Pareto tail fit"
    "lower-truncated log-normal comparator"

rather than language implying that Step 1 validates the Step-2 latent prior.

Typical usage
-------------

    python3 1-repertoire_characterization.py \
        --data_dir ./dataset_longitudinal \
        --out_dir ./results/1-repertoire_characterization \
        --file_sep $'\\t' \
        --pattern '^(?P<subject>\\d+)_(?P<time>\\d+)-(?P<replica>[12])\\.tsv$'

Compatible with Python 3.9+; requires NumPy, pandas and SciPy.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import erfcx, log_ndtr


PATTERN_DEFAULT = (
    r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])"
    r"(?:\.(?:tsv|csv|parquet|pq|feather|arrow))?$"
)

PARQUET_EXTS = {".parquet", ".pq"}
FEATHER_EXTS = {".feather", ".arrow"}

SCRIPT_VERSION = "v2-production-descriptive-heavy-tail-2026-09-09"


@dataclass
class PowerLawFitResult:
    success: bool
    message: str
    n_total: int
    n_tail: int
    gamma_hat: float
    xmin_hat: float
    ks_distance: float
    tail_fraction: float


@dataclass
class TruncatedLognormalFitResult:
    success: bool
    message: str
    loglik: float
    mu: float
    sigma: float
    beta: float
    psi: float
    at_powerlaw_boundary: bool
    boundary_loglik: float
    interior_loglik: float
    optimizer_success: bool
    optimizer_message: str


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def validate_positive_1d(x: np.ndarray) -> np.ndarray:
    """Return sorted, finite, strictly positive values."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x > 0]
    return np.sort(x)


def continuous_gamma_mle(x_tail: np.ndarray, xmin: float) -> float:
    """Continuous Pareto MLE for gamma at fixed xmin."""
    x_tail = np.asarray(x_tail, dtype=float)

    if x_tail.size == 0 or not np.isfinite(xmin) or xmin <= 0:
        return np.nan

    denominator = np.sum(np.log(x_tail / xmin))
    if not np.isfinite(denominator) or denominator <= 0:
        return np.nan

    return float(1.0 + x_tail.size / denominator)


def powerlaw_cdf_continuous(
    x: np.ndarray,
    xmin: float,
    gamma: float,
) -> np.ndarray:
    """CDF of continuous Pareto distribution above xmin."""
    x = np.asarray(x, dtype=float)

    if (
        not np.isfinite(xmin)
        or not np.isfinite(gamma)
        or xmin <= 0
        or gamma <= 1.0
    ):
        return np.full_like(x, np.nan, dtype=float)

    return 1.0 - np.power(x / xmin, 1.0 - gamma)


def ks_distance_continuous(
    x_tail_sorted: np.ndarray,
    xmin: float,
    gamma: float,
) -> float:
    """Kolmogorov-Smirnov distance between empirical and fitted tail CDFs."""
    x_tail_sorted = np.asarray(x_tail_sorted, dtype=float)
    n = x_tail_sorted.size

    if n == 0:
        return np.nan

    empirical_cdf = np.arange(1, n + 1, dtype=float) / n
    model_cdf = powerlaw_cdf_continuous(
        x_tail_sorted,
        xmin=xmin,
        gamma=gamma,
    )

    if not np.all(np.isfinite(model_cdf)):
        return np.nan

    return float(np.max(np.abs(empirical_cdf - model_cdf)))


def estimate_xmin_by_ks(
    x: np.ndarray,
    min_tail_size: int = 50,
    candidate_step: int = 1,
) -> Tuple[float, float, float, int]:
    """
    Select xmin by minimizing the KS distance over observed candidate values.

    Returns
    -------
    xmin_hat, gamma_hat, ks_hat, n_tail
    """
    x = validate_positive_1d(x)

    if x.size < min_tail_size:
        return np.nan, np.nan, np.nan, 0

    unique_values = np.unique(x)
    step = max(1, int(candidate_step))
    candidates = unique_values[::step]

    best_xmin = np.nan
    best_gamma = np.nan
    best_ks = np.inf
    best_n_tail = 0

    for xmin in candidates:
        x_tail = x[x >= xmin]
        n_tail = int(x_tail.size)

        if n_tail < min_tail_size:
            continue

        gamma_hat = continuous_gamma_mle(x_tail, xmin)
        if not np.isfinite(gamma_hat) or gamma_hat <= 1.0:
            continue

        ks = ks_distance_continuous(
            x_tail_sorted=x_tail,
            xmin=xmin,
            gamma=gamma_hat,
        )
        if not np.isfinite(ks):
            continue

        if ks < best_ks:
            best_xmin = float(xmin)
            best_gamma = float(gamma_hat)
            best_ks = float(ks)
            best_n_tail = n_tail

    if not np.isfinite(best_xmin):
        return np.nan, np.nan, np.nan, 0

    return best_xmin, best_gamma, best_ks, best_n_tail


def fit_powerlaw_on_frequencies(
    freqs: np.ndarray,
    min_tail_size: int = 50,
    candidate_step: int = 1,
) -> PowerLawFitResult:
    """Fit a KS-selected continuous Pareto tail to one repertoire."""
    x = validate_positive_1d(freqs)

    if x.size < min_tail_size:
        return PowerLawFitResult(
            success=False,
            message="Not enough positive frequencies ({})".format(x.size),
            n_total=int(x.size),
            n_tail=0,
            gamma_hat=np.nan,
            xmin_hat=np.nan,
            ks_distance=np.nan,
            tail_fraction=np.nan,
        )

    xmin_hat, gamma_hat, ks_hat, n_tail = estimate_xmin_by_ks(
        x=x,
        min_tail_size=min_tail_size,
        candidate_step=candidate_step,
    )

    if not np.isfinite(xmin_hat):
        return PowerLawFitResult(
            success=False,
            message="Failed KS-based xmin estimation",
            n_total=int(x.size),
            n_tail=0,
            gamma_hat=np.nan,
            xmin_hat=np.nan,
            ks_distance=np.nan,
            tail_fraction=np.nan,
        )

    return PowerLawFitResult(
        success=True,
        message="ok",
        n_total=int(x.size),
        n_tail=int(n_tail),
        gamma_hat=float(gamma_hat),
        xmin_hat=float(xmin_hat),
        ks_distance=float(ks_hat),
        tail_fraction=float(n_tail / max(len(x), 1)),
    )


def powerlaw_tail_loglik(
    x_tail: np.ndarray,
    xmin: float,
    gamma: float,
) -> float:
    """Log-likelihood of the continuous Pareto model."""
    x_tail = np.asarray(x_tail, dtype=float)

    if (
        x_tail.size == 0
        or not np.isfinite(xmin)
        or not np.isfinite(gamma)
        or xmin <= 0
        or gamma <= 1.0
    ):
        return np.nan

    loglik = (
        x_tail.size * np.log(gamma - 1.0)
        + (gamma - 1.0) * x_tail.size * np.log(xmin)
        - gamma * np.sum(np.log(x_tail))
    )

    return float(loglik)


def _log_quadratic_tail_normalizer(
    u: float,
    v: float,
) -> float:
    """
    Stable log normalizer for exp(-u*t - v*t^2), t >= 0, v > 0.

    This parameterization is equivalent to a lower-truncated normal model in
    log-frequency space.  Its v -> 0, u > 0 boundary is the exponential model
    in log-ratio space, i.e. the continuous Pareto model in frequency space.
    """
    if not np.isfinite(u) or not np.isfinite(v) or v <= 0.0:
        return np.inf

    sqrt_v = np.sqrt(v)

    if u >= 0.0:
        a = u / (2.0 * sqrt_v)
        scaled_erfc = erfcx(a)
        if not np.isfinite(scaled_erfc) or scaled_erfc <= 0.0:
            return np.inf
        return float(
            0.5 * np.log(np.pi)
            - np.log(2.0)
            - 0.5 * np.log(v)
            + np.log(scaled_erfc)
        )

    # For u < 0, the log-CDF form avoids overflow of erfcx at large negative
    # arguments.
    return float(
        0.5 * (np.log(np.pi) - np.log(v))
        + (u * u) / (4.0 * v)
        + log_ndtr(-u / np.sqrt(2.0 * v))
    )


def fit_truncated_lognormal_tail_mle(
    x_tail: np.ndarray,
    xmin: float,
) -> TruncatedLognormalFitResult:
    """
    Jointly fit the lower-truncated log-normal tail by maximum likelihood.

    Let y = log(x / xmin) >= 0.  A lower-truncated log-normal tail can be
    written as

        p(y | beta, psi) proportional to exp(-beta*y - psi*y^2),

    with psi > 0.  This is numerically more stable than optimizing mu and
    sigma directly.  The Pareto model is the boundary psi = 0, beta > 0.

    The routine optimizes the interior TLN model and explicitly compares it
    with the Pareto boundary.  If the interior does not improve the
    log-likelihood beyond numerical tolerance, the fitted family is recorded
    as lying on the Pareto boundary.
    """
    x_tail = np.asarray(x_tail, dtype=float)
    x_tail = x_tail[np.isfinite(x_tail)]
    x_tail = x_tail[x_tail > 0.0]

    if x_tail.size < 3 or not np.isfinite(xmin) or xmin <= 0.0:
        return TruncatedLognormalFitResult(
            success=False,
            message="invalid tail or xmin",
            loglik=np.nan,
            mu=np.nan,
            sigma=np.nan,
            beta=np.nan,
            psi=np.nan,
            at_powerlaw_boundary=False,
            boundary_loglik=np.nan,
            interior_loglik=np.nan,
            optimizer_success=False,
            optimizer_message="not run",
        )

    z = np.log(x_tail)
    zmin = float(np.log(xmin))
    y = z - zmin

    if np.any(y < -1e-10):
        return TruncatedLognormalFitResult(
            success=False,
            message="tail contains observations below xmin",
            loglik=np.nan,
            mu=np.nan,
            sigma=np.nan,
            beta=np.nan,
            psi=np.nan,
            at_powerlaw_boundary=False,
            boundary_loglik=np.nan,
            interior_loglik=np.nan,
            optimizer_success=False,
            optimizer_message="not run",
        )

    y = np.maximum(y, 0.0)
    scale = float(np.mean(y))

    if not np.isfinite(scale) or scale <= 1e-15:
        return TruncatedLognormalFitResult(
            success=False,
            message="degenerate log-ratio tail",
            loglik=np.nan,
            mu=np.nan,
            sigma=np.nan,
            beta=np.nan,
            psi=np.nan,
            at_powerlaw_boundary=False,
            boundary_loglik=np.nan,
            interior_loglik=np.nan,
            optimizer_success=False,
            optimizer_message="not run",
        )

    # Rescale y so its empirical mean is one.  This keeps optimization in a
    # stable numerical range across thresholds and repertoires.
    t = y / scale
    n = int(t.size)
    sum_t = float(np.sum(t))
    sum_t2 = float(np.sum(t * t))
    constant = -float(np.sum(z)) - n * np.log(scale)

    # Pareto boundary: exponential rate 1/scale in y-space.
    boundary_loglik = float(constant - sum_t)

    def negative_loglik(theta: np.ndarray) -> float:
        u = float(theta[0])
        log_v = float(theta[1])
        v = float(np.exp(log_v))
        log_z_norm = _log_quadratic_tail_normalizer(u=u, v=v)

        if not np.isfinite(log_z_norm):
            return 1e300

        loglik = (
            constant
            - u * sum_t
            - v * sum_t2
            - n * log_z_norm
        )

        return float(-loglik) if np.isfinite(loglik) else 1e300

    starts = [
        (1.0, -12.0),
        (1.0, -6.0),
        (0.5, -3.0),
        (0.0, -1.0),
        (-0.5, -1.0),
        (2.0, -3.0),
    ]
    bounds = [(-100.0, 100.0), (-20.0, 10.0)]

    successful_results = []
    all_messages = []

    for start in starts:
        result = minimize(
            negative_loglik,
            x0=np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": 2000,
                "ftol": 1e-12,
                "gtol": 1e-8,
            },
        )
        all_messages.append(str(result.message))
        if result.success and np.isfinite(result.fun):
            successful_results.append(result)

    if not successful_results:
        return TruncatedLognormalFitResult(
            success=False,
            message="joint TLN optimization failed",
            loglik=np.nan,
            mu=np.nan,
            sigma=np.nan,
            beta=np.nan,
            psi=np.nan,
            at_powerlaw_boundary=False,
            boundary_loglik=boundary_loglik,
            interior_loglik=np.nan,
            optimizer_success=False,
            optimizer_message=" | ".join(all_messages),
        )

    best = min(successful_results, key=lambda result: float(result.fun))
    interior_loglik = float(-best.fun)

    # Differences much smaller than the information-criterion penalty are
    # treated as numerical equality with the Pareto boundary.
    tolerance = 1e-8 * max(1.0, abs(boundary_loglik))

    if interior_loglik <= boundary_loglik + tolerance:
        beta = float(1.0 / scale)
        return TruncatedLognormalFitResult(
            success=True,
            message="power-law boundary",
            loglik=boundary_loglik,
            mu=np.nan,
            sigma=np.nan,
            beta=beta,
            psi=0.0,
            at_powerlaw_boundary=True,
            boundary_loglik=boundary_loglik,
            interior_loglik=interior_loglik,
            optimizer_success=True,
            optimizer_message=str(best.message),
        )

    u_hat = float(best.x[0])
    v_hat = float(np.exp(best.x[1]))
    beta_hat = float(u_hat / scale)
    psi_hat = float(v_hat / (scale * scale))
    sigma_hat = float(1.0 / np.sqrt(2.0 * psi_hat))
    mu_hat = float(zmin - beta_hat / (2.0 * psi_hat))

    return TruncatedLognormalFitResult(
        success=True,
        message="interior truncated log-normal",
        loglik=interior_loglik,
        mu=mu_hat,
        sigma=sigma_hat,
        beta=beta_hat,
        psi=psi_hat,
        at_powerlaw_boundary=False,
        boundary_loglik=boundary_loglik,
        interior_loglik=interior_loglik,
        optimizer_success=True,
        optimizer_message=str(best.message),
    )


def information_criteria(
    loglik: float,
    n: int,
    k: int,
) -> Tuple[float, float, float]:
    """Return AIC, small-sample corrected AIC, and BIC."""
    if not np.isfinite(loglik) or n <= 0 or k <= 0:
        return np.nan, np.nan, np.nan

    aic = float(2.0 * k - 2.0 * loglik)
    if n > k + 1:
        aicc = float(aic + (2.0 * k * (k + 1.0)) / (n - k - 1.0))
    else:
        aicc = np.nan
    bic = float(k * np.log(n) - 2.0 * loglik)
    return aic, aicc, bic


def threshold_scan_tail_models(
    freqs: np.ndarray,
    subject: int,
    time: int,
    replica: int,
    file_name: str,
    threshold_quantiles: np.ndarray,
    min_tail_size: int,
) -> pd.DataFrame:
    """
    Compare Pareto and jointly fitted lower-truncated log-normal tails.

    Raw log-likelihood differences are retained as diagnostics.  Model
    preference is based on AICc (primary) and BIC (sensitivity), because the
    TLN family has one additional free parameter and contains the Pareto model
    on its boundary.
    """
    x = validate_positive_1d(freqs)
    rows: List[Dict[str, object]] = []

    if x.size == 0:
        return pd.DataFrame()

    for q in threshold_quantiles:
        q = float(q)
        xmin = float(np.quantile(x, q))
        x_tail = x[x >= xmin]
        n_tail = int(x_tail.size)

        if n_tail < min_tail_size:
            continue

        gamma_hat = continuous_gamma_mle(x_tail, xmin)

        if np.isfinite(gamma_hat) and gamma_hat > 1.0:
            ll_pl = powerlaw_tail_loglik(
                x_tail=x_tail,
                xmin=xmin,
                gamma=gamma_hat,
            )
            ks_pl = ks_distance_continuous(
                x_tail_sorted=x_tail,
                xmin=xmin,
                gamma=gamma_hat,
            )
        else:
            ll_pl = np.nan
            ks_pl = np.nan

        tln_fit = fit_truncated_lognormal_tail_mle(
            x_tail=x_tail,
            xmin=xmin,
        )
        ll_tln = tln_fit.loglik if tln_fit.success else np.nan

        pl_aic, pl_aicc, pl_bic = information_criteria(
            loglik=ll_pl,
            n=n_tail,
            k=1,
        )
        tln_aic, tln_aicc, tln_bic = information_criteria(
            loglik=ll_tln,
            n=n_tail,
            k=2,
        )

        if np.isfinite(ll_pl) and np.isfinite(ll_tln):
            delta_ll = float(ll_pl - ll_tln)
            mean_delta_ll = float(delta_ll / n_tail)
        else:
            delta_ll = np.nan
            mean_delta_ll = np.nan

        if np.isfinite(pl_aicc) and np.isfinite(tln_aicc):
            delta_aicc = float(tln_aicc - pl_aicc)
            mean_delta_aicc = float(delta_aicc / n_tail)
        else:
            delta_aicc = np.nan
            mean_delta_aicc = np.nan

        if np.isfinite(pl_bic) and np.isfinite(tln_bic):
            delta_bic = float(tln_bic - pl_bic)
            mean_delta_bic = float(delta_bic / n_tail)
        else:
            delta_bic = np.nan
            mean_delta_bic = np.nan

        rows.append(
            {
                "subject": subject,
                "time": time,
                "replica": replica,
                "pair_id": "{}_{}".format(subject, time),
                "file": file_name,
                "threshold_quantile": q,
                "xmin_threshold": xmin,
                "n_clonotypes": int(x.size),
                "n_tail": n_tail,
                "tail_fraction": float(n_tail / max(len(x), 1)),
                "pl_gamma_hat": (
                    float(gamma_hat)
                    if np.isfinite(gamma_hat)
                    else np.nan
                ),
                "pl_loglik": (
                    float(ll_pl)
                    if np.isfinite(ll_pl)
                    else np.nan
                ),
                "pl_mean_loglik_per_clonotype": (
                    float(ll_pl / n_tail)
                    if np.isfinite(ll_pl) and n_tail > 0
                    else np.nan
                ),
                "pl_ks_distance": (
                    float(ks_pl)
                    if np.isfinite(ks_pl)
                    else np.nan
                ),
                "tln_fit_success": bool(tln_fit.success),
                "tln_fit_message": tln_fit.message,
                "tln_optimizer_success": bool(tln_fit.optimizer_success),
                "tln_optimizer_message": tln_fit.optimizer_message,
                "tln_at_powerlaw_boundary": bool(
                    tln_fit.at_powerlaw_boundary
                ),
                "tln_boundary_loglik": (
                    float(tln_fit.boundary_loglik)
                    if np.isfinite(tln_fit.boundary_loglik)
                    else np.nan
                ),
                "tln_interior_loglik": (
                    float(tln_fit.interior_loglik)
                    if np.isfinite(tln_fit.interior_loglik)
                    else np.nan
                ),
                "tln_loglik": (
                    float(ll_tln)
                    if np.isfinite(ll_tln)
                    else np.nan
                ),
                "tln_mean_loglik_per_clonotype": (
                    float(ll_tln / n_tail)
                    if np.isfinite(ll_tln) and n_tail > 0
                    else np.nan
                ),
                "tln_mu": (
                    float(tln_fit.mu)
                    if np.isfinite(tln_fit.mu)
                    else np.nan
                ),
                "tln_sigma": (
                    float(tln_fit.sigma)
                    if np.isfinite(tln_fit.sigma)
                    else np.nan
                ),
                "tln_beta": (
                    float(tln_fit.beta)
                    if np.isfinite(tln_fit.beta)
                    else np.nan
                ),
                "tln_psi": (
                    float(tln_fit.psi)
                    if np.isfinite(tln_fit.psi)
                    else np.nan
                ),
                "pl_aic": pl_aic,
                "pl_aicc": pl_aicc,
                "pl_bic": pl_bic,
                "tln_aic": tln_aic,
                "tln_aicc": tln_aicc,
                "tln_bic": tln_bic,
                "delta_loglik_pl_minus_tln": delta_ll,
                "mean_delta_loglik_pl_minus_tln": mean_delta_ll,
                "delta_aicc_tln_minus_pl": delta_aicc,
                "mean_delta_aicc_tln_minus_pl": mean_delta_aicc,
                "delta_bic_tln_minus_pl": delta_bic,
                "mean_delta_bic_tln_minus_pl": mean_delta_bic,
                "pl_preferred_aicc": (
                    bool(delta_aicc > 0.0)
                    if np.isfinite(delta_aicc)
                    else False
                ),
                "pl_preferred_bic": (
                    bool(delta_bic > 0.0)
                    if np.isfinite(delta_bic)
                    else False
                ),
            }
        )

    return pd.DataFrame(rows)


def read_repertoire_file(
    file_path: Path,
    sep: str,
) -> pd.DataFrame:
    """Read one repertoire and return one positive read-count row per aaSeqCDR3.

    CSV files are read with comma separation; TSV/TXT-like files use ``sep``.
    Duplicate amino-acid CDR3 rows are collapsed by summing readCount so that
    repertoire depth, clonotype number and frequencies are defined on the same
    clonotype unit used by downstream Step 2.
    """
    suffix = file_path.suffix.lower()

    if suffix in PARQUET_EXTS:
        df = pd.read_parquet(file_path)
    elif suffix in FEATHER_EXTS:
        df = pd.read_feather(file_path)
    elif suffix == ".csv":
        df = pd.read_csv(file_path)
    else:
        df = pd.read_csv(file_path, sep=sep)

    required = {"aaSeqCDR3", "readCount"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "{}: missing columns {}. Expected at least {}".format(
                file_path.name,
                sorted(missing),
                sorted(required),
            )
        )

    out = df[["aaSeqCDR3", "readCount"]].copy()
    out["readCount"] = pd.to_numeric(
        out["readCount"],
        errors="coerce",
    )
    out = out[np.isfinite(out["readCount"])].copy()
    out = out[out["readCount"] > 0].copy()
    out = out.dropna(subset=["aaSeqCDR3"]).copy()

    # Canonical clonotype unit: amino-acid CDR3. This mirrors Step 2 and
    # prevents duplicated MiXCR rows from inflating richness or splitting the
    # frequency assigned to the same aaSeqCDR3.
    out = (
        out.groupby("aaSeqCDR3", as_index=False, sort=False)["readCount"]
        .sum()
    )
    out["readCount"] = out["readCount"].astype(np.int64)

    return out


def index_repertoires(
    data_dir: Path,
    pattern_str: str,
) -> List[Tuple[int, int, int, Path]]:
    """Index one unique input file per subject/time/replica key.

    Multiple matching files for the same technical replicate are treated as an
    input ambiguity and raise an error instead of generating duplicate
    repertoire rows.
    """
    pattern = re.compile(pattern_str)
    index: Dict[Tuple[int, int, int], Path] = {}

    for file_path in data_dir.iterdir():
        if not file_path.is_file():
            continue

        match = pattern.match(file_path.name)
        if match is None:
            continue

        key = (
            int(match.group("subject")),
            int(match.group("time")),
            int(match.group("replica")),
        )
        if key in index:
            raise ValueError(
                "Multiple files map to subject/time/replica {}: {} and {}".format(
                    key, index[key].name, file_path.name
                )
            )
        index[key] = file_path

    return [
        (subject, time, replica, index[(subject, time, replica)])
        for subject, time, replica in sorted(index)
    ]


def characterize_one_repertoire(
    file_path: Path,
    file_sep: str,
    min_tail_size: int,
    candidate_step: int,
) -> Tuple[PowerLawFitResult, int, int, np.ndarray]:
    """Read, summarize and fit one repertoire."""
    df = read_repertoire_file(
        file_path=file_path,
        sep=file_sep,
    )

    depth = int(df["readCount"].sum())
    n_clonotypes = int(len(df))

    if depth <= 0:
        fit = PowerLawFitResult(
            success=False,
            message="Zero depth",
            n_total=n_clonotypes,
            n_tail=0,
            gamma_hat=np.nan,
            xmin_hat=np.nan,
            ks_distance=np.nan,
            tail_fraction=np.nan,
        )
        return fit, depth, n_clonotypes, np.array([], dtype=float)

    freqs = (
        df["readCount"].astype(float)
        / float(depth)
    ).to_numpy(dtype=float)

    fit = fit_powerlaw_on_frequencies(
        freqs=freqs,
        min_tail_size=min_tail_size,
        candidate_step=candidate_step,
    )

    return fit, depth, n_clonotypes, freqs


def summarize_threshold_scan(
    scan_df: pd.DataFrame,
) -> pd.DataFrame:
    """Create one threshold-scan summary row per repertoire."""
    summary_columns = [
        "subject",
        "time",
        "replica",
        "pair_id",
        "file",
        "scan_n_thresholds",
        "scan_tln_fit_success_fraction",
        "scan_tln_boundary_fraction",
        "scan_mean_delta_loglik",
        "scan_median_delta_loglik",
        "scan_mean_delta_loglik_per_clonotype",
        "scan_median_delta_loglik_per_clonotype",
        "scan_mean_delta_aicc_tln_minus_pl",
        "scan_median_delta_aicc_tln_minus_pl",
        "scan_frac_thresholds_pl_preferred_aicc",
        "scan_mean_delta_bic_tln_minus_pl",
        "scan_median_delta_bic_tln_minus_pl",
        "scan_frac_thresholds_pl_preferred_bic",
        "scan_median_tail_fraction",
    ]

    if scan_df.empty:
        return pd.DataFrame(columns=summary_columns)

    summary = (
        scan_df
        .groupby(
            ["subject", "time", "replica", "pair_id", "file"],
            as_index=False,
        )
        .agg(
            scan_n_thresholds=("threshold_quantile", "size"),
            scan_tln_fit_success_fraction=(
                "tln_fit_success",
                lambda x: float(pd.Series(x).mean()),
            ),
            scan_tln_boundary_fraction=(
                "tln_at_powerlaw_boundary",
                lambda x: float(pd.Series(x).mean()),
            ),
            scan_mean_delta_loglik=(
                "delta_loglik_pl_minus_tln",
                "mean",
            ),
            scan_median_delta_loglik=(
                "delta_loglik_pl_minus_tln",
                "median",
            ),
            scan_mean_delta_loglik_per_clonotype=(
                "mean_delta_loglik_pl_minus_tln",
                "mean",
            ),
            scan_median_delta_loglik_per_clonotype=(
                "mean_delta_loglik_pl_minus_tln",
                "median",
            ),
            scan_mean_delta_aicc_tln_minus_pl=(
                "delta_aicc_tln_minus_pl",
                "mean",
            ),
            scan_median_delta_aicc_tln_minus_pl=(
                "delta_aicc_tln_minus_pl",
                "median",
            ),
            scan_frac_thresholds_pl_preferred_aicc=(
                "pl_preferred_aicc",
                lambda x: float(pd.Series(x).mean()),
            ),
            scan_mean_delta_bic_tln_minus_pl=(
                "delta_bic_tln_minus_pl",
                "mean",
            ),
            scan_median_delta_bic_tln_minus_pl=(
                "delta_bic_tln_minus_pl",
                "median",
            ),
            scan_frac_thresholds_pl_preferred_bic=(
                "pl_preferred_bic",
                lambda x: float(pd.Series(x).mean()),
            ),
            scan_median_tail_fraction=(
                "tail_fraction",
                "median",
            ),
        )
    )

    # Explicit aliases preserve the sign convention in the column name while
    # retaining the historical summary fields for backward compatibility.
    summary["scan_mean_delta_loglik_pl_minus_tln"] = summary["scan_mean_delta_loglik"]
    summary["scan_median_delta_loglik_pl_minus_tln"] = summary["scan_median_delta_loglik"]
    summary["scan_mean_delta_loglik_per_clonotype_pl_minus_tln"] = summary[
        "scan_mean_delta_loglik_per_clonotype"
    ]
    summary["scan_median_delta_loglik_per_clonotype_pl_minus_tln"] = summary[
        "scan_median_delta_loglik_per_clonotype"
    ]

    return summary


def write_manifest(
    out_dir: Path,
    args: argparse.Namespace,
    threshold_quantiles: np.ndarray,
    n_repertoires: int,
) -> Path:
    """Write the production contract for descriptive Step-1 outputs."""
    payload = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "pipeline_step": 1,
        "analysis_role": "observed_repertoire_architecture_and_descriptive_heavy_tail_characterization",
        "n_repertoires_indexed": int(n_repertoires),
        "input_processing": {
            "clonotype_key": "aaSeqCDR3",
            "count_field": "readCount",
            "drop_nonfinite_counts": True,
            "drop_nonpositive_counts": True,
            "collapse_duplicate_aaSeqCDR3": "sum_readCount",
            "recompute_relative_frequency_from_readCount": True,
            "input_readFraction_used": False,
        },
        "primary_tail_fit": {
            "family": "continuous_Pareto",
            "xmin_selection": "KS_minimization_over_observed_frequency_candidates",
            "min_tail_size": int(args.min_tail_size),
            "candidate_step": int(args.candidate_step),
        },
        "threshold_scan": {
            "quantiles": [float(x) for x in threshold_quantiles],
            "min_tail_size": int(args.scan_min_tail_size),
            "comparator": "lower_truncated_lognormal_with_Pareto_boundary",
            "primary_preference_criterion": "AICc",
            "secondary_preference_criterion": "BIC",
            "raw_delta_loglik_definition": "logL_PL_minus_logL_TLN",
            "raw_delta_loglik_interpretation": "diagnostic_only_not_formal_model_preference",
            "formal_likelihood_ratio_test": False,
        },
        "pipeline_policy": {
            "qc_filtering": False,
            "sample_pass_fail_classification": False,
            "clonotype_filtering_by_tail_fit": False,
            "step1_gamma_passed_to_step2": False,
            "step1_parameterizes_step2_prior": False,
            "step1_validates_step2_prior": False,
        },
        "outputs": [
            "repertoire_characteristics.csv",
            "repertoire_tail_model_scan.csv",
            "repertoire_characterization_manifest.json",
        ],
    }
    path = Path(out_dir) / "repertoire_characterization_manifest.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Repertoire-level characterization and heavy-tail fitting "
            "without QC or filtering."
        )
    )

    parser.add_argument(
        "--data_dir",
        type=Path,
        required=True,
        help="Directory containing repertoire replicate files.",
    )
    parser.add_argument(
        "--out_dir",
        type=Path,
        required=True,
        help="Directory for repertoire characterization outputs.",
    )
    parser.add_argument(
        "--file_sep",
        type=str,
        default="\t",
        help="Input separator for text files. Default: tab.",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default=PATTERN_DEFAULT,
        help=(
            "Filename regex containing named subject, time and replica groups."
        ),
    )
    parser.add_argument(
        "--min_tail_size",
        type=int,
        default=50,
        help=(
            "Minimum number of clonotypes required in the KS-selected "
            "power-law tail. Default: 50."
        ),
    )
    parser.add_argument(
        "--candidate_step",
        type=int,
        default=1,
        help=(
            "Step used to thin observed candidate xmin values. Default: 1."
        ),
    )
    parser.add_argument(
        "--threshold_quantiles",
        type=str,
        default="0.70,0.75,0.80,0.85,0.90,0.95",
        help=(
            "Comma-separated empirical quantiles used for threshold-dependent "
            "tail-model comparison."
        ),
    )
    parser.add_argument(
        "--scan_min_tail_size",
        type=int,
        default=50,
        help=(
            "Minimum tail size for threshold-dependent model comparisons. "
            "Default: 50."
        ),
    )

    return parser


def main() -> None:
    args = build_argparser().parse_args()

    if not args.data_dir.exists():
        raise FileNotFoundError(
            "data_dir not found: {}".format(args.data_dir)
        )

    if args.min_tail_size < 2:
        raise ValueError("--min_tail_size must be >= 2")

    if args.scan_min_tail_size < 3:
        raise ValueError("--scan_min_tail_size must be >= 3")

    if args.candidate_step < 1:
        raise ValueError("--candidate_step must be >= 1")

    threshold_quantiles = np.unique(np.array(
        [
            float(value.strip())
            for value in str(args.threshold_quantiles).split(",")
            if value.strip()
        ],
        dtype=float,
    ))

    if (
        threshold_quantiles.size == 0
        or np.any(
            (threshold_quantiles <= 0)
            | (threshold_quantiles >= 1)
        )
    ):
        raise ValueError(
            "--threshold_quantiles must contain values strictly between 0 and 1"
        )

    out_dir = ensure_dir(args.out_dir)

    repertoires = index_repertoires(
        data_dir=args.data_dir,
        pattern_str=args.pattern,
    )

    print(
        "[Step1 {}] Found {} matching repertoire replicate files.".format(
            SCRIPT_VERSION, len(repertoires)
        )
    )

    if not repertoires:
        examples = sorted([p.name for p in args.data_dir.iterdir() if p.is_file()])[:10]
        raise ValueError(
            "No repertoire files matched --pattern. Pattern: {}. First files: {}".format(
                args.pattern, examples
            )
        )

    repertoire_rows: List[Dict[str, object]] = []
    scan_frames: List[pd.DataFrame] = []

    for subject, time, replica, file_path in repertoires:
        pair_id = "{}_{}".format(subject, time)

        print(
            "Processing {}, replicate {}: {}".format(
                pair_id,
                replica,
                file_path.name,
            )
        )

        try:
            fit, depth, n_clonotypes, freqs = characterize_one_repertoire(
                file_path=file_path,
                file_sep=args.file_sep,
                min_tail_size=args.min_tail_size,
                candidate_step=args.candidate_step,
            )

            repertoire_rows.append(
                {
                    "subject": subject,
                    "time": time,
                    "replica": replica,
                    "pair_id": pair_id,
                    "file": file_path.name,
                    "depth": depth,
                    "n_clonotypes": n_clonotypes,
                    "pl_fit_success": fit.success,
                    "pl_fit_message": fit.message,
                    "pl_gamma_hat": fit.gamma_hat,
                    "pl_xmin_hat": fit.xmin_hat,
                    "pl_ks_distance": fit.ks_distance,
                    "pl_n_tail": fit.n_tail,
                    "pl_tail_fraction": fit.tail_fraction,
                }
            )

            freqs_pos = validate_positive_1d(freqs)

            if freqs_pos.size >= args.scan_min_tail_size:
                scan_df = threshold_scan_tail_models(
                    freqs=freqs_pos,
                    subject=subject,
                    time=time,
                    replica=replica,
                    file_name=file_path.name,
                    threshold_quantiles=threshold_quantiles,
                    min_tail_size=args.scan_min_tail_size,
                )

                if not scan_df.empty:
                    scan_frames.append(scan_df)

        except Exception as exc:
            repertoire_rows.append(
                {
                    "subject": subject,
                    "time": time,
                    "replica": replica,
                    "pair_id": pair_id,
                    "file": file_path.name,
                    "depth": np.nan,
                    "n_clonotypes": np.nan,
                    "pl_fit_success": False,
                    "pl_fit_message": "ERROR: {}".format(exc),
                    "pl_gamma_hat": np.nan,
                    "pl_xmin_hat": np.nan,
                    "pl_ks_distance": np.nan,
                    "pl_n_tail": np.nan,
                    "pl_tail_fraction": np.nan,
                }
            )

            print(
                "WARNING: failed to characterize {}: {}".format(
                    file_path.name,
                    exc,
                )
            )

    repertoire_df = pd.DataFrame(repertoire_rows)

    if scan_frames:
        scan_df = pd.concat(
            scan_frames,
            ignore_index=True,
        )
    else:
        scan_df = pd.DataFrame()

    scan_summary = summarize_threshold_scan(scan_df)

    if not repertoire_df.empty:
        if not scan_summary.empty:
            repertoire_df = repertoire_df.merge(
                scan_summary,
                on=[
                    "subject",
                    "time",
                    "replica",
                    "pair_id",
                    "file",
                ],
                how="left",
                validate="one_to_one",
            )
        else:
            for column in [
                "scan_n_thresholds",
                "scan_tln_fit_success_fraction",
                "scan_tln_boundary_fraction",
                "scan_mean_delta_loglik",
                "scan_median_delta_loglik",
                "scan_mean_delta_loglik_per_clonotype",
                "scan_median_delta_loglik_per_clonotype",
                "scan_mean_delta_aicc_tln_minus_pl",
                "scan_median_delta_aicc_tln_minus_pl",
                "scan_frac_thresholds_pl_preferred_aicc",
                "scan_mean_delta_bic_tln_minus_pl",
                "scan_median_delta_bic_tln_minus_pl",
                "scan_frac_thresholds_pl_preferred_bic",
                "scan_median_tail_fraction",
                "scan_mean_delta_loglik_pl_minus_tln",
                "scan_median_delta_loglik_pl_minus_tln",
                "scan_mean_delta_loglik_per_clonotype_pl_minus_tln",
                "scan_median_delta_loglik_per_clonotype_pl_minus_tln",
            ]:
                repertoire_df[column] = np.nan

        repertoire_df = repertoire_df.sort_values(
            ["subject", "time", "replica", "file"]
        ).reset_index(drop=True)

    if not scan_df.empty:
        scan_df = scan_df.sort_values(
            [
                "subject",
                "time",
                "replica",
                "threshold_quantile",
            ]
        ).reset_index(drop=True)

    characteristics_path = (
        out_dir
        / "repertoire_characteristics.csv"
    )
    scan_path = (
        out_dir
        / "repertoire_tail_model_scan.csv"
    )

    repertoire_df.to_csv(
        characteristics_path,
        index=False,
    )
    scan_df.to_csv(
        scan_path,
        index=False,
    )

    manifest_path = write_manifest(
        out_dir=out_dir,
        args=args,
        threshold_quantiles=threshold_quantiles,
        n_repertoires=len(repertoires),
    )

    print("\nSaved repertoire characterization outputs:")
    print(" - {}".format(characteristics_path))
    print(" - {}".format(scan_path))
    print(" - {}".format(manifest_path))

    if not repertoire_df.empty:
        n_total = int(len(repertoire_df))
        n_fit_success = int(
            repertoire_df["pl_fit_success"]
            .fillna(False)
            .astype(bool)
            .sum()
        )

        print(
            "\nPrimary power-law fits successful: {}/{}".format(
                n_fit_success,
                n_total,
            )
        )

    if not scan_df.empty:
        n_scan = int(len(scan_df))
        n_tln_success = int(
            scan_df["tln_fit_success"]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        n_pl_aicc = int(
            scan_df["pl_preferred_aicc"]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        n_pl_bic = int(
            scan_df["pl_preferred_bic"]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        max_raw_delta = float(
            pd.to_numeric(
                scan_df["delta_loglik_pl_minus_tln"],
                errors="coerce",
            ).max()
        )

        print(
            "Joint TLN fits successful: {}/{}".format(
                n_tln_success,
                n_scan,
            )
        )
        print(
            "PL preferred by AICc: {}/{}".format(
                n_pl_aicc,
                n_scan,
            )
        )
        print(
            "PL preferred by BIC: {}/{}".format(
                n_pl_bic,
                n_scan,
            )
        )
        print(
            "Maximum raw (logL_PL - logL_TLN): {:.6g}".format(
                max_raw_delta,
            )
        )

        if np.isfinite(max_raw_delta) and max_raw_delta > 1e-5:
            print(
                "WARNING: positive raw likelihood difference detected; "
                "inspect TLN convergence and boundary matching."
            )


if __name__ == "__main__":
    main()
