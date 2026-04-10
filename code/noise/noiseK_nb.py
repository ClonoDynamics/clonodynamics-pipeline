#!/usr/bin/env python3
"""
noiseK_nb.py (v2)

Negative Binomial noise model with power-law prior for TCR repertoire data.
This module provides likelihood-based inference of sampling and experimental
noise from pairs of technical replicates.

This version removes all plotting/visualization dependencies and is safe to
use in headless/HPC environments. Numerical behavior and API are identical
to the original implementation.

Requirements:
- Python >= 3.9
- numpy
- pandas
- scipy (optimize, special)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp


@dataclass
class FitResult:
    success: bool
    message: str
    log_likelihood: float
    params: Dict[str, float]
    per_clone_logP: pd.DataFrame


def _log_nb_pmf(k: np.ndarray, mu: np.ndarray, r: float) -> np.ndarray:
    return (
        gammaln(k + r)
        - gammaln(r)
        - gammaln(k + 1)
        + r * np.log(r / (r + mu))
        + k * np.log(mu / (r + mu))
    )


def fit_noiseK_nb_powerlaw(
    merged: pd.DataFrame,
    depths: List[int],
    count_cols: List[str],
    gamma_init: float = 1.6,
    k_init: float = 50.0,
    grid_size: int = 500,
) -> FitResult:
    c1 = merged[count_cols[0]].to_numpy(dtype=int)
    c2 = merged[count_cols[1]].to_numpy(dtype=int)

    depth1, depth2 = depths
    N_total = 0.5 * (depth1 + depth2)

    fmin = 1.0 / max(N_total, 1.0)

    f_grid = np.logspace(np.log10(fmin), 0, grid_size)

    def neg_log_likelihood(theta):
        gamma, kappa = theta
        if gamma <= 0 or kappa <= 0:
            return np.inf

        log_prior = -gamma * np.log(f_grid)
        log_prior -= logsumexp(log_prior)

        mu1 = f_grid * depth1
        mu2 = f_grid * depth2

        logP1 = _log_nb_pmf(c1[:, None], mu1[None, :], kappa)
        logP2 = _log_nb_pmf(c2[:, None], mu2[None, :], kappa)

        log_joint = logP1 + logP2 + log_prior
        log_like = logsumexp(log_joint, axis=1)

        return -np.sum(log_like)

    res = minimize(
        neg_log_likelihood,
        x0=np.array([gamma_init, k_init]),
        method="L-BFGS-B",
        bounds=[(1e-6, None), (1e-6, None)],
    )

    gamma_hat, k_hat = res.x
    ll = -res.fun

    log_prior = -gamma_hat * np.log(f_grid)
    log_prior -= logsumexp(log_prior)

    mu1 = f_grid * depth1
    mu2 = f_grid * depth2

    logP1 = _log_nb_pmf(c1[:, None], mu1[None, :], k_hat)
    logP2 = _log_nb_pmf(c2[:, None], mu2[None, :], k_hat)

    log_joint = logP1 + logP2 + log_prior
    per_clone_logP = logsumexp(log_joint, axis=1)

    out = merged.copy()
    out["logP"] = per_clone_logP

    params = {
        "gamma": float(gamma_hat),
        "k": float(k_hat),
        "fmin": float(fmin),
        "N_total_mean": float(N_total),
        "N_total_min": float(min(depth1, depth2)),
        "N_total_max": float(max(depth1, depth2)),
    }

    return FitResult(
        success=res.success,
        message=res.message,
        log_likelihood=float(ll),
        params=params,
        per_clone_logP=out,
    )
