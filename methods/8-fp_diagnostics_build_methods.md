# Supplementary Methods — Fokker–Planck closed-vs-open diagnostics in log-frequency space (`8-fp_diagnostics_build.py`)

## Overview

This script constructs diagnostic grids for comparing empirical clonotype dynamics with the stationary solution of a one-dimensional Fokker–Planck equation in log-frequency space,

\[
x = \ln f.
\]

Its purpose is not to fit new dynamical models, but to combine previously estimated drift and diffusion functions with the empirical distribution of clonotype abundances and thereby evaluate whether the observed dynamics are compatible with a closed stationary process or instead show signatures of open, non-conservative behavior.

It reads:

- a drift grid \(b(x)\) from a fitted drift model,
- a binned diffusion estimate \(D(x)\),
- an empirical trajectory table used to reconstruct the observed density \(p_{\mathrm{emp}}(x)\),

and returns a unified grid containing:

- the model drift \(b(x)\),
- the interpolated diffusion \(D(x)\),
- the empirical density \(p_{\mathrm{emp}}(x)\),
- the closed stationary density \(p_{\mathrm{closed}}(x)\),
- the density ratio \(p_{\mathrm{emp}}/p_{\mathrm{closed}}\),
- the empirical probability current \(J(x)\),
- the source term \(S(x)=dJ/dx\).

This documentation is based on the uploaded script `8-fp_diagnostics_build.py`. fileciteturn9file0

## Conceptual role in the workflow

Upstream steps provide:

1. a fitted deterministic drift function \(b(x)\),
2. a diffusion estimate \(D(x)\) from residual fluctuations,
3. a trajectory table giving empirical occupancy of state space.

This script combines these ingredients to test two related questions:

1. **Closed stationary consistency**  
   If the repertoire behaved as a stationary one-dimensional diffusion with drift \(b(x)\) and diffusion \(D(x)\), what stationary density would it imply?

2. **Open-process diagnostics**  
   Given the empirical density \(p_{\mathrm{emp}}(x)\), what probability current \(J(x)\) and source term \(S(x)\) would be required to sustain it under the same drift–diffusion field?

In this sense, the script bridges

\[
\text{drift fit} + \text{diffusion estimate} + \text{empirical occupancy}
\;\longrightarrow\;
\text{stationary FP diagnostics}.
\]

## Inputs

The script requires three CSV files.

### 1. `drift_function_grid.csv`

This file contains a grid of \(x\) values and one or more fitted drift functions. The relevant model-specific drift column must be one of:

- `b_OU_linear`
- `b_tanh_saturating`
- `b_rational_saturating`
- `b_hinge_plateau_smooth`

The script also expects filtering variables such as:

- `dt`
- `obs_class`
- optionally `weight_col`
- optionally `weight_cap_q`

The drift grid is filtered to one specific dynamical slice defined by command-line arguments.

### 2. `diffusion_binned.csv`

This file contains binned estimates of the variance of residual fluctuations as a function of abundance. It is filtered by:

- `dt`
- `obs_class`
- `model`
- optionally `weight_col`
- optionally `weight_cap_q`

The abundance coordinate is identified from one of:

- `x0_mid`
- `x_mid`
- `x`
- `x_center`

The script then extracts one variance proxy from the following priority order:

1. `var_mad`
2. `mean_resid2`
3. `mad` (converted to `mad^2`)

### 3. `trajectories_long.csv`

This file provides empirical occupancy of state space. The script extracts the log-frequency coordinate from:

1. `log_freq`
2. `x`
3. `log(freq)` if only `freq` is available

Optionally, it can restrict the empirical density to rows with pointwise `observable == True`.

## State-space representation

All calculations are performed in log-frequency space:

\[
x = \ln f.
\]

This choice is natural because:

- clone frequencies span several orders of magnitude,
- upstream drift and diffusion fits are already expressed in \(x\)-space,
- stationary Fokker–Planck equations are more interpretable when multiplicative abundance changes become additive.

The diagnostic grid is therefore defined over a one-dimensional axis \(x\), inherited from the fitted drift grid.

## Filtering logic

The script identifies one consistent dynamical slice across drift and diffusion inputs by matching:

- a single \(dt\),
- a single `obs_class`,
- a single drift model,
- and, if specified, exact weight metadata.

Formally, the diagnostics are computed for one selected combination

\[
(dt,\ \mathrm{obs\_class},\ \mathrm{model},\ \mathrm{weighting}).
\]

This is important because the drift and diffusion estimates may differ substantially across time lag, detectability class, weighting scheme, and fitted drift model. The script therefore avoids mixing incompatible slices.

## Drift extraction

After filtering, the script extracts:

- a grid \(x_1,\dots,x_n\),
- a model-specific drift function \(b(x)\).

The drift grid is sorted by \(x\), and a minimum number of points is required in order to compute stable numerical derivatives later.

If the requested drift model is \(m\), then the script uses the corresponding column:

\[
b(x) = b_m(x).
\]

No further smoothing of the drift is applied inside this script; it assumes that the upstream fit has already produced a suitable grid representation.

## Diffusion construction

The diffusion input is not necessarily provided directly as \(D(x)\), but rather as a variance proxy for displacements over lag \(dt\). The script converts this to diffusion using

\[
D(x) = \frac{\mathrm{Var}(\Delta x \mid x, dt)}{2\,dt}.
\]

Depending on the available column, the variance proxy is defined as:

### Option 1 — `var_mad`

A robust variance estimate already computed upstream:

\[
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathrm{var\_mad}(x).
\]

### Option 2 — `mean_resid2`

The conditional mean squared residual:

\[
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathbb{E}[r^2 \mid x].
\]

### Option 3 — `mad`

If only the median absolute deviation is available, the script uses

\[
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathrm{mad}(x)^2.
\]

After constructing \(D(x)\) on the diffusion bins, the script interpolates it onto the drift grid \(x_{\mathrm{grid}}\). If interpolation leaves missing values at the edges, it fills them by one-dimensional interpolation across the available finite points.

Thus the final diffusion field is

\[
D_{\mathrm{on}}(x_i), \qquad i=1,\dots,n,
\]

defined on the same grid as the drift.

## Empirical density \(p_{\mathrm{emp}}(x)\)

The empirical density is constructed from the trajectory table by histogramming state-space occupancy onto the drift grid.

### Choice of empirical variable

If not overridden by `--x_col_traj`, the script uses:

- `log_freq`, if present,
- else `x`,
- else \(\ln(\mathrm{freq})\).

Only finite values are retained.

### Optional observability filter

If `--traj_observable_only` is set, the empirical density is built only from rows satisfying

\[
\mathrm{observable} = \mathrm{True}.
\]

This is a pointwise filter on state occupancy and should not be confused with transition classes such as TT, TF, FT, FF.

### Histogram construction on a predefined grid

The script defines histogram bin edges from the drift grid midpoints. Let the grid be

\[
x_1 < x_2 < \cdots < x_n.
\]

Then the histogram edges are built as midpoint boundaries between adjacent grid points, with extrapolated outer edges at the ends. Counts are converted to a density by dividing by total counts and bin widths.

Thus the empirical density satisfies approximately

\[
\int p_{\mathrm{emp}}(x)\,dx = 1.
\]

### Optional smoothing

The empirical density can be smoothed in index-space by Gaussian convolution with standard deviation `smooth_p_sigma_pts` measured in grid points. After smoothing, the density is re-normalized by numerical integration.

## Closed stationary density \(p_{\mathrm{closed}}(x)\)

For a one-dimensional Fokker–Planck equation with drift \(b(x)\) and diffusion \(D(x)\), the stationary zero-current solution satisfies

\[
J(x)=0.
\]

Under the convention used in this script, the resulting closed stationary density is

\[
p_{\mathrm{closed}}(x)
\propto
\frac{1}{D(x)}
\exp\!\left(\int^x \frac{b(y)}{D(y)}\,dy\right).
\]

### Numerical implementation

The script computes:

1. a safe diffusion field

\[
D_{\mathrm{safe}}(x) = \max(D(x), \varepsilon),
\]

with \(\varepsilon = 10^{-12}\),

2. the integrand

\[
I'(x) = \frac{b(x)}{D_{\mathrm{safe}}(x)},
\]

3. its cumulative trapezoidal integral

\[
I(x) = \int^x \frac{b(y)}{D_{\mathrm{safe}}(y)}\,dy,
\]

4. the unnormalized density

\[
p_{\mathrm{unn}}(x)
=
\frac{1}{D_{\mathrm{safe}}(x)}
\exp(I(x)).
\]

For numerical stability, the maximum of \(I(x)\) is subtracted before exponentiation. Finally, the density is normalized using trapezoidal integration:

\[
p_{\mathrm{closed}}(x)
=
\frac{p_{\mathrm{unn}}(x)}
{\int p_{\mathrm{unn}}(x)\,dx}.
\]

## Empirical current \(J(x)\)

The script next computes the current implied by the empirical density \(p_{\mathrm{emp}}(x)\) under the same drift and diffusion fields:

\[
J(x) = b(x)\,p_{\mathrm{emp}}(x) - \frac{d}{dx}\bigl[D(x)\,p_{\mathrm{emp}}(x)\bigr].
\]

This quantity measures local imbalance between deterministic transport and diffusive redistribution.

The derivative is computed numerically by `np.gradient` on the drift grid. Because numerical derivatives amplify noise, the raw current is smoothed by Gaussian convolution with width `smooth_J_sigma_pts` measured in grid points:

\[
J_s(x) = \mathcal{G}_\sigma * J(x).
\]

The smoothed current is what is written to output as `J_emp`.

## Source term \(S(x)\)

The source term is defined as the derivative of the smoothed current:

\[
S(x) = \frac{dJ_s(x)}{dx}.
\]

In a strictly closed stationary system with zero current, one expects

\[
J(x) = 0, \qquad S(x)=0.
\]

Nonzero \(J(x)\) or \(S(x)\) indicates that the empirical density cannot be explained purely as a closed stationary balance under the supplied drift and diffusion fields.

As with the current, the raw derivative is smoothed again before saving:

\[
S_s(x) = \mathcal{G}_\sigma * \frac{dJ_s(x)}{dx}.
\]

The saved column `S_emp` corresponds to this smoothed source profile.

## Density ratio

The script also computes the pointwise ratio

\[
\frac{p_{\mathrm{emp}}(x)}{p_{\mathrm{closed}}(x)}.
\]

This ratio directly identifies abundance regions where the empirical occupancy is enriched relative to the closed stationary prediction, depleted relative to it, or broadly consistent with it.

## Summary metrics: JS and KS

To compare the empirical and closed stationary densities globally, the script computes two summary metrics.

### Jensen–Shannon divergence

The densities are converted to discrete probability masses on the grid by multiplying by local grid spacing \(dx\), then normalized. The Jensen–Shannon divergence is

\[
\mathrm{JS}(P,Q)
=
\frac{1}{2}\mathrm{KL}(P\|M)
+
\frac{1}{2}\mathrm{KL}(Q\|M),
\qquad
M=\frac{P+Q}{2}.
\]

### Kolmogorov–Smirnov distance

The script also computes

\[
\mathrm{KS}(P,Q)
=
\max_x
\left|
F_P(x) - F_Q(x)
\right|.
\]

Both metrics are stored in the summary output.

## Output files

The script writes three files.

### 1. `fp_<tag>.grid.csv`

This is the main diagnostic grid. It contains:

- `x`
- `b`
- `D`
- `D_var_source`
- `p_emp`
- `p_closed`
- `ratio_pemp_over_pclosed`
- `J_emp`
- `S_emp`

### 2. `fp_<tag>.summary.json`

This JSON file stores metadata and summary statistics, including:

- file paths
- selected `dt`, `obs_class`, and model
- weighting filters
- observability settings
- smoothing parameters
- variance source
- JS divergence
- KS distance
- path to the diagnostic grid CSV

### 3. `fp_<tag>.summary.txt`

This is a plain-text key=value version of the same summary metadata.

## Tag construction and reproducibility

The output tag is built from the requested parameter combination, typically including:

- \(dt\)
- `obs_class`
- drift model
- optional weight metadata
- optional observability restriction

This ensures that diagnostic grids from different configurations do not overwrite one another and remain traceable to the exact analysis settings used.

## Numerical conventions and safeguards

Several safeguards are built into the implementation.

### Safe diffusion floor

Diffusion is clipped below by a small positive constant:

\[
D(x) \ge 10^{-12}.
\]

### Interpolation gap filling

If diffusion interpolation leaves missing values on the drift grid, the script fills them by interpolation from nearby finite points, provided enough finite support exists.

### Smoothing in grid-index space

Gaussian smoothing is performed in index space, not in physical \(x\)-space. Therefore the smoothing width is measured in grid points, not in units of \(\ln f\).

### Normalization by trapezoidal integration

All densities are normalized by numerical integration using the trapezoidal rule.

## Statistical interpretation

This script does not fit a new Fokker–Planck model. Instead, it asks whether the previously fitted drift and diffusion functions are jointly consistent with the observed abundance distribution.

### Closed-process interpretation

If

\[
p_{\mathrm{emp}}(x) \approx p_{\mathrm{closed}}(x),
\qquad
J(x) \approx 0,
\qquad
S(x) \approx 0,
\]

then the observed density is broadly compatible with a closed stationary diffusion.

### Open-process interpretation

If instead one finds substantial deviations such as

\[
p_{\mathrm{emp}}(x) \neq p_{\mathrm{closed}}(x),
\qquad
J(x)\neq 0,
\qquad
S(x)\neq 0,
\]

then the empirical density behaves as if additional mechanisms are present, such as influx or efflux of probability mass, nonstationarity, misspecified drift, misspecified diffusion, or other violations of the closed one-dimensional approximation.

## Limitations

Several limitations should be kept in mind.

- The diagnostics depend directly on upstream drift and diffusion estimates.
- The method assumes a one-dimensional diffusion in \(x=\ln f\).
- The empirical density depends on histogram construction, smoothing, and observability filtering.
- Current and source estimates involve numerical derivatives and are noise-sensitive.
- This script produces point estimates only and does not propagate bootstrap uncertainty.

## Practical summary

In practical terms, `8-fp_diagnostics_build.py` performs the following sequence:

1. select one drift model and one dynamical slice,
2. reconstruct \(D(x)\) from a residual-based variance proxy,
3. build the empirical state-density \(p_{\mathrm{emp}}(x)\),
4. compute the closed stationary prediction \(p_{\mathrm{closed}}(x)\),
5. derive the empirical current \(J(x)\),
6. derive the source term \(S(x)\),
7. summarize agreement between empirical and closed stationary densities by JS and KS.

This makes the script a central step for evaluating whether the inferred clonotype dynamics are compatible with a closed stationary Fokker–Planck description.

## Conclusion

This script synthesizes fitted drift, residual-based diffusion, and empirical state occupancy into a unified Fokker–Planck diagnostic framework. It provides the quantitative ingredients needed to compare closed stationary predictions with empirical repertoire organization and to identify abundance regions in which the observed system behaves as if it were open, nonstationary, or otherwise inconsistent with a closed one-dimensional diffusion process.
