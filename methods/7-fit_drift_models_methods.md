# Supplementary Methods — Drift model fitting and diffusion-from-residuals (`7-fit_drift_models_analysis_and_residuals.py`)

## Overview

This script performs **parametric fitting of drift functions** on transition-level data at a fixed lag \(dt\), and derives **diffusion estimates from model residuals**. It is **analysis-only** (no plots), producing tabular outputs for downstream diagnostics and figures.

Core steps:
1. Filter transitions (by \(dt\), class, weights)
2. Fit multiple drift models \( \hat b(x_0) \) to displacements \( \Delta x \)
3. Compare models via AIC/BIC
4. Compute per-transition residuals and squared residuals
5. Estimate **conditional diffusion** \(D(x_0)\) from residuals (binned)
6. (Optional) Cluster bootstrap for parameter and model-comparison uncertainty

Based on script: fileciteturn6file0

---

## Input

A transitions CSV with at least:
- `dt`, `x0`, `dx`
Optional:
- `obs_class` (for filtering)
- weight column (e.g. `w_clone`)
- identifiers (`subject`, `aaSeqCDR3`, `time0`, `time1`)

Filtering:
\[
dt = \text{--dt}, \qquad \text{obs\_class} = \text{--filter\_class (e.g. TT)}
\]

Finite rows only are retained.

---

## Models

The script fits a **set of drift models** \( \hat b(x) \approx \mathbb{E}[\Delta x \mid x_0=x] \):

### 1. OU-linear
\[
\hat b(x) = a + b x
\]
Parameters: \(a,b\).  
Equilibrium: \(x^* = -a/b\).  
Timescale: \(\tau = -1/b\) (for \(b<0\)).

### 2. Tanh-saturating
\[
\hat b(x) = A \tanh\!\left(\frac{x^* - x}{\Delta}\right)
\]

### 3. Rational-saturating
\[
\hat b(x) = A \frac{z}{\sqrt{1+z^2}}, \quad z=\frac{x^*-x}{\Delta}
\]

### 4. Hinge–plateau (smooth)
Smooth interpolation between linear regime and plateau:
\[
\hat b(x) = (a+bx)\,s(x) + (-A)\,[1-s(x)], \quad s(x)=1-\sigma\!\left(\frac{x-x_c}{\Delta}\right)
\]

---

## Estimation

### Loss function

All models are fit by minimizing weighted SSE:
\[
\mathrm{SSE} = \sum_i w_i (\Delta x_i - \hat b(x_{0,i}))^2
\]

- Unweighted: \(w_i=1\)
- Weighted: sanitized and clipped weights (quantile cap `--weight_cap_q`), then rescaled to mean 1.

Robust optimization uses `least_squares(..., loss="soft_l1")` for nonlinear models.

### Standardization (OU-linear)

For numerical stability:
\[
z = \frac{x-\mu}{\sigma}
\]
Fit in \(z\)-space, then map back to \(x\)-space.

---

## Model comparison

Information criteria (Gaussian residual approximation):
\[
\mathrm{AIC} = n\ln(\mathrm{SSE}/n) + 2k,\qquad
\mathrm{BIC} = n\ln(\mathrm{SSE}/n) + k\ln n
\]

Outputs:
- `fit_summary.csv` (per-model parameters, SSE, AIC/BIC)
- `model_compare.csv` (ΔAIC/ΔBIC vs baseline)

---

## Fitted values and residuals

For each model:
\[
\hat{\Delta x}_i = \hat b(x_{0,i}),\qquad
r_i = \Delta x_i - \hat{\Delta x}_i,\qquad
r_i^2
\]

Saved in:
- `fitted_values.csv`

---

## Diffusion from residuals

Assuming:
\[
\Delta x = b(x_0)\,dt + \epsilon, \quad \mathbb{E}[\epsilon]=0,\ \mathrm{Var}(\epsilon)=2D(x_0)dt
\]

Estimate diffusion via residuals:

### Mean-based
\[
\mathbb{E}[r^2 \mid x_0] \approx 2D(x_0)dt
\]

### Robust (MAD)
\[
\mathrm{MAD} = \mathrm{median}(|r - \mathrm{median}(r)|),\quad
\sigma \approx 1.4826\,\mathrm{MAD},\quad
\mathrm{Var}_{\mathrm{MAD}} = \sigma^2
\]

---

## Binning

Quantile bins over \(x_0\) with constraints:
- `nbins_diff`
- `min_per_bin_diff`

Per-bin outputs:
- median residual
- MAD
- robust variance (`var_mad`)
- mean \(r^2\)

Saved in:
- `diffusion_binned.csv`

Long-format residuals:
- `diffusion_residuals_long.csv`

---

## Drift grid

A regular grid \(x\in[q_{lo},q_{hi}]\) is built:
\[
x \in [Q_{q_{lo}}(x_0),\, Q_{q_{hi}}(x_0)]
\]

Predictions \( \hat b(x) \) for each model are saved in:
- `drift_function_grid.csv`

---

## Bootstrap (optional)

Cluster bootstrap over:
\[
(\text{subject}, \text{aaSeqCDR3})
\]

Outputs:
- `bootstrap_reps.csv`
- `bootstrap_summary.csv`

Statistics:
- parameter medians and CI
- ΔAIC/ΔBIC distributions

---

## Outputs (summary)

- `fit_summary.csv`
- `model_compare.csv`
- `fitted_values.csv`
- `drift_function_grid.csv`
- `binned_medians.csv`
- `diffusion_residuals_long.csv`
- `diffusion_binned.csv`
- (optional) bootstrap files

---

## Interpretation

- **Drift models** capture deterministic trends (mean reversion, saturation)
- **Residuals** capture stochastic fluctuations
- **Diffusion estimates** reveal heteroscedasticity \(D(x)\)

Key diagnostic:
- If OU is insufficient → nonlinear models reduce SSE and AIC
- Residual variance decreasing with \(x\) → abundance-dependent noise

---

## Limitations

- Fixed \(dt\) analysis
- Gaussian residual assumption in IC
- Residual-based diffusion assumes correct drift model
- Weighting may affect variance estimates

---

## Role in pipeline

\[
\text{transitions} \rightarrow \text{drift fit} \rightarrow \text{residuals} \rightarrow \text{diffusion}
\]

---

## Conclusion

This script provides a **joint estimation of deterministic (drift) and stochastic (diffusion) components** of clonotype dynamics, enabling comparison of competing models and quantification of abundance-dependent variability.
