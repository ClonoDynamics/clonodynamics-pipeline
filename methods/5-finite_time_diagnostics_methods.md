# Supplementary Methods — Model-agnostic local diagnostics and finite-time inference (`5-finite_time_diagnostics.py`)

## Overview

This script performs **model-agnostic characterization of clonotype dynamics** from a transition dataset. It implements two complementary analysis blocks:

1. **Local diagnostics**: estimation of conditional displacement statistics and variance structure as a function of initial abundance and time lag.
2. **Finite-time inference**: extraction of effective drift and diffusion parameters from transition statistics across multiple time intervals.

All computations are performed in **natural logarithm space (ln)**, regardless of the input representation.

The script operates on transition-level data generated from clonotype trajectories and provides the statistical foundation for subsequent dynamical modeling.

This documentation is based on the uploaded script `5-finite_time_diagnostics.py`. fileciteturn4file0

---

## Input Data Structure

The input is a transition table with at least the following columns:

- `subject`: subject identifier  
- `aaSeqCDR3`: clonotype identifier  
- `dt`: time lag between observations  
- `x0`: log-frequency at initial time  
- `x1`: log-frequency at final time  
- `dx`: displacement \( x_1 - x_0 \)  
- `obs_class`: detectability class  

Optional:

- `w`: transition weight  

---

## Log-space normalization

The script accepts input in either:

- natural log (`ln`)
- base-10 log (`log10`)

If input is in base-10:

\[
x^{(\ln)} = x^{(\log_{10})} \cdot \ln(10)
\]

All outputs are expressed in **natural log space**:

- \( \ln f_0 \)
- \( \ln f_1 \)
- \( \Delta x \)

Additionally, the geometric mean abundance is defined as:

\[
\ln f_{\mathrm{geo}} = \frac{1}{2}(\ln f_0 + \ln f_1)
\]

---

## Filtering and selection

Transitions are filtered according to:

### Time constraints

\[
dt > dt_{\min}, \quad dt \le dt_{\max}
\]

Optional exact selection:

\[
dt \in \{dt_1, dt_2, ...\}
\]

### Detectability

- `--only_TT`: restrict to TT transitions  
- or specify classes explicitly  

---

## Weighting

If a weight column is provided:

\[
w_i > 0
\]

weights are incorporated in all estimators.

If `weight_col = "1"`:

\[
w_i = 1
\]

---

## Binning strategy

Two variables are used:

- \( \ln f_0 \) (initial abundance)
- \( \ln f_{\mathrm{geo}} \) (mean abundance)

Bins are defined using:

### Quantile binning

\[
\text{edges} = \mathrm{quantile}(x, q)
\]

ensuring approximately equal sample sizes per bin.

### Fixed binning

User-supplied bin edges can be used for consistency across datasets.

---

## Local diagnostics

### Drift estimation

For each bin:

\[
b(x_0, dt) \approx \mathrm{median}(\Delta x \mid x_0, dt)
\]

Additionally:

\[
P(\Delta x > 0 \mid x_0, dt)
\]

is computed.

Bootstrap resampling provides confidence intervals.

---

### Mean squared displacement (MSD)

Global:

\[
\mathrm{MSD}(dt) = \mathbb{E}[\Delta x^2 \mid dt]
\]

Binned:

\[
\mathrm{MSD}(x_0, dt), \quad \mathrm{MSD}(f_{\mathrm{geo}}, dt)
\]

---

### Variance estimation via MAD

Robust variance is estimated using:

\[
\mathrm{MAD} = \mathrm{median}(|\Delta x - \mathrm{median}(\Delta x)|)
\]

Converted to standard deviation:

\[
\sigma \approx 1.4826 \cdot \mathrm{MAD}
\]

Variance:

\[
\mathrm{Var} \approx \sigma^2
\]

---

## Finite-time inference

### Zero-crossing point

The drift curve satisfies:

\[
b(x^*) = 0
\]

The zero-crossing \( x^* \) is identified by sign change.

Interpolation:

\[
x^* = x_i - \frac{b(x_i)}{b(x_{i+1}) - b(x_i)} (x_{i+1} - x_i)
\]

---

### Local slope estimation

A local linear fit is performed:

\[
b(x) \approx a + s(x - x^*)
\]

Slope:

\[
s = \frac{db}{dx}
\]

---

### Characteristic timescale

\[
\tau = \frac{dt}{-s}
\]

Valid only when:

\[
s < 0
\]

---

### Diffusion coefficient

\[
D = \frac{\mathrm{Var}(\Delta x)}{2dt}
\]

---

## Bootstrap procedure

Cluster bootstrap is used with clustering over:

\[
(\text{subject}, \text{clonotype})
\]

This preserves temporal dependencies.

Confidence intervals are computed as:

\[
[\alpha/2, 1 - \alpha/2]
\]

---

## Outputs

### Local diagnostics

- `drift_by_dt_xbins.csv`
- `msd_by_dt.csv`
- `msd_by_dt_xbins.csv`
- `msd_by_dt_freqbins.csv`
- `dt_freqbin_counts.csv`

### Finite-time

- `finite_time_summary.csv`
- `finite_time_by_dt_xbins_boot.csv`

---

## Statistical interpretation

The system is approximated as:

\[
x(t + dt) = x(t) + b(x)dt + \sqrt{2D(x)dt}\,\eta
\]

with:

\[
\eta \sim \mathcal{N}(0,1)
\]

Estimators:

\[
b(x) \approx \frac{\mathbb{E}[\Delta x]}{dt}
\]

\[
D(x) \approx \frac{\mathrm{Var}(\Delta x)}{2dt}
\]

---

## Bias and dependence

### Non-independence

Transitions share observations:

\[
(t_0 \to t_2), (t_0 \to t_3)
\]

### Heterogeneity

\[
P(dt), P(x), P(\text{class})
\]

are non-uniform.

### Implication

Bootstrap at cluster level is required.

---

## Role in pipeline

This script performs:

\[
\text{transitions} \rightarrow \text{moments} \rightarrow \text{dynamical inference}
\]

---

## Conclusion

This module provides a **model-independent statistical characterization** of clonotype dynamics and extracts effective drift and diffusion structure from empirical transition data.
