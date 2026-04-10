# Supplementary Methods — Full-pipeline bootstrap of FP flux and source (`9-boot_full_pipeline_flux.py`)

## Overview

This script implements a **full uncertainty propagation framework** for Fokker–Planck diagnostics by combining:

- drift refitting  
- diffusion re-estimation  
- empirical density reconstruction  
- flux and source computation  

within a **cluster bootstrap scheme**.

Unlike previous steps (point estimates), this script quantifies **confidence intervals for J(x) and S(x)** by resampling the entire pipeline.

Based on script: fileciteturn10file0

---

## Conceptual pipeline

Each bootstrap replicate performs:

\[
	ext{transitions} 
ightarrow b(x) 
ightarrow D(x) 
ightarrow p(x) 
ightarrow J(x), S(x)
\]

This ensures that **all sources of uncertainty are propagated jointly**.

---

## Bootstrap scheme

### Cluster resampling

Data are resampled at the clonotype level:

\[
(	ext{subject}, 	ext{aaSeqCDR3})
\]

This preserves temporal structure and avoids pseudoreplication.

---

## Step-by-step algorithm

For each bootstrap replicate \(r\):

### 1. Resample transitions

\[
\{(x_0, \Delta x)\}_r
\]

### 2. Refit drift

Estimate:

\[
b_r(x)
\]

using:
- OU_linear  
- hinge_plateau_smooth  

---

### 3. Estimate diffusion

Two modes:

#### Default (recommended)

\[
r_i = \Delta x_i - b_r(x_{0,i})
\]

\[
D_r(x) = rac{\mathrm{Var}(r \mid x)}{2dt}
\]

#### Optional

Use precomputed residuals.

---

### 4. Resample trajectories → empirical density

\[
p_r(x)
\]

Supports:

- ignore  
- filter (observable only)  
- weight (weighted density)

---

### 5. Compute flux

\[
J_r(x) = b_r(x)p_r(x) - rac{d}{dx}[D_r(x)p_r(x)]
\]

---

### 6. Compute source

\[
S_r(x) = rac{dJ_r(x)}{dx}
\]

---

## Diffusion estimation

Residual-based:

### Robust estimator

\[
\mathrm{Var} pprox (1.4826 \cdot \mathrm{MAD})^2
\]

### Alternative

\[
\mathrm{Var} = \mathbb{E}[r^2]
\]

---

## Empirical density

Estimated via histogram on fixed grid:

\[
p(x) = rac{	ext{counts}}{N \cdot \Delta x}
\]

Optional Gaussian smoothing applied.

---

## Bulk region

To avoid edge artifacts, metrics are computed on:

- trimmed grid (remove extremes)  
- density threshold filter  

---

## Global metrics

Computed on each bootstrap:

- \( \int |J(x)| dx \)
- \( \int J(x)^2 dx \)
- \( \int J(x) dx \)
- \( \int |S(x)| dx \)
- \( \int S(x)^2 dx \)

---

## Confidence intervals

For each quantity:

\[
	ext{mean},\quad 	ext{CI}_{lpha}
\]

computed across bootstrap replicates.

---

## Outputs

### Curve summaries

- J_boot_summary.csv  
- S_boot_summary.csv  

### Global metrics

- global_boot_summary.csv  

### Drift / diffusion summaries

- drift_boot_summary.csv  
- diffusion_boot_summary.csv  

### Raw bootstrap arrays

- J_boot_reps.npy  
- S_boot_reps.npy  

### Figures

- flux_J_with_CI.png  
- source_S_with_CI.png  

### Report

- report.txt  

---

## Observable handling

Three modes:

- ignore → use all data  
- filter → keep observable only  
- weight → weighted density  

---

## Interpretation

This step answers:

> Are flux and source significantly different from zero?

Key diagnostic:

- If CI of \( \int J(x) dx \) includes 0 → consistent with closed system  
- Otherwise → evidence for open dynamics  

---

## Advantages

- Full uncertainty propagation  
- Model comparison under noise  
- Robust to sampling variability  

---

## Limitations

- Computationally intensive  
- Depends on upstream model choices  
- Bootstrap assumes cluster independence  

---

## Role in workflow

Final validation step:

\[
	ext{Drift + Diffusion} 
ightarrow 	ext{FP diagnostics} 
ightarrow 	ext{Bootstrap validation}
\]

---

## Conclusion

This script provides a **rigorous statistical validation layer** for Fokker–Planck diagnostics, allowing quantitative assessment of whether observed repertoire dynamics are compatible with a closed stochastic system or require additional open-process mechanisms.
