# Supplementary Methods — OU finite-time validation (`6-finite_time_fit.py`)

## Overview

This script performs **finite-time validation of an Ornstein–Uhlenbeck (OU) model** using empirical estimates obtained from transition data.

It does not estimate drift/diffusion from raw data; instead, it:

1. Takes as input summary statistics
2. Calibrates an OU process at a reference timescale
3. Predicts variance across time
4. Compares predictions to empirical data

Based on script: fileciteturn5file0

---

## OU model

\[
dx = -\lambda x \, dt + \sqrt{2D} \, dW_t
\]

---

## Calibration

\[
\lambda = -\text{slope}(dt_{\mathrm{ref}})
\]

\[
D = \lambda \cdot \frac{\mathrm{Var}(dt_{\mathrm{ref}})}{1 - e^{-2\lambda dt_{\mathrm{ref}}}}
\]

---

## Predictions

\[
\mathrm{Var}(dt) = \frac{D}{\lambda}(1 - e^{-2\lambda dt})
\]

\[
D_{\mathrm{app}} = \frac{\mathrm{Var}(dt)}{2dt}
\]

---

## Outputs

- ou_pred_vs_empirical.csv  
- ou_fit_meta.json  

---

## Interpretation

Agreement with OU:

- variance saturates  
- apparent diffusion decreases  

---

## Role

Validation step linking empirical dynamics to OU theory.
