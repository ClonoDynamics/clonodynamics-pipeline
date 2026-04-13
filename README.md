# TCR Dynamics Noise-Aware Pipeline

## Overview

This repository implements a complete, noise-aware computational pipeline to quantify in vivo dynamics of human T cell receptor (TCR) clonotypes from longitudinal repertoire sequencing data. The framework integrates experimental noise modelling, trajectory reconstruction, and stochastic inference to link clonotype-level fluctuations to the global statistical organization of the repertoire.

The pipeline is designed to:

* infer experimental noise using technical replicates
* reconstruct clonotype trajectories across time
* quantify finite-time dynamics (drift and diffusion)
* fit nonlinear models of clonal regulation
* reconstruct stationary clone-size distributions via Fokker–Planck formalism
* propagate uncertainty through bootstrap analysis

---

## Repository Structure

```
code/
 ├── noise/
 │    ├── 1-powerlaw_qc_analysis.py
 │    ├── 2-noise_pipeline.py
 │    └── noiseK_nb.py
 │
 ├── dynamics/
 │    ├── 3-trajectory_builder.py
 │    ├── 4-build_transition_patterns.py
 │    ├── 5-finite_time_diagnostics.py
 │    ├── 6-finite_time_fit.py
 │    ├── 7-fit_drift_models_analysis_and_residuals.py
 │    ├── 8-fp_diagnostics_build.py
 │    └── 9-boot_full_pipeline_flux.py
 │
 ├── plotting/
 │    ├── 1-plot_powerlaw_qc.py
 │    ├── 2-plot_diagnostics_qc.py
 │    ├── 4-plot_transitions_alpha_summary.py
 │    ├── 5-plot_finite_time_diagnostics.py
 │    ├── 6-plot_finite_time_fit.py
 │    ├── 7-plot_fit_drift_models.py
 │    ├── 8-plot_fp_diagnostics.py
 │    └── 9-plot_fp_boot_full_pipeline.py
 │
 └── utils/
      ├── build_multipanel.py
      ├── extract_dynamical_summary.py
      ├── fp_sensitivity_diagnostics.py
      ├── plot_compare_finite_time_conditions.py
      └── spikein_with_background_v2.py

methods/
 ├── step-by-step methodological descriptions (.md)
 └── supplementary_methods_raw.pdf

results/
 ├── finite-time diagnostics and fits
 ├── drift inference outputs
 ├── stationary distributions (Fokker–Planck)
 ├── bootstrap analyses
 ├── alpha sensitivity diagnostics
 └── dynamical summaries

supplementary_data/
 ├── Supplementary_Data_S1.csv
 ├── Supplementary_Data_S2.csv
 └── Supplementary_Data_S3.csv
```

---

## Input Data

The pipeline operates on clonotype-level TCR repertoire data with the following fields:

* `aaSeqCDR3` — clonotype identifier
* `readCount` — number of sequencing reads
* `readFraction` — relative frequency

Samples should include **technical replicates** and follow a consistent naming convention:

```
subject_time-replicate
```

---

## Workflow

### 1. Noise modelling and QC

```bash
python code/noise/1-powerlaw_qc_analysis.py
python code/noise/2-noise_pipeline.py
```

### 2. Trajectory reconstruction

```bash
python code/dynamics/3-trajectory_builder.py
```

### 3. Transition definition

```bash
python code/dynamics/4-build_transition_patterns.py
```

### 4. Finite-time dynamical characterization

```bash
python code/dynamics/5-finite_time_diagnostics.py
python code/dynamics/6-finite_time_fit.py
```

### 5. Drift inference (nonlinear models)

```bash
python code/dynamics/7-fit_drift_models_analysis_and_residuals.py
```

### 6. Fokker–Planck reconstruction

```bash
python code/dynamics/8-fp_diagnostics_build.py
```

### 7. Full-pipeline bootstrap

```bash
python code/dynamics/9-boot_full_pipeline_flux.py
```

---

## Outputs

All results are saved in the `results/` directory as CSV files, including:

* drift functions (`drift_function_grid.csv`)
* diffusion estimates (`diffusion_binned.csv`)
* stationary distributions (`stationary_closed.csv`)
* bootstrap summaries and flux diagnostics

Supplementary datasets used in the manuscript are provided in:

```
supplementary_data/
```

---

## Reproducibility

The pipeline is fully modular and reproducible. Each step produces intermediate outputs that can be inspected independently. All scripts are designed to be run sequentially, but individual stages can be executed separately for validation or extension.

---

## Biological Interpretation

This framework enables quantitative inference of clonal dynamics under homeostatic constraints, revealing:

* mean-reverting behaviour of T cell clones
* nonlinear regulation of clonal expansion and contraction
* abundance-dependent fluctuations
* near-stationary organization of the repertoire

---

## Requirements

* Python ≥ 3.9
* numpy
* pandas
* scipy
* matplotlib
* scikit-learn

---

## Data Availability

Processed data and all analysis outputs required to reproduce the results are included in this repository. Raw sequencing data are not included due to ethical restrictions but are available upon reasonable request.

---

## Citation

If you use this pipeline, please cite:

> [Manuscript reference]

---

## Contact

Camillo Palmieri
For questions or collaborations, please contact the corresponding author.
