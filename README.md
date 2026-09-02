# ClonoDynamics

## Noise-aware inference of longitudinal clonotype dynamics from repertoire sequencing

ClonoDynamics is a modular computational framework for reconstructing and characterizing longitudinal T cell receptor (TCR) clonotype dynamics from repertoire sequencing data with technical replicates.

The current workflow is built around a strict separation between:

1. **measurement-noise and latent-state inference**;
2. **trajectory and transition construction**;
3. **pseudo-longitudinal technical benchmarking**;
4. **genuine longitudinal dynamical analysis**;
5. **longitudinal-vs-pseudo comparison**;
6. **sensitivity analysis of temporal structure, observation-domain inclusion, and operational observation thresholds**.

A central design principle is that **no longitudinal dynamical model is imposed during latent-frequency inference**. Paired technical replicates are first used to estimate an uncertainty-aware latent clonotype state at each subject–timepoint. Longitudinal trajectories and finite-time transitions are constructed only afterwards and are analyzed empirically.

ClonoDynamics is intended for RepSeq studies in which technical replicate structure is available and sampling noise, stochastic dropout, conditioning geometry, observation-domain definitions, and operational observation thresholds must be distinguished from genuine temporal structure.

For routine use, the complete software is exposed through a single top-level launcher:

```bash
python3 code/dynamics/clonodynamics.py
```

The launcher separates **Analysis** from **Figures**. Analysis is delegated to four branch-specific orchestrators, whereas plotting is delegated to a dedicated plot-only orchestrator; the launcher contains no scientific analysis itself.

---

## Scope

The current pipeline is designed to:

- characterize repertoire abundance distributions and heavy-tail structure;
- infer latent clonotype frequencies from paired technical replicates;
- quantify posterior uncertainty, detectability, and replicate consistency;
- reconstruct longitudinal clonotype trajectories;
- build abundance-resolved finite-time transition datasets;
- benchmark conditioning geometry and same-replicate coupling using pseudo-longitudinal controls;
- estimate replicate-decoupled genuine longitudinal forward drift;
- quantify abundance-dependent fluctuation magnitude across finite temporal lags;
- compare genuine longitudinal dynamics with pseudo-longitudinal technical baselines;
- test temporal scaling of clonotype fluctuations;
- assess calendar-position and interval-composition sensitivity;
- evaluate observation-domain sensitivity at a fixed operational threshold;
- test robustness to alternative operational observation thresholds.

The current primary workflow does **not** require a predefined stochastic differential equation, Fokker–Planck stationary reconstruction, or nonlinear drift model.

---

# Pipeline architecture

The repository contains 16 numbered public analysis steps.

Steps 1–6 define the shared repertoire-to-transition layer. Step 6 is descriptive and is not itself a computational prerequisite for the downstream branches.

```text
RAW REPERTOIRE TABLES
        |
        +-------------------------> STEP 1
        |                           repertoire characterization
        |
        v
     STEP 2
noise-aware latent inference
        |
        +-------------------------> STEP 3
        |                           posterior-quality characterization
        |
        v
     STEP 4
latent trajectories
        |
        v
     STEP 5
latent transitions
        |
        +-------------------------> STEP 6
        |                           transition-dataset characterization
        |
        +--------------------------------------------------------------------------+
        |                                                                          |
        |                                                                          |
        | PSEUDO-LONGITUDINAL                                                      | GENUINE LONGITUDINAL
        |                                                                          |
        +--> STEP 7                                                                 +--> STEP 9
        |    conditioning / ordering benchmark                                     |    replicate-decoupled
        |                                                                          |    forward drift
        +--> STEP 8                                                                 |
             xstar representation assessment                                       +--> STEP 11
                                                                                   |    fluctuation dynamics
                                                                                   |        |
                                                                                   |        v
                                                                                   |    STEP 13
                                                                                   |    temporal scaling
                                                                                   |        |
                                                                                   |        v
                                                                                   |    STEP 14
                                                                                   |    interval-position structure
                                                                                   |
                                                                                   +--> STEP 15
                                                                                   |    fixed-threshold
                                                                                   |    observation-domain sensitivity
                                                                                   |
                                                                                   +--> STEP 16
                                                                                        operational observation-
                                                                                        threshold robustness


CROSS-DATASET COMPARISONS

STEP 9  + STEP 7  -----------------> STEP 10
                                      longitudinal vs pseudo
                                      forward drift

STEP 11 + STEP 8 + pseudo STEP 4 --> STEP 12
                                      longitudinal vs pseudo
                                      fluctuation baseline
```

The numbering follows the manuscript analysis sequence, but the computational graph is **not a simple linear chain**.

Important consequences:

- Step 1 and Step 3 are characterization/diagnostic layers.
- Step 4 reads Step-2 latent estimates directly.
- Step 5 uses Step-4 trajectories and can also use Step-2 full posterior files.
- Step 8 does not consume Step-7 output files.
- Step 10 is a cross-dataset comparison and does not feed Step 11.
- Step 12 is a cross-dataset comparison and does not feed Step 13.
- Step 13 fits genuine longitudinal Step-11 fluctuation data; the Step-12 pseudo baseline is **not subtracted** beforehand.
- Step 14 uses Step-5 transitions together with Step-11 and Step-13 outputs.
- Step 15 keeps the reference operational threshold fixed at `alpha = 0.05` and changes the included transition domain using the current Step-9/11/13 implementations.
- Step 16 varies the operational threshold across `0.10`, `0.05`, `0.025`, and `0.01`, reclassifies endpoint T/F states, and reuses the canonical Step-9/11/13 analyses.
- Steps 15 and 16 are parallel sensitivity branches: Step 16 does not consume Step-15 outputs.

---

# Analysis steps

| Step | Script | Branch | Main role |
|---|---|---|---|
| 1 | `1-repertoire_characterization.py` | shared | replicate-level repertoire characterization and heavy-tail fitting |
| 2 | `2-noise_aware_latent_inference.py` | shared | noise-aware latent-frequency inference from paired technical replicates |
| — | `noiseK_latent.py` | shared dependency | statistical engine used by Step 2 |
| 3 | `3-latent_posterior_quality_characterization.py` | shared diagnostic | posterior precision, detectability, entropy, and replicate consistency |
| 4 | `4-latent_trajectory_construction.py` | shared | standardized clonotype trajectories |
| 5 | `5-latent_transition_construction.py` | shared | finite-time clonotype transitions |
| 6 | `6-latent_transition_dataset_characterization.py` | shared diagnostic | structural characterization of the transition dataset |
| 7 | `7-pseudo_longitudinal_conditioning_benchmark.py` | pseudo | conditioning geometry, replicate decoupling, pseudo-time ordering null |
| 8 | `8-pseudo_longitudinal_representation_assessment.py` | pseudo | xstar-conditioned displacement representation assessment |
| 9 | `9-replicate_decoupled_forward_drift.py` | longitudinal | genuine replicate-decoupled forward drift |
| 10 | `10-longitudinal_vs_pseudo_forward_drift.py` | comparison | longitudinal forward drift vs pseudo ordering null |
| 11 | `11-longitudinal_fluctuation_dynamics.py` | longitudinal | abundance-resolved finite-time fluctuation statistics |
| 12 | `12-longitudinal_vs_pseudo_fluctuation_baseline.py` | comparison | longitudinal fluctuation magnitude vs pseudo technical baseline |
| 13 | `13-temporal_fluctuation_scaling.py` | longitudinal | finite-lag temporal-scaling model comparison |
| 14 | `14-interval_position_structure.py` | longitudinal sensitivity | calendar-position and anchored-lag sensitivity |
| 15 | `15-detectability_boundary_sensitivity.py` | longitudinal sensitivity | fixed-threshold TT vs boundary-crossing observation-domain sensitivity |
| 16 | `16-analyze_observation_threshold_robustness.py` | longitudinal sensitivity | robustness to alternative operational observation thresholds |


Step-specific figures are generated by separate plot-only scripts under `code/plotting/`. The Step-16 plotter is `16-plot_observation_threshold_robustness.py`.

---

# Step 2 statistical model

For clonotype latent frequency \(f\), technical replicate counts are modeled as:

\[
C_r \mid f,N_r,\kappa
\sim
\mathrm{NB}(\mu=fN_r,\mathrm{size}=\kappa),
\]

with:

\[
\mathrm{Var}(C_r\mid f)
=
\mu+\frac{\mu^2}{\kappa}.
\]

Latent frequency is represented on a log-spaced grid with a **power-law-shaped discrete prior mass**. Pair-specific model parameters are estimated independently for each technical-replicate pair.

Posterior inference provides:

- latent-frequency mean, median, mode, and credible intervals;
- latent log-frequency uncertainty;
- replicate-specific posterior summaries;
- replicate-posterior overlap diagnostics;
- posterior-predictive detectability;
- full posterior matrices for downstream uncertainty propagation.

Two quantities should not be conflated:

- `p_detect_state` is a posterior-predictive detectability probability;
- `observable` is an empirical Step-2 classification based on fitted-model `logP` and the selected observability `alpha`.

The primary workflow uses `alpha = 0.05`. Step 16 revisits only the operational T/F classification at alternative `alpha` values; it does not alter the latent-frequency model or refit the latent posterior.

No longitudinal dynamics are imposed at this stage.

---

# Input data

ClonoDynamics starts from **clonotype-level repertoire tables**, not raw FASTQ files.

For Steps 1–2, each repertoire table must contain at least:

| Column | Meaning |
|---|---|
| `aaSeqCDR3` | amino-acid CDR3 clonotype identifier |
| `readCount` | clonotype sequencing read count |

Technical replicates are identified from filenames.

Default naming convention:

```text
subject_time-replicate.ext
```

Example:

```text
1_1-1.tsv
1_1-2.tsv
1_2-1.tsv
1_2-2.tsv
```

The shared input readers support the text and columnar formats used by the pipeline, including TSV/CSV and Parquet.

---

# Repository structure

Recommended public layout:

```text
ClonoDynamics/
│
├── code/
│   ├── dynamics/
│   │   │
│   │   ├── 1-repertoire_characterization.py
│   │   ├── 2-noise_aware_latent_inference.py
│   │   ├── noiseK_latent.py
│   │   ├── 3-latent_posterior_quality_characterization.py
│   │   ├── 4-latent_trajectory_construction.py
│   │   ├── 5-latent_transition_construction.py
│   │   ├── 6-latent_transition_dataset_characterization.py
│   │   │
│   │   ├── 7-pseudo_longitudinal_conditioning_benchmark.py
│   │   ├── 8-pseudo_longitudinal_representation_assessment.py
│   │   │
│   │   ├── 9-replicate_decoupled_forward_drift.py
│   │   ├── 10-longitudinal_vs_pseudo_forward_drift.py
│   │   ├── 11-longitudinal_fluctuation_dynamics.py
│   │   ├── 12-longitudinal_vs_pseudo_fluctuation_baseline.py
│   │   ├── 13-temporal_fluctuation_scaling.py
│   │   ├── 14-interval_position_structure.py
│   │   ├── 15-detectability_boundary_sensitivity.py
│   │   ├── 16-analyze_observation_threshold_robustness.py
│   │   │
│   │   ├── clonodynamics.py
│   │   ├── clonodynamics_orchestrator_steps_1_6.py
│   │   ├── clonodynamics_orchestrator_pseudo_reference.py
│   │   ├── clonodynamics_orchestrator_longitudinal_dynamics.py
│   │   ├── clonodynamics_orchestrator_longitudinal_vs_pseudo.py
│   │   └── clonodynamics_orchestrator_plotting.py
│   │
│   └── plotting/
│       ├── 1-plot_repertoire_characterization.py
│       ├── 3-plot_latent_posterior_quality_characterization.py
│       ├── 7-plot_pseudo_longitudinal_conditioning_benchmark.py
│       ├── 8-plot_pseudo_longitudinal_representation_assessment.py
│       ├── 9-plot_replicate_decoupled_forward_drift.py
│       ├── 10-plot_longitudinal_vs_pseudo_forward_drift.py
│       ├── 11-plot_longitudinal_fluctuation_dynamics.py
│       ├── 12-plot_longitudinal_vs_pseudo_fluctuation_baseline.py
│       ├── 13-plot_temporal_fluctuation_scaling.py
│       ├── 14-plot_interval_position_structure.py
│       ├── 15-plot_detectability_boundary_sensitivity.py
│       └── 16-plot_observation_threshold_robustness.py
│
├── methods/
│   └── methodological documentation
│
├── example_repseq/
│   └── lightweight example repertoire tables
│
├── Example_README.md
└── README.md
```

`clonodynamics.py` is the recommended top-level entry point. The five orchestrator filenames shown above—four analysis orchestrators plus the plotting orchestrator—are the canonical public names used by the launcher. Development suffixes such as `_v2` or `_fixed` should be removed before public release so that the GitHub interface has stable filenames.

---

# Execution and orchestration

## Master launcher

The recommended entry point for routine use is:

```bash
python3 code/dynamics/clonodynamics.py
```

The launcher contains **no scientific analysis and no Step-specific analytical logic**. It separates two responsibilities:

- **Analysis** dispatches to four dependency-aware scientific orchestrators;
- **Figures** dispatches to a dedicated plot-only orchestrator.

```text
                              clonodynamics.py
                                     |
                         +-----------+-----------+
                         |                       |
                         v                       v
                      ANALYSIS                 FIGURES
                         |                       |
                         |                       v
                         |       clonodynamics_orchestrator_plotting.py
                         |                       |
                         |              dedicated Step plotters
                         |
             +-----------+-----------+-----------+-----------+
             |                       |                       |
             v                       v                       v
       ORCHESTRATOR A          ORCHESTRATOR B          ORCHESTRATOR C
          Steps 1-6               Steps 7-8        Steps 9,11,13,14,15,16
             |                       |                       |
             +-----------------------+-----------+-----------+
                                                 |
                                                 v
                                          ORCHESTRATOR D
                                             Steps 10,12
```

Interactive execution first asks whether to enter **Analysis** or **Figures**. The Analysis menu exposes:

```text
[1] Repertoire -> latent transitions       Steps 1-6
[2] Build pseudo-longitudinal reference    Steps 7-8
[3] Analyze genuine longitudinal dynamics  Steps 9,11,13,14,15,16
[4] Compare longitudinal vs pseudo         Steps 10,12
```

The Figures menu exposes the available plotters for Steps `1`, `3`, and `7-16`. For each selected plotter, the plotting orchestrator asks for the corresponding analysis directory and writes figures to:

```text
<input-dir>/figures/
```

Plotters read already-generated analysis outputs. They do not rebuild transitions, change operational classifications, or refit the scientific analyses.

Direct analysis dispatch remains available:

```bash
python3 code/dynamics/clonodynamics.py --workflow shared
python3 code/dynamics/clonodynamics.py --workflow pseudo
python3 code/dynamics/clonodynamics.py --workflow longitudinal
python3 code/dynamics/clonodynamics.py --workflow compare
```

The equivalent explicit form is:

```bash
python3 code/dynamics/clonodynamics.py --mode analysis --workflow longitudinal
```

Direct figure dispatch is also supported. For example, to plot Step 16:

```bash
python3 code/dynamics/clonodynamics.py \
  --mode figures \
  --plot step16 \
  --input-dir /path/to/16-observation_threshold_robustness
```

Arguments following `--` are forwarded unchanged to the selected child orchestrator or plotter. The `--dry-dispatch` option resolves and prints the immediate child command without launching it.

## Advanced/direct orchestrator access

The analysis orchestrators and plotting orchestrator remain directly executable for development, automation, selective reruns, and reproducibility.

## A. Shared repertoire-to-transition workflow

```text
clonodynamics_orchestrator_steps_1_6.py
```

Handles:

```text
Steps 1–6
```

for either:

```text
longitudinal
pseudo-longitudinal
```

Supported execution modes:

- all Steps 1–6;
- one Step;
- any contiguous Step range.

When execution begins after Step 1, the orchestrator resolves the required upstream artifacts automatically.

```bash
python3 code/dynamics/clonodynamics_orchestrator_steps_1_6.py
```

---

## B. Pseudo-longitudinal reference workflow

```text
clonodynamics_orchestrator_pseudo_reference.py
```

Handles:

```text
Step 7
Step 8
Steps 7–8
```

and resolves the required pseudo Step-4/5 inputs automatically.

```bash
python3 code/dynamics/clonodynamics_orchestrator_pseudo_reference.py
```

---

## C. Genuine longitudinal dynamics workflow

```text
clonodynamics_orchestrator_longitudinal_dynamics.py
```

Handles:

```text
Step 9
Step 11
Step 13
Step 14
Step 15
Step 16
```

Supported workflows include:

- complete longitudinal dynamics;
- core dynamics;
- fluctuation chain;
- sensitivity analyses;
- one Step;
- custom Step subsets.

Dependencies produced earlier in the same selected workflow are passed directly to downstream Steps. Step 16 is available in the complete and sensitivity workflows and remains parallel to Step 15.

```bash
python3 code/dynamics/clonodynamics_orchestrator_longitudinal_dynamics.py
```

---

## D. Longitudinal-vs-pseudo comparison workflow

```text
clonodynamics_orchestrator_longitudinal_vs_pseudo.py
```

Handles:

```text
Step 10
Step 12
Steps 10 and 12
```

It requires two independent upstream result trees:

```text
longitudinal results root
pseudo-longitudinal results root
```

and writes comparison outputs into a third, separate result tree.

When metadata are available, the orchestrator checks cross-dataset provenance, including:

- dataset type;
- observability alpha;
- common upstream script hashes.

```bash
python3 code/dynamics/clonodynamics_orchestrator_longitudinal_vs_pseudo.py
```

---

## E. Plotting workflow

```text
clonodynamics_orchestrator_plotting.py
```

The plotting orchestrator manages plot-only scripts for Steps `1`, `3`, and `7-16`. It accepts an analysis-output directory, resolves the matching plotter under `code/plotting/`, and writes the resulting figures under `<input-dir>/figures/`.

```bash
python3 code/dynamics/clonodynamics_orchestrator_plotting.py
```

The plotting layer is deliberately separated from scientific analysis so that rerendering a figure cannot alter the analysis tables on which it is based.

---

# Recommended complete workflow

For routine execution, the complete analysis can be driven from the single master launcher:

```bash
python3 code/dynamics/clonodynamics.py
```

The four stages below correspond directly to its four Analysis choices. Figures can be generated afterwards from the same launcher without rerunning the analyses.

## 1. Process the genuine longitudinal dataset

From the master launcher select:

```text
[1] Repertoire -> latent transitions
    Steps 1-6
```

Then select:

```text
dataset type: longitudinal
```

This builds the shared latent-state, trajectory, and transition layers for the genuine longitudinal dataset.

---

## 2. Process the pseudo-longitudinal dataset and build the pseudo reference

First run the shared branch again from the master launcher:

```text
[1] Repertoire -> latent transitions
    Steps 1-6
```

and select:

```text
dataset type: pseudo-longitudinal
```

Then return to the master launcher and select:

```text
[2] Build pseudo-longitudinal reference
    Steps 7-8
```

This constructs the conditioning/ordering benchmark and the xstar-conditioned representation assessment used by the cross-dataset comparisons.

---

## 3. Run genuine longitudinal dynamics

From the master launcher select:

```text
[3] Analyze genuine longitudinal dynamics
    Steps 9,11,13,14,15,16
```

The longitudinal orchestrator can run the complete branch or a selected subset, while resolving required upstream Step-4/5 and intermediate Step-11/13 dependencies. The complete branch includes both P15 fixed-threshold domain sensitivity and P16 operational-threshold robustness.

---

## 4. Compare longitudinal and pseudo results

After the longitudinal and pseudo branches are available, select:

```text
[4] Compare longitudinal vs pseudo
    Steps 10,12
```

The comparison orchestrator reads the two independent upstream result trees and writes Step-10/12 outputs into a third, separate comparison-results tree.

---

## Direct orchestrator execution

The master launcher is a convenience entry point, not a replacement for the four scientific orchestrators or the plotting orchestrator. Each orchestrator can still be called directly:

```bash
python3 code/dynamics/clonodynamics_orchestrator_steps_1_6.py
python3 code/dynamics/clonodynamics_orchestrator_pseudo_reference.py
python3 code/dynamics/clonodynamics_orchestrator_longitudinal_dynamics.py
python3 code/dynamics/clonodynamics_orchestrator_longitudinal_vs_pseudo.py
python3 code/dynamics/clonodynamics_orchestrator_plotting.py
```

Direct execution is useful for development, automation, selective reruns, and reproducibility workflows.

---

# What each downstream branch tests

## Pseudo-longitudinal branch

### Step 7 — conditioning and pseudo-time benchmark

Evaluates:

- `x0`, `xmid`, and `xstar` conditioning geometry;
- same-replicate vs replicate-decoupled forward estimators;
- AB and BA cross-fit folds;
- pseudo-time ordering and exact reversal.

The pseudo ordering ensemble provides the technical null used by Step 10.

### Step 8 — representation assessment

Fixes `xstar` conditioning and compares displacement representations such as:

```text
latent_frequency
observed_frequency
log_count
```

The latent-frequency representation is used for downstream fluctuation analysis.

Step 8 does not automatically declare a statistical "winner"; it provides a structured representation assessment.

---

## Genuine longitudinal branch

### Step 9 — forward drift

Primary estimand:

\[
E[\Delta x_{\mathrm{latent}}\mid x_0]
\]

with replicate-decoupled AB/BA conditioning and symmetric `cross_combined` aggregation.

### Step 11 — fluctuation dynamics

Primary analysis:

```text
representation = latent
conditioning   = xstar
mode           = TT
```

Primary metric:

\[
\mathrm{Var}(\Delta x_{\mathrm{latent}}\mid x_\star).
\]

MSD is complementary; MAD is descriptive.

### Step 13 — temporal scaling

Tests finite-lag models across an internally comparable abundance core.

Candidate model classes include:

- constant;
- signed linear;
- non-negative incremental linear;
- zero-origin linear;
- power-law;
- saturating.

The analysis asks whether fluctuations show systematic temporal accumulation over the observed finite lags.

### Step 14 — interval-position structure

Tests whether lag-dependent patterns may reflect:

- calendar interval position;
- subject composition;
- common-start anchoring;
- common-end anchoring.

### Step 15 — fixed-threshold observation-domain sensitivity

Step 15 keeps the primary operational threshold fixed at `alpha = 0.05` and compares the primary TT estimands with boundary-crossing alternatives:

```text
Forward:
TT vs T0PLUS = TT + TF

Fluctuations:
TT vs NON_FF = TT + TF + FT

Temporal scaling:
TT vs NON_FF
```

Step 15 does not relabel endpoint states or redefine the operational threshold. It changes the included transition population using the T/F labels already assigned under the reference threshold.

### Step 16 — operational observation-threshold robustness

Step 16 varies the numerical operational threshold while leaving latent-frequency inference unchanged. The canonical threshold sweep is:

```text
alpha = 0.10, 0.05, 0.025, 0.01
reference alpha = 0.05
```

For each threshold, endpoint states are reclassified as T or F and transitions are reassigned to TT, TF, FT, or FF. The analyzer then reuses the canonical Step-9, Step-11, and Step-13 implementations to quantify:

- pooled transition-class composition and TT retention relative to `alpha = 0.05`;
- replicate-decoupled Step-9 forward-drift profiles and differences from the reference threshold;
- Step-11 fluctuation curves and differences from the reference threshold;
- Step-13 temporal-model comparisons within threshold-specific support;
- the common abundance core available across all analyzed thresholds.

Step 16 is a parallel sensitivity branch and does not consume Step-15 outputs. Its analysis script and plotter are:

```text
code/dynamics/16-analyze_observation_threshold_robustness.py
code/plotting/16-plot_observation_threshold_robustness.py
```

The analyzer can run the complete sweep with `--alphas 0.10 0.05 0.025 0.01`, `--reference-alpha 0.05`, and `--run-all`.

---

# Cross-dataset comparisons

## Step 10 — forward drift vs pseudo ordering null

Inputs:

```text
longitudinal Step 9
pseudo Step 7
```

The datasets retain their native absolute abundance bins.

Clonotype observations are not re-binned across datasets. Already-estimated curves are compared only on common absolute `x0` support.

Primary representation:

```text
cross_combined
```

---

## Step 12 — fluctuation magnitude vs pseudo technical baseline

Inputs:

```text
longitudinal Step 11
pseudo Step 8
pseudo Step 4 trajectories
```

Primary comparison:

\[
\mathrm{Var}(\Delta x_{\mathrm{latent}}\mid x_\star).
\]

Complementary comparison:

\[
E[\Delta x_{\mathrm{latent}}^2\mid x_\star].
\]

Pseudo-time configurations provide an empirical technical fluctuation baseline on shared absolute `xstar` support.

This baseline should not be interpreted as an exact decomposition of technical and biological variance.

The pseudo baseline is **not subtracted** before Step-13 temporal-scaling analysis.

---

# Principal intermediate outputs

Important shared outputs include:

```text
Step 1
    repertoire_characteristics.csv
    repertoire_tail_model_scan.csv

Step 2
    latent_params_by_pair.csv
    latent_pair_qc.csv
    per_clone_latent_subject.parquet
    per_clone_latent_logP/*.posterior.npz

Step 3
    clonotype_latent_posterior_qc.parquet

Step 4
    latent_trajectories_long.parquet

Step 5
    latent_transitions.parquet
```

Important downstream outputs include:

```text
Step 9
    forward_drift_by_bin.csv

Step 11
    10_binned_dynamics_long.parquet

Step 13
    05_model_comparison_by_bin.csv
    07_interpretation_aid.csv

Step 14
    10_analysis_summary.csv

Step 15
    20_step15_summary.csv

Step 16 — selected cross-alpha analysis/plot tables
    01_class_composition_pooled.csv
    02_tt_retention.csv
    03_p09_forward_plot_long.csv
    04_p09_forward_difference_vs_reference.csv
    05_p11_core_curves_long.csv
    06_p11_difference_vs_reference_long.csv
    07_p13_core_model_comparison.csv
    08_cross_alpha_core_coverage.csv

Step 10
    09_comparison_summary.csv

Step 12
    10_comparison_summary.csv
```

Each analysis step also writes configuration, diagnostic, or manifest files appropriate to that analysis. Plot-only outputs are kept under the selected analysis directory in a separate `figures/` subdirectory.

---

# Typical result-tree organization

```text
project/
│
├── dataset_longitudinal/
│
├── dataset_longitudinal_clonodynamics_results/
│   └── longitudinal/
│       ├── 1-repertoire_characterization/
│       ├── 2-noise_aware_latent_inference/
│       ├── 3-latent_posterior_quality_characterization/
│       ├── 4-latent_trajectory_construction/
│       ├── 5-latent_transition_construction/
│       ├── 6-latent_transition_dataset_characterization/
│       ├── 9-replicate_decoupled_forward_drift/
│       ├── 11-longitudinal_fluctuation_dynamics/
│       ├── 13-temporal_fluctuation_scaling/
│       ├── 14-interval_position_structure/
│       ├── 15-detectability_boundary_sensitivity/
│       ├── 16-observation_threshold_robustness/
│       └── logs/
│
├── dataset_pseudo/
│
├── dataset_pseudo_clonodynamics_results/
│   └── pseudo-longitudinal/
│       ├── 1-repertoire_characterization/
│       ├── 2-noise_aware_latent_inference/
│       ├── 3-latent_posterior_quality_characterization/
│       ├── 4-latent_trajectory_construction/
│       ├── 5-latent_transition_construction/
│       ├── 6-latent_transition_dataset_characterization/
│       ├── 7-pseudo_longitudinal_conditioning_benchmark/
│       ├── 8-pseudo_longitudinal_representation_assessment/
│       └── logs/
│
└── clonodynamics_longitudinal_vs_pseudo_results/
    ├── 10-longitudinal-vs-pseudo-forward-drift/
    ├── 12-longitudinal-vs-pseudo-fluctuation-baseline/
    └── logs/
```

A plotted Step may additionally contain a `figures/` subdirectory generated by the plotting orchestrator.

---

# Reproducibility

ClonoDynamics is designed around explicit intermediate states and dependency-aware reruns.

The master launcher is intentionally stateless with respect to scientific inference: it selects either the Analysis or Figures layer and delegates execution. Provenance and reproducibility are owned by the selected scientific orchestrator; plotters consume existing outputs without modifying them.

The orchestrators record, where applicable:

- selected workflow;
- input paths;
- dataset type;
- observability alpha;
- for Step 16, the analyzed alpha set and reference alpha;
- Python executable;
- analysis parameters;
- script hashes;
- timestamps;
- Step-specific logs.

Individual analysis Steps can still be executed directly for testing, sensitivity analysis, or methodological development.

For routine use, `clonodynamics.py` is the recommended entry point. The analysis orchestrators validate upstream dependencies, isolate output trees, record provenance, and run the scientific Steps; the plotting orchestrator renders their outputs separately.

---

# Statistical interpretation

ClonoDynamics distinguishes several layers of inference.

### Measurement model

Technical replicates are used to estimate uncertainty-aware latent clonotype frequencies.

This stage models measurement noise without imposing temporal dynamics.

### Forward drift

Replicate decoupling limits direct same-replicate coupling between the conditioning coordinate and measured displacement.

### Fluctuation magnitude

Variance and MSD are estimated conditionally on latent abundance using common absolute abundance geometry.

### Pseudo-longitudinal reference

Pseudo-time ordering provides an empirical reference for structure that can arise without genuine longitudinal ordering.

### Temporal scaling

Finite-lag temporal models are fit only to the genuine longitudinal data within an internally comparable abundance core.

Absence of systematic positive accumulation should not be interpreted as proof of complete temporal independence.

### Observation domain and operational threshold

TT-only analyses are estimand-specific. Step 15 holds `alpha = 0.05` fixed and changes which transition classes enter the estimand; including boundary-crossing transitions therefore changes the analyzed population rather than simply correcting TT estimates.

Step 16 addresses a different question. It varies `alpha`, recalculates endpoint T/F labels, and measures how the Step-9, Step-11, and Step-13 conclusions change as the operational observation boundary moves. Fixed-threshold domain sensitivity and threshold sensitivity are therefore reported as distinct, parallel analyses.

---

# Example data

A lightweight example dataset can be distributed under:

```text
example_repseq/
```

using the same technical-replicate naming convention:

```text
1_1-1.tsv
1_1-2.tsv
1_2-1.tsv
1_2-2.tsv
...
```

The example data are intended for software testing and workflow demonstration and are not a substitute for the full study dataset.

---

# Requirements

Recommended environment:

```text
Python >= 3.9
```

Core packages used across the current pipeline include:

```text
numpy
pandas
scipy
polars
pyarrow
matplotlib
```

Individual modules may require additional packages.

For public release, a pinned environment file or `requirements.txt` should accompany the repository.

---

# Data availability

ClonoDynamics operates on processed clonotype-level repertoire tables.

Study-specific raw sequencing accessions, processed datasets, and archived analysis outputs should be linked here using the final public accession numbers and repository DOIs.

```text
Raw sequencing:
    [NCBI SRA accession]

Processed repertoire and analysis data:
    [Zenodo / Figshare / institutional repository DOI]

Source repository:
    https://github.com/ClonoDynamics/tcr-dynamics-noise-pipeline

Software release:
    [GitHub release / Zenodo DOI]
```

---

# Citation

If you use ClonoDynamics, please cite the associated methods paper:

```text
ClonoDynamics: noise-aware inference of longitudinal clonotype dynamics
from repertoire sequencing

[full citation / DOI to be added]
```

---

# Contact

Camillo Palmieri

For questions about the software, reproducibility, or methodological implementation, please use the GitHub issue tracker or contact the corresponding author.
