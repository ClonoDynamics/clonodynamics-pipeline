# ClonoDynamics

## Replicate-resolved inference and validation of longitudinal T-cell receptor repertoire dynamics

**ClonoDynamics** is a modular computational framework for reconstructing, quantifying, and validating longitudinal T-cell receptor (TCR) clonotype dynamics from repertoire sequencing data with paired technical replicates.

The current framework is built around a strict separation between:

1. **repertoire characterization and measurement-noise-aware state inference**;
2. **trajectory and transition construction**;
3. **pseudo-longitudinal technical-reference benchmarking**;
4. **genuine longitudinal forward, fluctuation, and temporal analyses**;
5. **support-conditioned synthetic positive controls**;
6. **calendar/composition and operational-observation sensitivity analyses**;
7. **read-only publication plotting from finalized analysis outputs**.

A central design principle is that **technical-replicate information is used without collapsing all downstream dynamics into a single latent trajectory**. Latent abundance remains an uncertainty-aware state representation and conditioning coordinate, whereas the primary longitudinal dynamical estimands preserve replicate-resolved observed measurements.

For routine use, the whole repository is exposed through a single interactive launcher:

```bash
python3 code/orchestrator_clonodynamics.py
```

The launcher contains the validated orchestration engines for analysis and figures. The scientific scripts themselves remain modular under `code/01_core` through `code/05_plotting`.

---

# Scientific scope

ClonoDynamics is designed for longitudinal RepSeq studies in which technical replicate structure is available and where measurement noise, stochastic dropout, conditioning geometry, observation-domain definitions, and operational observation thresholds must be distinguished from genuine temporal structure.

The current workflow can:

- characterize repertoire abundance distributions and heavy-tail structure;
- infer uncertainty-aware clonotype states from paired technical replicates;
- quantify posterior uncertainty and replicate-posterior agreement;
- construct multirepresentation longitudinal clonotype trajectories;
- assemble finite-time transition datasets across temporal lags;
- build a pseudo-longitudinal technical reference from exchangeable technical measurements;
- quantify replicate-decoupled observed forward dynamics;
- quantify replicate-consistent longitudinal fluctuations through cross-replicate covariance;
- compare genuine longitudinal dynamics with pseudo-longitudinal technical nulls;
- test finite-lag temporal scaling on a fixed abundance core and complete-case cohort;
- validate temporal-signal recovery with support-conditioned synthetic positive controls;
- test calendar-position and interval-composition sensitivity;
- test operational observation-domain sensitivity at a fixed threshold;
- test robustness across alternative operational observation thresholds;
- regenerate publication figures without rerunning scientific inference.

The framework does **not** require a predefined stochastic differential equation, Fokker–Planck stationary reconstruction, or a parametric nonlinear drift law.

---

# Documentation

The root `README.md` provides the conceptual overview, repository structure, and execution entry points for ClonoDynamics. Detailed scientific and computational methods are maintained under:

```text
docs/methods/
```

The methods documentation is organized into three complementary documents:

1. **[Primary analysis: Steps 1–13](docs/methods/01_primary_analysis.md)**  
   Repertoire characterization, paired-replicate state inference, multirepresentation trajectory and transition construction, pseudo-longitudinal technical calibration, observed replicate-decoupled forward dynamics, cross-replicate fluctuation dynamics, longitudinal-versus-pseudo comparisons, and temporal scaling.

2. **[Synthetic positive controls](docs/methods/02_positive_controls.md)**  
   Support-conditioned semi-synthetic validation, frozen empirical calibration, oracle dose construction, synthetic repertoire generation, recovery through the production pipeline, verification, and recovery summaries.

3. **[Robustness and sensitivity controls: Steps 14–16](docs/methods/03_robustness_controls.md)**  
   Calendar-position and interval-composition controls, fixed-threshold operational observation-domain sensitivity, and numerical observation-threshold robustness.

A compact index is available at **[`docs/methods/README.md`](docs/methods/README.md)**.

## Documentation and manuscript versioning

The files under `docs/methods/` are the **living scientific documentation** of the current ClonoDynamics software. They may evolve when repository paths, orchestration, implementation details, or explanatory text are improved.

The exact methods corresponding to a submitted or published manuscript should instead be frozen together with the matching software version using a tagged release, for example:

```text
v1.0.0-paper
```

This separation allows the public repository to remain maintainable while preserving an immutable software-and-methods record for the reported study.


---

# Repository architecture

The public repository is organized into five computational blocks plus a dedicated scientific-documentation layer.

```text
ClonoDynamics/
|
+-- README.md
|
+-- code/
|   +-- orchestrator_clonodynamics.py
|   |
|   +-- 01_core/
|   +-- 02_pseudo_reference/
|   +-- 03_validation/
|   +-- 04_controls/
|   `-- 05_plotting/
|
`-- docs/
    `-- methods/
        +-- README.md
        +-- 01_primary_analysis.md
        +-- 02_positive_controls.md
        `-- 03_robustness_controls.md
```

The master orchestrator is the only orchestration file required directly under `code/`. The `docs/methods/` directory provides the extended scientific documentation of the corresponding analysis blocks.

## 01_core

Genuine longitudinal repertoire processing and transition construction.

```text
01_core/
+-- 1-repertoire_characterization.py
+-- 2-multirepresentation_clonotype_state_inference.py
+-- noiseK_latent.py
+-- 3-latent_state_and_observation_model_diagnostics.py
+-- 4-multirepresentation_trajectory_assembly.py
+-- 5-longitudinal_transition_assembly.py
`-- 6-transition_support_characterization.py
```

## 02_pseudo_reference

Pseudo-longitudinal technical-reference construction.

```text
02_pseudo_reference/
+-- 00-generate_pseudo_design.py
+-- pseudo_1x12_build_upstream_multirepresentation_streaming.py
+-- 7-pseudo_forward_technical_null_compact.py
+-- pseudo_1x12_build_step8_transition_bank.py
`-- 8-pseudo_fluctuation_conditioning_validation_pairbank.py
```

## 03_validation

Primary genuine-longitudinal validation analyses plus synthetic positive controls.

```text
03_validation/
+-- 9-observed_replicate_decoupled_forward_drift.py
+-- 10-longitudinal_vs_pseudo_forward_null.py
+-- 11-cross_replicate_fluctuation_dynamics.py
+-- 12-longitudinal_vs_pseudo_cross_replicate_fluctuations.py
+-- 13-temporal_fluctuation_scaling.py
|
`-- positive_controls/
    +-- 01_freeze_empirical_calibration.py
    +-- 02_calibrate_oracle_doses.py
    +-- 03_generate_support_conditioned_repertoires.py
    +-- 04_summarize_temporal_recovery.py
    +-- 05_plot_synthetic_validation.py
    `-- verify_workflow.py
```

## 04_controls

Final robustness and sensitivity analyses.

```text
04_controls/
+-- 14-interval_position_structure.py
+-- 15-detectability_boundary_sensitivity.py
`-- 16-analyze_observation_threshold_robustness.py
```

## 05_plotting

Read-only publication plotting.

```text
05_plotting/
+-- 01_core/
|   +-- plot_step01_repertoire_characterization.py
|   +-- plot_step03_latent_state_diagnostics.py
|   `-- plot_step06_transition_support.py
|
+-- 02_pseudo_reference/
|   `-- plot_steps07_08_pseudo_reference.py
|
+-- 03_validation/
|   +-- plot_steps09_10_forward_dynamics.py
|   +-- plot_steps11_12_fluctuation_dynamics.py
|   +-- plot_step13_temporal_scaling.py
|   `-- positive_controls/
|       `-- plot_positive_controls.py
|
`-- 04_controls/
    +-- plot_steps14_16_robustness_summary.py
    `-- details/
        +-- plot_step14_interval_position_details.py
        +-- plot_step15_observation_domain_details.py
        `-- plot_step16_threshold_robustness_details.py
```

The three detailed Step-14/15/16 plotters are retained for inspection and supplementary rendering. The standard figure workflow uses the integrated Steps-14–16 robustness plotter.

---

# Computational graph

The analysis is modular and is **not** a simple 1→16 linear chain.

```text
                           GENUINE LONGITUDINAL DATA
                                      |
                                      v
                                01_core
                             Steps 1 → 6
                                      |
                      longitudinal_transitions.parquet
                                      |
                +---------------------+---------------------+
                |                                           |
                v                                           v
       03_validation Step 9                        03_validation Step 11
     observed AB/BA forward                    cross-replicate fluctuation
                |                                           |
                v                                  +--------+--------+
             Step 10                               v                 v
       real vs pseudo forward                   Step 12            Step 13
                ^                          real vs pseudo       temporal scaling
                |                             fluctuation           |
                |                                  ^                |
                |                                  |                v
                |                                  |              Step 14
                |                                  |       interval/composition controls
                |                                  |                |
                |                                  |                v
                |                                  |              Step 15
                |                                  |      fixed-threshold domain sensitivity
                |                                  |                |
                |                                  |                v
                |                                  |              Step 16
                |                                  |       cross-threshold robustness
                |                                  |
                |                                  |
                `--------------- 02_pseudo_reference ---------------+
                              Steps 7–8 technical reference
```

A separate support-conditioned synthetic validation branch reuses the production pipeline:

```text
final empirical core/validation results
              |
              v
freeze empirical calibration
              |
              v
calibrate oracle accumulation doses
              |
              v
generate R0p00 / R0p25 / R0p50 / R1p00
              |
              v
for each scenario:
Step 2 → Step 4 → Step 5 → Step 11 → Step 13
              |
              v
verify → summarize → plot
```

Important consequences:

- Step 1 and Step 3 are characterization/diagnostic layers.
- Step 4 consumes the Step-2 multirepresentation clonotype-state table.
- Step 5 consumes Step-4 trajectories and can propagate Step-2 full posterior uncertainty.
- Step 6 characterizes the transition universe but is not a computational prerequisite for Steps 9–13.
- Step 7 and Step 8 are pseudo-reference analyses and do not represent biological longitudinal time.
- Step 10 compares Step 9 with pseudo Step 7.
- Step 12 compares Step 11 with pseudo Step 8.
- Step 12 does **not** feed Step 13.
- Step 13 fits genuine longitudinal Step-11 results without subtracting the pseudo baseline.
- Step 14 inherits the Step-13 abundance core and complete-case cohort and analyzes interval-resolved Step-11 quantities.
- Step 15 keeps the operational threshold fixed at `alpha = 0.05` and changes the included operational observation domain.
- Step 16 varies the numerical operational threshold and reuses the validated Step-15 implementation while keeping latent-state inference fixed.
- Positive controls are a separate validation branch and are not pseudo-longitudinal randomizations.

---

# Analysis modules

- **Step 1 — core:** `1-repertoire_characterization.py` — repertoire-level descriptive characterization and heavy-tail analysis.
- **Step 2 — core:** `2-multirepresentation_clonotype_state_inference.py` — replicate-resolved latent/state inference and operational observability.
- **Core dependency:** `noiseK_latent.py` — count-noise / latent-state statistical engine used by Step 2.
- **Step 3 — core diagnostic:** `3-latent_state_and_observation_model_diagnostics.py` — posterior uncertainty, replicate agreement, and observation-model diagnostics.
- **Step 4 — core:** `4-multirepresentation_trajectory_assembly.py` — longitudinal assembly of latent and observed replicate-resolved states.
- **Step 5 — core:** `5-longitudinal_transition_assembly.py` — generic finite-time longitudinal transition table.
- **Step 6 — core diagnostic:** `6-transition_support_characterization.py` — coverage, observation classes, and estimand-support characterization.
- **Step 7 — pseudo:** `7-pseudo_forward_technical_null_compact.py` — pseudo technical null for replicate-decoupled forward structure.
- **Step 8 — pseudo:** `8-pseudo_fluctuation_conditioning_validation_pairbank.py` — pseudo technical reference for replicate-consistent fluctuation covariance.
- **Step 9 — validation:** `9-observed_replicate_decoupled_forward_drift.py` — genuine observed replicate-decoupled AB/BA forward dynamics.
- **Step 10 — validation:** `10-longitudinal_vs_pseudo_forward_null.py` — genuine forward dynamics versus pseudo technical null.
- **Step 11 — validation:** `11-cross_replicate_fluctuation_dynamics.py` — replicate-consistent longitudinal fluctuation covariance.
- **Step 12 — validation:** `12-longitudinal_vs_pseudo_cross_replicate_fluctuations.py` — genuine cross-replicate covariance versus pseudo technical null.
- **Step 13 — validation:** `13-temporal_fluctuation_scaling.py` — finite-lag temporal scaling on a fixed abundance core.
- **Step 14 — controls:** `14-interval_position_structure.py` — calendar-position, anchoring, and subject-composition controls.
- **Step 15 — controls:** `15-detectability_boundary_sensitivity.py` — fixed-threshold operational observation-domain sensitivity.
- **Step 16 — controls:** `16-analyze_observation_threshold_robustness.py` — cross-threshold robustness of operational observation-domain conclusions.

---

# Core statistical representation

## Count-noise model and latent state

At each subject–timepoint, paired technical replicates are modeled independently of any longitudinal dynamical law.

For latent clonotype frequency \(f\), replicate counts are modeled using an overdispersed count model of the form

$$
C_r \mid f,N_r,\kappa
\sim
\mathrm{NB}(\mu=fN_r,\mathrm{size}=\kappa),
$$

with

$$
\mathrm{Var}(C_r\mid f) = \mu+\frac{\mu^2}{\kappa}.
$$

Latent frequency is represented on a discrete log-frequency grid. Pair-specific parameters and clonotype-level posterior quantities are inferred from the technical-replicate pair.

The state layer stores, among other quantities:

- latent-frequency summaries;
- latent log-frequency uncertainty;
- replicate-specific observed abundance;
- replicate-specific posterior summaries;
- replicate-posterior overlap;
- endpoint observation-model quantities;
- operational observation state;
- optional full posterior distributions for downstream uncertainty propagation.

Two concepts must remain distinct:

- **posterior/model-based detectability**, which is a property of the fitted observation model;
- **operational observability**, defined from the endpoint `p_value` and a chosen threshold `alpha`.

The primary reference threshold is

```text
alpha = 0.05
```

but this threshold does not alter the latent-frequency fit. Step 16 varies only the operational T/F classification.

No longitudinal dynamical model is imposed during Step 2.

---

# Primary longitudinal estimands

The current pipeline deliberately separates **conditioning representation** from **replicate-resolved dynamical measurement**.

## Step 9 — observed replicate-decoupled forward dynamics

The primary forward estimator is built at a one-week lag from reciprocal cross-replicate folds.

### AB fold

```text
conditioning coordinate:
    x_observed_rep1_t0

displacement:
    dx_observed_rep2

support:
    forward_ab_eligible
```

### BA fold

```text
conditioning coordinate:
    x_observed_rep2_t0

displacement:
    dx_observed_rep1

support:
    forward_ba_eligible
```

The two folds are combined **after abundance binning** with exact equal weight:

$$
\mathrm{AB/BA\ combined} = 0.5\,\mathrm{AB}+0.5\,\mathrm{BA}.
$$

No row-count weighting is used.

This geometry prevents the same observed replicate from simultaneously defining both the conditioning coordinate and the measured displacement.

---

## Step 11 — replicate-consistent fluctuation dynamics

The primary fluctuation estimand is

$$
\mathrm{cross_cov} =
\mathrm{Cov}
\left(
\Delta x^{(1)}_{\mathrm{obs}},
\Delta x^{(2)}_{\mathrm{obs}}
\mid
x_{\mathrm{mid,latent}},
\Delta t,
\mathrm{common4}
\right).
$$

where:

- `dx_observed_rep1` and `dx_observed_rep2` remain replicate-specific;
- `xmid_latent` is the transition-centred latent conditioning coordinate;
- `common4` defines the four-measure support required to evaluate both replicate displacements.

Complementary quantities are

$$
\mathrm{same\_var\_mean} = \frac{\operatorname{Var}\!\left(\Delta x^{(1)}_{\mathrm{obs}}\right) + \operatorname{Var}\!\left(\Delta x^{(2)}_{\mathrm{obs}}\right)}{2}.
$$

and

$$
\mathrm{replicate\_specific\_excess} = \mathrm{same\_var\_mean} - \mathrm{cross_cov}.
$$

`cross_cov` is signed and is never clipped at zero.

The primary Step-11 analysis does **not** impose a TT-only operational-class filter.

---

## Step 13 — temporal scaling

Step 13 evaluates genuine longitudinal `cross_cov` over lags of 1–5 weeks.

The primary analysis uses:

```text
metric:
    cross_cov

abundance support:
    fixed Step-11 abundance core

cohort:
    biological subjects with complete support across all analyzed lags

primary estimator:
    equal-bin within subject
    then equal-subject across the cohort
```

The primary signed descriptor is

$$
M(\Delta t)
=
K
+
D(\Delta t-1),
$$

where \(D\) is the temporal slope per week.

Model comparison and bootstrap uncertainty are evaluated without subtracting the pseudo technical baseline.

---

# Pseudo-longitudinal technical reference

The pseudo reference is not a second biological longitudinal cohort.

It originates from **one biological source sample** processed into 12 independent technical measurements:

```text
one blood sample
    ↓
PBMC isolation
    ↓
12 aliquots
    ↓
independent RNA extraction / RT / library preparation / sequencing
    ↓
12 exchangeable technical repertoire measurements
```

The 12 measurements are repeatedly arranged as:

```text
6 pseudo-timepoints × 2 technical replicates
```

The production ensemble contains:

```text
2,000 pseudo-time configurations
```

## Pair bank

Twelve measurements define

$$
\binom{12}{2}=66
$$

unique unordered technical-replicate pairs.

Because pair-specific Step-2 latent inference depends only on the pair and not on the later pseudo-time assignment, the 66 pair fits are computed once and reused across the 2,000 pseudo configurations.

## Step 7

Step 7 evaluates the pseudo technical null of the observed replicate-decoupled forward estimator.

The same AB/BA logic used downstream in genuine longitudinal data is retained, and the two folds are combined with exact 0.5/0.5 weight after binning.

The pseudo configuration is the randomization unit. It is **not** a biological subject.

## Step 8

Step 8 evaluates technical-reference fluctuation structure using:

```text
support:
    common4

primary conditioning:
    xmid_latent

primary metric:
    cross_cov
```

The current implementation uses a reusable pair-transition bank so that the same pair-vs-pair transition block is not recomputed independently in all 2,000 pseudo configurations.

Randomization intervals across configurations are technical-reference intervals, not biological confidence intervals.

---

# Longitudinal-versus-pseudo comparisons

## Step 10 — forward technical-null comparison

Inputs:

```text
genuine longitudinal Step 9
pseudo-reference Step 7
```

Two comparison geometries are reported.

### Absolute abundance

Genuine and pseudo curves are compared only over shared native observed-log-frequency support.

Pseudo missing bins are not interpolated or imputed.

### Percentile geometry

The comparison is repeated in within-unit abundance-percentile space:

- within subject/fold for the genuine cohort;
- within pseudo configuration/fold for the technical reference.

Pseudo configurations remain randomization units and are never treated as biological replicates.

---

## Step 12 — cross-replicate fluctuation technical-null comparison

Inputs:

```text
genuine longitudinal Step 11
pseudo-reference Step 8
```

The primary comparison uses matched absolute `xmid_latent` support and the same signed `cross_cov` estimand.

Pseudo native bins are not imputed or bridged across missing support.

The comparison is a technical-null benchmark. It does **not** define the temporal accumulation analysis used in Step 13.

---

# Synthetic positive controls

The positive-control workflow is separate from the pseudo-longitudinal technical reference.

Its purpose is to test whether the production pipeline can recover a known injected temporal signal under empirical observation/support constraints.

The empirical analysis is not rebuilt. The validation branch freezes the finalized real-data calibration, including:

- empirical source/support structure;
- Step-2 parameters;
- Step-11 abundance grid;
- Step-13 abundance core;
- Step-13 complete-case subject cohort.

Four production scenarios are generated:

```text
R0p00
R0p25
R0p50
R1p00
```

For each scenario:

```text
synthetic repertoires
    ↓
production Step 2
    ↓
production Step 4
    ↓
production Step 5
    ↓
production Step 11
    ↓
production Step 13
```

The hidden simulation ground truth is stored separately and is never passed to Step 2 or the downstream production pipeline.

Before summarization, `verify_workflow.py` checks that the synthetic analysis preserves the frozen empirical grid, core, support contract, complete-case design, and Step-13 reconstruction.

---

# Robustness and sensitivity controls

## Step 14 — interval position and subject composition

Step 14 asks whether the temporal result from Step 13 could be materially driven by which calendar intervals and subjects contribute at different lags.

It inherits:

- the Step-13 fixed abundance core;
- the Step-13 complete-case subject universe.

It analyzes interval-resolved Step-11 quantities and evaluates:

- synchronized calendar-position structure;
- common-start anchored temporal profiles;
- common-end anchored temporal profiles;
- available-subject analyses;
- matched-all-lags subject analyses;
- abundance-resolved anchored sensitivity.

Step 14 does not redefine the Step-13 primary cohort or abundance core.

---

## Step 15 — fixed-threshold observation-domain sensitivity

Step 15 keeps

```text
alpha = 0.05
```

fixed and asks whether the conclusions change when the operational observation domain is restricted.

Operational endpoint classes are defined by the Step-5 endpoint `p_value` values:

```text
T(alpha) = 1[p_endpoint < alpha]

TT
TF
FT
FF
```

These are **operational observation states**, not biological presence/absence states.

Forward sensitivity evaluates:

```text
PRIMARY
T0PLUS = TT + TF
TT
TF
```

Fluctuation sensitivity evaluates:

```text
PRIMARY
NON_FF = TT + TF + FT
TT
```

Temporal sensitivity is evaluated on a matched common operational abundance core while preserving the frozen Step-13/14 subject universe.

---

## Step 16 — operational-threshold robustness

Step 16 varies the numerical operational threshold:

```text
alpha = 0.10
alpha = 0.05
alpha = 0.025
alpha = 0.01
```

with:

```text
reference alpha = 0.05
```

Only the operational T/F classification changes.

Step 16:

- does not refit latent abundance;
- does not change the Step-9 forward estimand;
- does not change the Step-11 fluctuation estimand;
- does not replace the Step-13 temporal estimator;
- reuses the validated Step-15 analysis implementation across thresholds;
- uses paired biological-subject bootstrap draws across alpha;
- constructs common cross-alpha forward and temporal support.

The finalized alpha=0.05 Step-15 output can also be used as an external replication audit of the Step-16 reference-threshold branch.

---

# Input data

ClonoDynamics starts from processed clonotype-level repertoire tables, not raw FASTQ files.

For core Steps 1–2, each repertoire table must contain at least:

| Column | Meaning |
|---|---|
| `aaSeqCDR3` | amino-acid CDR3 clonotype identifier |
| `readCount` | clonotype sequencing read count |

`readFraction` may be present but is not required if it can be reconstructed from `readCount`.

## Longitudinal filenames

The default longitudinal naming convention is:

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

where the final `1` or `2` identifies the technical replicate.

The current core readers support the text and columnar formats used by the pipeline, including TSV/CSV and Parquet.

---

# Master execution

The recommended entry point is:

```bash
python3 code/orchestrator_clonodynamics.py
```

The master is self-contained and embeds the validated orchestration engines for:

```text
core
pseudo
validation
positive_controls
controls
figures
```

Separate runtime `orchestrator_core.py`, `orchestrator_pseudo.py`, etc. files are not required.

The master first explains the project and then asks whether to enter:

```text
[1] Analysis
[2] Figures
[3] Quit
```

---

## Analysis menu

The analysis interface exposes five blocks.

```text
[1] core
[2] pseudo
[3] validation
[4] positive_controls
[5] controls
[6] all analysis blocks
```

Multiple blocks can also be selected when appropriate.

### Core

Needs:

```text
genuine longitudinal repertoire dataset
destination for core results
reference operational alpha
```

Produces Steps 1–6 directly under the selected results root.

### Pseudo reference

Needs:

```text
12 pseudo technical repertoire measurements
destination for pseudo results
reference operational alpha
```

Builds:

- the pseudo design;
- the 66-pair state bank;
- the 2,000 compact pseudo configurations;
- Step 7;
- the reusable Step-8 transition bank;
- Step 8.

### Validation

Needs:

```text
completed core results
completed pseudo-reference results
destination for validation results
```

Runs Steps 9–13.

Pseudo Step 7 is required for Step 10.

Pseudo Step 8 is required only for Step 12. If pseudo Step 8 is incomplete, the remaining genuine-data branches can still proceed.

### Positive controls

Needs:

```text
corrected longitudinal dataset
completed core results
completed validation results
destination for positive-control results
```

Runs empirical calibration freeze, oracle calibration, synthetic scenario generation, production recovery pipelines, verification, summarization, and plotting.

### Controls

Needs:

```text
completed longitudinal/validation results
destination for control results
```

Runs Steps 14–16.

---

# Resume, repair, and fresh execution

The orchestrators use conservative existing-output policies.

## Resume

`resume` preserves compatible completed output and continues only where needed.

Typical behavior:

```text
complete step  → skip
missing step   → run
interrupted resumable work → reuse compatible checkpoints
```

## Fresh

`fresh` deliberately deletes the outputs owned by that analysis block and reruns the block from the beginning.

Raw input datasets are not deleted.

## Repair

The expensive pseudo-reference upstream additionally supports a repair/resume mode.

This is used when metadata or some configurations are incomplete but the expensive reusable pair bank and compact pseudo configurations should be preserved whenever compatible.

It is specifically intended to avoid unnecessary regeneration of the 66 pair fits and 2,000 pseudo configurations.

---

# Recommended result trees

## Genuine longitudinal results

If the selected core output directory is:

```text
dataset_longitudinal_results/
```

the core writes directly into that root:

```text
dataset_longitudinal_results/
+-- 00_orchestrator_config.json
+-- 00_orchestrator_step_manifest.csv
+-- logs/
|   `-- core/
|
+-- 1-repertoire_characterization/
+-- 2-multirepresentation_clonotype_state_inference/
+-- 3-latent_state_and_observation_model_diagnostics/
+-- 4-multirepresentation_trajectory_assembly/
+-- 5-longitudinal_transition_assembly/
`-- 6-transition_support_characterization/
```

The validation and control blocks can add their outputs to the same root:

```text
dataset_longitudinal_results/
+-- ...
+-- 9-observed_replicate_decoupled_forward_drift/
+-- 10-longitudinal_vs_pseudo_forward_null/
+-- 11-cross_replicate_fluctuation_dynamics/
+-- 12-longitudinal_vs_pseudo_cross_replicate_fluctuations/
+-- 13-temporal_fluctuation_scaling/
+-- 14-interval_position_structure/
+-- 15-detectability_boundary_sensitivity/
+-- 16-observation_threshold_robustness/
`-- positive_controls/
```

No `alpha_*` subdirectory is created by the current core orchestrator.

---

## Pseudo-reference results

A typical pseudo result root is:

```text
dataset_pseudo_results/
+-- source_measurements.tsv
+-- pseudo_configurations.tsv
+-- pseudo_configurations_wide.tsv
+-- 00_pseudo_1x12_config.json
|
+-- pair_bank_multirepresentation/
|   +-- pair_bank_manifest.tsv
|   +-- 00_pairbank_signature.json
|   `-- 2-clonotype_state_inference/
|
+-- ensemble_results_v3/
|   +-- C000001/step7_cache/step7_dt1.parquet
|   +-- ...
|   `-- C002000/step7_cache/step7_dt1.parquet
|
+-- 7-pseudo_forward_technical_null_compact/
+-- step8_transition_bank_v1/
`-- 8-pseudo_fluctuation_conditioning_validation_pairbank/
```

The 2,000 configuration directories retain only the compact Step-7 cache after production cleanup.

---

# Principal outputs

Selected production outputs include:

```text
Step 1
    repertoire_characteristics.csv
    repertoire_tail_model_scan.csv

Step 2
    latent_params_by_pair.csv
    latent_pair_qc.csv
    per_clone_latent_subject.parquet
    per_clone_latent_logP/

Step 3
    clonotype_latent_posterior_qc.parquet
    clonotype_latent_posterior_qc_metric_summary.csv
    clonotype_latent_posterior_qc_by_abundance.csv

Step 4
    multirepresentation_trajectories_long.parquet

Step 5
    longitudinal_transitions.parquet

Step 6
    01_transition_dataset_summary.csv
    08_estimand_support_by_dt.csv
    11_estimand_readiness_audit.csv
    12_structural_consistency_audit.csv

Step 7
    02_ensemble_metric_summary.csv
    04_curve_randomization_envelope.csv

Step 8
    02_ensemble_dt_summary.csv
    03_configuration_binned_metrics_long.csv
    04_binned_randomization_envelope.csv
    06_temporal_slope_randomization_summary.csv
    09_support_by_configuration_dt.csv

Step 9
    02_primary_forward_by_bin.csv
    03_matched_same_vs_cross_by_bin.csv
    04_forward_global_summary.csv
    05_subject_forward_by_bin.csv
    06_loso_forward_by_bin.csv
    07_AB_BA_concordance.csv
    09_bootstrap_linear_descriptors.csv

Step 10
    01_absolute_abundance/
    02_percentile_geometry/

Step 11
    02_fluctuation_by_bin_dt.csv
    03_global_fluctuation_by_dt.csv
    05_subject_interval_fluctuation_by_bin.csv
    08_subject_bin_sufficient.parquet

Step 12
    01_shared_support_by_dt.csv
    02_pointwise_cross_covariance_comparison.csv
    03_pseudo_configuration_matched_support_metrics.csv
    04_empirical_tests_by_dt.csv
    06_joint_comparison.csv
    07b_stable_native_bins_by_dt.csv

Step 13
    01_core_bin_selection.csv
    02_complete_case_subjects.csv
    04_cohort_core_metric_by_dt.csv
    05_temporal_slope_summary.csv
    08_binwise_temporal_slopes.csv
    09_model_comparison.csv
    10_joint_subject_bootstrap.csv
    11_step14_contract.csv

Step 14
    06_core_interval_position_tests.csv
    07_binwise_interval_position_tests.csv
    08_anchored_core_profiles.csv
    09_anchored_core_slope_bootstrap.csv
    10_anchored_subject_slopes.csv
    13_step15_contract.csv

Step 15
    01_operational_class_composition_by_dt.csv
    03_forward_domain_by_bin.csv
    07_fluctuation_domain_by_bin_dt.csv
    11_temporal_common_core.csv
    12_temporal_domain_by_dt.csv
    13_temporal_domain_slope_summary.csv

Step 16
    02_class_composition_pooled_by_alpha.csv
    03_tt_retention_by_alpha.csv
    04_forward_domain_cross_alpha_by_bin.csv
    05_forward_slope_cross_alpha.csv
    08_fluctuation_domain_cross_alpha_by_bin_dt.csv
    11_cross_alpha_temporal_common_core.csv
    12_temporal_profiles_cross_alpha.csv
    13_temporal_slopes_cross_alpha.csv
```

Each analysis directory also contains configuration, signature, audit, provenance, or manifest files appropriate to that step.

---

# Figures

Figure generation is separate from scientific analysis.

From the master:

```bash
python3 code/orchestrator_clonodynamics.py
```

choose:

```text
Figures
```

The figure interface then allows interactive selection of:

```text
core
pseudo
validation
positive_controls
controls
controls_details
```

The canonical plotting workflow is read-only with respect to the finalized analysis tables.

A typical figure root is:

```text
figures/
+-- 01_core/
+-- 02_pseudo_reference/
+-- 03_validation/
`-- 04_controls/
```

The controls block uses the integrated Steps-14–16 robustness plotter by default.

The individual detailed Step-14/15/16 plotters are optional and are not part of the standard canonical figure run.

---

# Reproducibility and provenance

The scientific definitions, estimands, support rules, and block-specific computational procedures are documented in [`docs/methods/`](docs/methods/README.md). Reproducibility therefore relies on the combination of **code + recorded provenance + matching methods documentation**, rather than on script names alone.


ClonoDynamics is designed around explicit intermediate states and dependency-aware reruns.

Depending on the block, the workflow records:

- input paths;
- output paths;
- operational threshold;
- random seeds;
- bootstrap settings;
- pseudo-configuration counts;
- code/script hashes;
- analysis signatures;
- Python executable and version;
- timestamps;
- command lines;
- step-specific logs;
- resumable checkpoint state.

The master embeds the validated orchestration engines and verifies their embedded source hashes before execution.

The individual scientific scripts remain directly executable for development, testing, and methodological inspection.

---

# Interpretation boundaries

ClonoDynamics distinguishes several layers of inference.

## Latent state

The latent representation is a model-based uncertainty-aware representation of clonotype abundance.

It should not be interpreted as a universally preferred direct dynamical measurement.

## Technical replicate measurements

Observed technical replicates remain distinct in the primary forward and fluctuation estimands.

## Operational observation state

T/F and TT/TF/FT/FF are operational classifications defined relative to an observation threshold.

They are not biological appearance/disappearance, birth/extinction, or presence/absence states.

## Pseudo-longitudinal reference

Pseudo-time is an arbitrary technical ordering of exchangeable measurements.

It is a technical null, not biological time.

## Positive controls

Synthetic ground truth is used only to generate and evaluate controlled perturbations.

It is not exposed to the production inference pipeline.

## Temporal scaling

Finite-lag temporal structure is estimated only from genuine longitudinal data.

The pseudo technical baseline is used for comparison, not subtracted from Step 13.

---

# Example data

A lightweight example dataset may be distributed under:

```text
example_repseq/
```

using the same longitudinal technical-replicate naming convention.

Example:

```text
1_1-1.tsv
1_1-2.tsv
1_2-1.tsv
1_2-2.tsv
```

Example data should be treated as workflow-demonstration material and not as a substitute for the full study dataset.

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

A public release should include a pinned environment specification or `requirements.txt`.

---

# Data availability

ClonoDynamics operates on processed clonotype-level repertoire tables.

Final public releases should provide links to:

```text
Raw sequencing:
    [NCBI SRA accession]

Processed repertoire and analysis data:
    [Zenodo / Figshare / institutional repository DOI]

Source repository:
    https://github.com/ClonoDynamics/clonodynamics-pipeline

Software release:
    [GitHub release / Zenodo DOI]
```

---

# Related projects

The methodological output of ClonoDynamics is intended to support downstream quantitative representations of repertoire dynamics.

A separate downstream repository is planned for:

**ClonoDynamics Reference Space**

```text
clonodynamics-reference-space
```

Its scope is distinct from the present repository:

1. quantitative encoding of healthy repertoire dynamics;
2. construction of a healthy physiological dynamic reference space;
3. projection and characterization of perturbed repertoire dynamics relative to that frozen reference.

Keeping this downstream project in a separate repository preserves a clear distinction between:

```text
ClonoDynamics
    repertoire data → validated longitudinal dynamics

ClonoDynamics Reference Space
    validated dynamics → physiological reference geometry → perturbation analysis
```

---

# Citation

If you use ClonoDynamics, please cite the associated methods paper:

```text
ClonoDynamics: replicate-resolved inference and validation
of longitudinal T-cell receptor repertoire dynamics

[full citation / DOI to be added]
```

A `CITATION.cff` file should be included in the public repository once the final citation is available.

---

# Contact

Camillo Palmieri, PhD
cpalmieri@unicz.it

For questions about the software, reproducibility, or methodological implementation, please use the GitHub issue tracker or contact the corresponding author.
