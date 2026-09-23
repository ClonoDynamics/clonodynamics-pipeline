# ClonoDynamics output reference

## Purpose

This document is a practical reference for the principal outputs produced by the current ClonoDynamics workflow.

It is not a replacement for the scientific methods in [`methods/`](methods/README.md).

For each analysis block, the emphasis here is on:

- canonical files;
- row or object granularity;
- scientific role;
- principal downstream consumers.

Diagnostic, audit, manifest, and plotting-source files may exist in addition to the files listed here.

---

# 1. General output conventions

## 1.1 Analysis roots

The genuine longitudinal workflow can use one shared results root, for example:

```text
dataset_longitudinal_results/
```

Core, validation, positive-control, and control outputs may coexist below that root.

The pseudo-reference workflow normally uses a separate root, for example:

```text
dataset_pseudo_results/
```

Publication figures are written to a separate figure root, for example:

```text
figures/
```

## 1.2 Provenance files

Analysis directories may contain files such as:

```text
00_run_config.json
00_run_signature.json
00_input_manifest.csv
00_analysis_summary.csv
00_pipeline_manifest.csv
```

The exact set is step-specific.

These files are part of the provenance record and should be retained with the scientific outputs.

## 1.3 Logs

Orchestrator logs are written below the selected result root.

For the core block, the canonical log location is:

```text
logs/core/
```

Other blocks retain their own logs or manifests as defined by their embedded orchestration engine.

---

# 2. Core results: Steps 1–6

## Step 1 — repertoire characterization

Canonical directory:

```text
1-repertoire_characterization/
```

### `repertoire_characteristics.csv`

**Granularity**

```text
one row per repertoire / technical replicate
```

**Role**

Contains repertoire-level descriptive properties, including depth, richness, and heavy-tail fit summaries.

**Use**

- descriptive cohort characterization;
- Step-1 figures;
- positive-control empirical calibration metadata.

It is not used to filter the primary downstream cohort by tail-fit quality.

### `repertoire_tail_model_scan.csv`

**Granularity**

```text
repertoire × prespecified tail threshold
```

**Role**

Stores the Pareto versus lower-truncated-log-normal tail-model comparison across threshold quantiles.

**Use**

- descriptive robustness analysis;
- Step-1 model-comparison figures.

---

## Step 2 — multirepresentation clonotype-state inference

Canonical directory:

```text
2-multirepresentation_clonotype_state_inference/
```

### `latent_params_by_pair.csv`

**Granularity**

```text
one row per successfully fitted technical-replicate pair
```

**Role**

Stores pair-specific fitted latent-model parameters.

**Use**

- model provenance;
- empirical calibration for positive controls;
- pair-level diagnostic inspection.

### `latent_pair_qc.csv`

**Granularity**

```text
one row per attempted or fitted replicate pair
```

**Role**

Documents pair-level inference/QC status and fitted-pair metadata.

### `per_clone_latent_subject.parquet`

**Granularity**

```text
subject × time × clonotype
```

for successfully fitted complete replicate pairs.

**Role**

Canonical Step-2 state table.

Contains:

- latent-state summaries;
- replicate-resolved observed quantities;
- operational observation statistics;
- posterior/detectability diagnostics;
- fields required by Step 4.

**Primary downstream consumer**

```text
Step 4
```

### `per_clone_latent_logP/`

**Object type**

Directory containing per-pair full posterior caches.

**Role**

Optional full-posterior endpoint uncertainty representation.

**Primary downstream consumer**

```text
Step 5 posterior propagation
```

The primary Step-9 and Step-11 observed estimands do not require these posterior caches.

---

## Step 3 — latent-state and observation-model diagnostics

Canonical directory:

```text
3-latent_state_and_observation_model_diagnostics/
```

### `clonotype_latent_posterior_qc.parquet`

**Granularity**

```text
clonotype-state diagnostic row
```

**Role**

Detailed posterior-quality and observation-model diagnostic table.

### `clonotype_latent_posterior_qc_metric_summary.csv`

**Role**

Summary of key diagnostic metrics across the analyzed state table.

### `clonotype_latent_posterior_qc_by_abundance.csv`

**Role**

Abundance-resolved diagnostic summaries.

**Use**

Step-3 figures and methodological diagnostics.

Step 3 is not a filtering stage for the primary downstream analyses.

---

## Step 4 — multirepresentation trajectory assembly

Canonical directory:

```text
4-multirepresentation_trajectory_assembly/
```

### `multirepresentation_trajectories_long.parquet`

**Granularity**

```text
subject × clonotype × time
```

**Role**

Canonical longitudinal state table carrying latent, observed replicate-resolved, uncertainty, and observation-model fields.

**Primary downstream consumer**

```text
Step 5
```

---

## Step 5 — longitudinal transition assembly

Canonical directory:

```text
5-longitudinal_transition_assembly/
```

### `longitudinal_transitions.parquet`

**Granularity**

```text
subject × clonotype × ordered interval (t0, t1)
```

**Role**

Canonical finite-time transition universe for downstream ClonoDynamics analyses.

Contains, among other fields:

```text
t0
t1
dt

x0_latent
x1_latent
xmid_latent
xstar_latent
dx_latent

x_observed_rep1_t0
x_observed_rep2_t0
dx_observed_rep1
dx_observed_rep2

forward_ab_eligible
forward_ba_eligible
common4

p_value_t0
p_value_t1
reference operational class
```

**Primary downstream consumers**

```text
Step 6
Step 9
Step 11
Step 15
Step 16
```

### `transition_assembly_report.md`

**Role**

Human-readable transition-assembly summary.

### `transition_assembly_manifest.json`

**Role**

Machine-readable Step-5 construction/provenance record.

---

## Step 6 — transition-support characterization

Canonical directory:

```text
6-transition_support_characterization/
```

### `01_transition_dataset_summary.csv`

**Role**

Global transition-universe summary.

### `08_estimand_support_by_dt.csv`

**Granularity**

```text
temporal lag
```

**Role**

Summarizes support for downstream estimator-specific domains such as AB, BA, and `common4`.

### `11_estimand_readiness_audit.csv`

**Role**

Checks whether the fields required by the intended downstream estimands are structurally available on their declared support.

### `12_structural_consistency_audit.csv`

**Role**

Checks internal consistency of transition fields, support flags, endpoint positivity, lag definitions, and operational annotations.

Step 6 is descriptive/readiness-oriented and does not estimate longitudinal dynamics.

---

# 3. Pseudo-reference outputs: Steps 7–8

The pseudo-reference result root is independent of the genuine longitudinal result root.

## Reusable upstream pseudo objects

### Pseudo design tables

The pseudo workflow stores design/provenance objects defining the 12 source measurements and the 2,000 randomized configurations.

### Pair bank

Canonical reusable pair-bank root:

```text
pair_bank_multirepresentation/
```

Important objects include:

```text
pair_bank_manifest.tsv
00_pairbank_signature.json
```

The pair bank represents the 66 unique unordered pairs of the 12 technical source measurements.

### Compact pseudo ensemble

Canonical compact ensemble root:

```text
ensemble_results_v3/
```

Each completed configuration retains its compact Step-7 cache rather than a full duplicated upstream pipeline.

---

## Step 7 — pseudo forward technical null

Canonical directory:

```text
7-pseudo_forward_technical_null_compact/
```

### `02_ensemble_metric_summary.csv`

**Role**

Across-configuration summaries of Step-7 forward technical-null descriptors.

### `04_curve_randomization_envelope.csv`

**Granularity**

```text
representation / section × abundance bin
```

summarized across pseudo configurations.

**Role**

Canonical pseudo randomization envelope used for Step-7 figures and Step-10 comparison.

**Primary downstream consumer**

```text
Step 10
```

---

## Step 8 — pseudo fluctuation technical reference

Canonical directory:

```text
8-pseudo_fluctuation_conditioning_validation_pairbank/
```

### `02_ensemble_dt_summary.csv`

**Granularity**

```text
pseudo lag
```

**Role**

Across-configuration pseudo-lag summary of fluctuation quantities.

### `03_configuration_binned_metrics_long.csv`

**Granularity**

```text
pseudo configuration × pseudo lag × abundance bin
```

**Role**

Configuration-resolved binned pseudo fluctuation metrics.

**Primary downstream consumer**

```text
Step 12
```

### `04_binned_randomization_envelope.csv`

**Role**

Across-configuration randomization envelope by pseudo lag and abundance bin.

**Use**

- Step-8 figures;
- pseudo-reference visualization;
- Step-12 support/comparison context.

### `06_temporal_slope_randomization_summary.csv`

**Role**

Pseudo-lag slope distribution summary.

Pseudo lag is not biological time.

### `09_support_by_configuration_dt.csv`

**Granularity**

```text
pseudo configuration × pseudo lag
```

**Role**

Support audit used by plotting and downstream pseudo-reference validation.

### Reusable Step-8 transition bank

Canonical reusable bank root:

```text
step8_transition_bank_v1/
```

This bank prevents recomputation of identical pair-vs-pair transition blocks across pseudo configurations.

---

# 4. Validation outputs: Steps 9–13

## Step 9 — observed replicate-decoupled forward dynamics

Canonical directory:

```text
9-observed_replicate_decoupled_forward_drift/
```

### `02_primary_forward_by_bin.csv`

**Granularity**

```text
forward representation/fold × abundance bin
```

**Role**

Canonical binned AB, BA, and equal-weight combined forward profiles.

**Primary downstream consumer**

```text
Step 10
```

### `03_matched_same_vs_cross_by_bin.csv`

**Role**

Matched `common4` same-measure versus cross-replicate coupling control.

### `04_forward_global_summary.csv`

**Role**

Global descriptive summary of the forward analysis.

### `05_subject_forward_by_bin.csv`

**Granularity**

```text
subject × fold/representation × abundance bin
```

**Role**

Subject-resolved forward profiles.

### `06_loso_forward_by_bin.csv`

**Granularity**

```text
leave-one-subject-out run × abundance bin
```

**Role**

LOSO robustness diagnostic.

### `07_AB_BA_concordance.csv`

**Role**

AB-versus-BA agreement summary.

### `09_bootstrap_linear_descriptors.csv`

**Role**

Subject-bootstrap distributions for compact linear descriptors such as slope and zero crossing.

---

## Step 10 — longitudinal versus pseudo forward null

Canonical directory:

```text
10-longitudinal_vs_pseudo_forward_null/
```

Two principal result subtrees are:

```text
01_absolute_abundance/
02_percentile_geometry/
```

### `01_absolute_abundance/`

**Role**

Compares finalized longitudinal Step-9 forward structure with the pseudo Step-7 technical null on shared native observed-abundance support.

Pseudo missing bins are not imputed.

### `02_percentile_geometry/`

**Role**

Performs the complementary within-unit abundance-percentile comparison.

Longitudinal uncertainty is biological subject-cluster uncertainty.

Pseudo intervals describe the randomization distribution across pseudo configurations.

---

## Step 11 — cross-replicate fluctuation dynamics

Canonical directory:

```text
11-cross_replicate_fluctuation_dynamics/
```

### `02_fluctuation_by_bin_dt.csv`

**Granularity**

```text
temporal lag × abundance bin
```

**Role**

Canonical cohort-level abundance-resolved fluctuation profiles.

Principal quantities include:

```text
cross_cov
same_var_mean
replicate_specific_excess
```

`cross_cov` is signed.

**Primary downstream consumers**

```text
Step 12
Step 13
```

### `03_global_fluctuation_by_dt.csv`

**Granularity**

```text
temporal lag
```

**Role**

Lag-level global fluctuation summaries.

### `05_subject_interval_fluctuation_by_bin.csv`

**Granularity**

```text
subject × physical interval × abundance bin
```

**Role**

Interval-resolved fluctuation estimates.

**Primary downstream consumer**

```text
Step 14
```

### `08_subject_bin_sufficient.parquet`

**Granularity**

```text
subject × temporal lag × abundance bin
```

with sufficient-statistic fields.

**Role**

Compact subject-level representation for pooled reconstruction and downstream temporal/comparison analyses.

---

## Step 12 — longitudinal versus pseudo cross-replicate fluctuation comparison

Canonical directory:

```text
12-longitudinal_vs_pseudo_cross_replicate_fluctuations/
```

### `01_shared_support_by_dt.csv`

**Role**

Matched longitudinal/pseudo abundance support by lag.

### `02_pointwise_cross_covariance_comparison.csv`

**Granularity**

```text
lag × matched abundance point/bin
```

**Role**

Pointwise comparison of signed longitudinal and pseudo cross-replicate covariance.

### `03_pseudo_configuration_matched_support_metrics.csv`

**Granularity**

```text
pseudo configuration × lag
```

**Role**

Configuration-resolved pseudo metrics evaluated over the matched support.

### `04_empirical_tests_by_dt.csv`

**Granularity**

```text
lag
```

**Role**

Lag-specific empirical longitudinal-versus-pseudo comparison summaries.

### `06_joint_comparison.csv`

**Role**

Joint comparison summary across the analyzed lags/support.

### `07b_stable_native_bins_by_dt.csv`

**Role**

Stable pseudo-native support retained for matched comparison.

Step 12 does not feed Step 13 and is not subtracted from the genuine longitudinal covariance.

---

## Step 13 — temporal fluctuation scaling

Canonical directory:

```text
13-temporal_fluctuation_scaling/
```

### `01_core_bin_selection.csv`

**Role**

Records the fixed abundance core selected for temporal analysis.

### `02_complete_case_subjects.csv`

**Role**

Records the biological-subject cohort with complete support over the selected core and lags.

### `04_cohort_core_metric_by_dt.csv`

**Granularity**

```text
metric/estimator × temporal lag
```

**Role**

Canonical fixed-core cohort temporal profile.

### `05_temporal_slope_summary.csv`

**Role**

Primary and sensitivity signed temporal slope summaries with bootstrap uncertainty.

### `08_binwise_temporal_slopes.csv`

**Granularity**

```text
abundance bin
```

**Role**

Secondary abundance-resolved temporal slope analysis.

### `09_model_comparison.csv`

**Role**

Finite-lag model-comparison summary.

### `10_joint_subject_bootstrap.csv`

**Granularity**

```text
bootstrap draw × metric/estimator
```

**Role**

Joint subject-cluster bootstrap output propagated across lags.

### `11_step14_contract.csv`

**Role**

Frozen contract passed to Step 14, including the Step-13 abundance core and subject universe required for the interval/composition control.

---

# 5. Positive-control outputs

Canonical positive-control root:

```text
positive_controls/
```

The workflow contains calibration, oracle, scenario, verification, summary, and figure outputs.

## Calibration

The calibration directory freezes the finalized empirical design required to generate support-conditioned synthetic controls.

It contains, among other objects:

- eligible sampling design;
- empirical pair parameters;
- support targets;
- subject-level reference repertoires;
- empirical Step-11 abundance geometry;
- Step-13 core and cohort contract.

## Oracle

Canonical oracle outputs include:

```text
00_qbio_oracle_calibration.json
01_scenario_table.csv
02_oracle_temporal_profiles.csv
03_search_trace.csv
04_manifest.json
```

`01_scenario_table.csv` is the authoritative record of calibrated scenario parameters.

## Scenario roots

Canonical scenarios are:

```text
R0p00
R0p25
R0p50
R1p00
```

Each scenario contains:

```text
repertoires/
ground_truth/
analysis outputs for Steps 2, 4, 5, 11, and 13
```

`ground_truth/` is not an input to production state inference.

## Verification

Canonical verification outputs document agreement with the frozen empirical support and matched bootstrap design.

## Summary

Canonical summary files are:

```text
00_validation_summary.csv
01_temporal_profiles.csv
02_incremental_recovery_summary.csv
03_recovery_fit_summary.csv
04_validation_metadata.json
```

### `00_validation_summary.csv`

**Role**

Scenario-level compact validation summary.

### `01_temporal_profiles.csv`

**Role**

Oracle and recovered temporal profiles used for the validation figures.

### `02_incremental_recovery_summary.csv`

**Role**

Paired baseline-adjusted recovery across positive doses.

### `03_recovery_fit_summary.csv`

**Role**

Descriptive oracle-versus-recovered recovery-fit summaries.

---

# 6. Control outputs: Steps 14–16

## Step 14 — interval-position structure

Canonical directory:

```text
14-interval_position_structure/
```

### `06_core_interval_position_tests.csv`

**Role**

Core-level synchronized calendar-position permutation tests.

### `07_binwise_interval_position_tests.csv`

**Role**

Abundance-resolved interval-position tests.

### `08_anchored_core_profiles.csv`

**Role**

Common-start and common-end lag profiles over the frozen Step-13 core.

### `09_anchored_core_slope_bootstrap.csv`

**Role**

Subject-cluster bootstrap distributions for anchored temporal slopes.

### `10_anchored_subject_slopes.csv`

**Role**

Subject-specific anchored slope diagnostics.

### `13_step15_contract.csv`

**Role**

Frozen support/cohort contract inherited by Step 15.

---

## Step 15 — fixed-threshold observation-domain sensitivity

Canonical directory:

```text
15-detectability_boundary_sensitivity/
```

### `01_operational_class_composition_by_dt.csv`

**Role**

TT/TF/FT/FF composition at the fixed reference threshold.

### `03_forward_domain_by_bin.csv`

**Role**

Forward profiles for the primary and operationally restricted domains.

### `07_fluctuation_domain_by_bin_dt.csv`

**Role**

Abundance- and lag-resolved fluctuation quantities for primary, NON_FF, and TT domains.

### `11_temporal_common_core.csv`

**Role**

Matched operational abundance core used for Step-15 temporal sensitivity.

### `12_temporal_domain_by_dt.csv`

**Role**

Lag profiles by operational domain on the matched sensitivity support.

### `13_temporal_domain_slope_summary.csv`

**Role**

Signed temporal slope summaries by operational domain.

---

## Step 16 — observation-threshold robustness

Canonical directory:

```text
16-observation_threshold_robustness/
```

Alternative compatible directory names may exist in historical result trees, but the public repository should use one canonical name consistently.

### `02_class_composition_pooled_by_alpha.csv`

**Role**

Pooled operational-class composition across analyzed alpha values.

### `03_tt_retention_by_alpha.csv`

**Role**

Retention of TT support relative to the reference threshold.

### `04_forward_domain_cross_alpha_by_bin.csv`

**Role**

Cross-alpha forward profiles on harmonized support.

### `05_forward_slope_cross_alpha.csv`

**Role**

Cross-alpha forward slope summaries.

### `08_fluctuation_domain_cross_alpha_by_bin_dt.csv`

**Role**

Cross-alpha fluctuation-domain profiles.

### `11_cross_alpha_temporal_common_core.csv`

**Role**

Common operational abundance core valid across all analyzed thresholds and lags.

### `12_temporal_profiles_cross_alpha.csv`

**Role**

Cross-alpha temporal profiles on the common support.

### `13_temporal_slopes_cross_alpha.csv`

**Role**

Cross-alpha temporal slope summaries and paired reference-threshold comparisons.

---

# 7. Figure outputs

Figure generation is read-only with respect to the finalized analysis results.

A typical figure root is:

```text
figures/
+-- 01_core/
+-- 02_pseudo_reference/
+-- 03_validation/
`-- 04_controls/
```

## Core figures

Canonical plot groups:

```text
Step 1
Step 3
Step 6
```

## Pseudo figures

Canonical plot group:

```text
Steps 7–8
```

## Validation figures

Canonical plot groups:

```text
Steps 9–10
Steps 11–12
Step 13
positive controls
```

## Controls figures

The standard controls figure workflow uses the integrated Steps-14–16 robustness plotter.

The individual Step-14, Step-15, and Step-16 detail plotters are optional and are not part of the standard canonical figure run.

Plotters may create subdirectories such as:

```text
main_figure/
supplementary_figure/
figure_source_data/
```

depending on the figure block.

---

# 8. Which outputs should downstream projects consume?

Downstream projects should prefer finalized, documented canonical outputs rather than private intermediate caches.

Typical stable interfaces are:

| Purpose | Preferred ClonoDynamics output |
|---|---|
| clonotype state representation | Step 2 `per_clone_latent_subject.parquet` |
| longitudinal state trajectories | Step 4 `multirepresentation_trajectories_long.parquet` |
| generic finite-time transitions | Step 5 `longitudinal_transitions.parquet` |
| replicate-decoupled forward profile | Step 9 `02_primary_forward_by_bin.csv` |
| abundance-resolved fluctuation profile | Step 11 `02_fluctuation_by_bin_dt.csv` |
| subject-level fluctuation information | Step 11 `08_subject_bin_sufficient.parquet` |
| fixed temporal abundance core | Step 13 `01_core_bin_selection.csv` |
| complete-case temporal cohort | Step 13 `02_complete_case_subjects.csv` |
| cohort temporal profile | Step 13 `04_cohort_core_metric_by_dt.csv` |
| temporal slope summary | Step 13 `05_temporal_slope_summary.csv` |

A downstream project such as **ClonoDynamics Reference Space** should explicitly record which of these interfaces it consumes and from which ClonoDynamics release.
