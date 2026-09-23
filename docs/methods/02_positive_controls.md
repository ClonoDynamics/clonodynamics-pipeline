# Synthetic positive controls

## Support-conditioned temporal-signal recovery in ClonoDynamics

This document describes the semi-synthetic positive-control framework used to test whether the finalized ClonoDynamics longitudinal analysis can recover a prescribed temporal accumulation signal under the empirical observation/support structure.

The mathematical specification below preserves the **reported-study validation design**. Numerical cohort sizes, selected core bins, scenario parameters, seeds, and recovery summaries are therefore study-specific frozen quantities rather than general requirements of ClonoDynamics. The current GitHub execution section that follows Part I reflects the reorganized repository and unified master launcher.

For routine execution, use:

```bash
python3 code/orchestrator_clonodynamics.py
```

then select **Analysis → positive_controls**.

Current code mapping:

- production Steps 2, 4 and 5: `code/01_core/`
- production Steps 11 and 13: `code/03_validation/`
- positive-control builders, verifier, summarizer and dedicated plotter: `code/03_validation/positive_controls/`

[Back to the ClonoDynamics README](../../README.md)

---

# Part I. Mathematical model and recovery analysis

## 1. Objective and empirical conditioning

We assessed whether the longitudinal analysis could recover a prescribed increase in clonotype displacement variability with elapsed time under the observation structure of the empirical study. The controls were semi-synthetic: underlying frequency trajectories and read counts were generated, whereas the sampling design, replicate-specific sequencing depths and clonotype detection patterns were inherited from the empirical data. Consequently, the experiment evaluated recovery conditional on observed support, rather than jointly simulating biological abundance and the process by which clonotypes become detectable.

The empirical design comprised 10 participants sampled at six nominal weekly time points, with two technical replicates per available visit. The 114 available repertoires formed 57 complete replicate pairs. Eligibility for downstream analysis was recorded separately from physical sample availability: 56 pairs were eligible after state inference, defining 130 within-subject sampling-time intervals. A pair excluded because of unsuccessful inference was not interpreted as a biologically absent visit. Its raw repertoire files remained available for construction of the subject-level reference, but it did not contribute synthetic count repertoires to the analysis-eligible design.

Let $s$ index subjects, $i$ clonotypes, $t\in\{1,\ldots,6\}$ nominal sampling times and $r\in\{A,B\}$ technical replicates. All logarithms are natural. We distinguish empirical counts $C^{\mathrm{emp}}_{sir}(t)$, simulated frequencies $f^{(d)}_{si}(t)$ in scenario $d$, and synthetic counts $C^{(d)}_{sir}(t)$. Superscript $d=0$ denotes the non-accumulating plateau control.

## 2. Subject-level reference repertoires and frozen calibration

For each subject, the reference universe $\mathcal{I}_s$ was the union of amino-acid CDR3 clonotypes observed in all available empirical repertoire files. Repeated `aaSeqCDR3` rows were aggregated within each repertoire by summing their read counts. No never-observed clonotypes were introduced, and temporal order was not used to construct reference frequencies.

If $\mathcal{R}_s$ denotes the available empirical repertoires of subject $s$, $m_s=|\mathcal{R}_s|$ and $N_{sr}(t)$ their sequencing depths, the reference was

$$
\bar f_{si}
=\frac{1}{m_s}\sum_{(t,r)\in\mathcal{R}_s}
\frac{C^{\mathrm{emp}}_{sir}(t)}{N_{sr}(t)},
\qquad
f^{\mathrm{ref}}_{si}
=\frac{\bar f_{si}}{\sum_{j\in\mathcal{I}_s}\bar f_{sj}}.
$$

An unobserved clonotype contributed zero to a repertoire-specific frequency. This construction gave equal weight to repertoires, rather than weighting them by sequencing depth. The pooled-count reference exported by the calibration builder was a separate quantity and was not the reference used for this dose series.

The frozen calibration recorded the analysis-eligibility mask, replicate depths, clonotype-positive sets and pair-specific dispersion parameters, together with the empirical abundance-bin grid and temporal-analysis core. The fitted Step-2 parameter `k` was retained as $\kappa_{st}$, the Gamma shape parameter used in synthetic count allocation. Empirical temporal estimates were retained as documentary fields only; they were not targets for selecting the cumulative signal or adjusting its recovery.

## 3. Simulated frequency trajectories

For each subject and clonotype, scenario-specific log weights were generated as

$$
\ell^{(d)}_{si}(t)
=\log f^{\mathrm{ref}}_{si}
+\sigma_{\mathrm{fast}}z_{si}(t)
+\sqrt{q_d}\,W_{si}(t),
$$

where $z_{si}(t)$ are independent standard-normal visit-specific shocks before frequency normalization. The unit random walk started at the first nominal time with $W_{si}(1)=0$ and evolved according to

$$
W_{si}(t_j)
=W_{si}(t_{j-1})
+\sqrt{t_j-t_{j-1}}\,\varepsilon_{si}(t_j),
\qquad
\varepsilon_{si}(t_j)\sim\mathcal{N}(0,1).
$$

The fast shocks and random-walk increments were generated separately. The fast component was fixed at $\sigma_{\mathrm{fast}}=0.20$ in every scenario; $q_d$, corresponding to `q_bio`, controlled the variance of cumulative log-weight increments per week. Frequencies were normalized independently at each subject and time point:

$$
f^{(d)}_{si}(t)
=\frac{\exp\{\ell^{(d)}_{si}(t)\}}
{\sum_{j\in\mathcal{I}_s}\exp\{\ell^{(d)}_{sj}(t)\}}.
$$

The implementation used a stable log-softmax calculation. Trajectories were generated at all six nominal times, including times that were not eligible for count generation. Only eligible visits were subsequently used to create observed synthetic repertoires.

### 3.1 Common random numbers

The master simulation seed was `20260918`. Fast shocks and unit random-walk increments were determined by subject and time using NumPy `SeedSequence` inputs `[master, subject, 117031, 501, time]` and `[master, subject, 117031, 502, time]`, respectively. The same random fields and reference-clonotype ordering were retained across scenarios. Changing $q_d$ therefore rescaled the same unit random-walk path rather than drawing a new path; the plateau and positive-dose scenarios shared identical fast shocks.

Count-generation streams were initialized separately using `[master, subject, 910247, 500]`. Within each subject, eligible visits were processed in time order and the two replicates sequentially. Separate Gamma and multinomial variates were used for the two replicates. Reusing a count-stream initialization across doses does not imply that every subsequent count-generation variate is matched across scenarios; the explicitly matched random fields are the fast shocks and unit random-walk paths.

## 4. Oracle definition and dose calibration

The known simulated frequencies provided an oracle reference computed without using counts or downstream recovery estimates. For every analysis-eligible interval $(t_0,t_1)$, the true displacement was

$$
\Delta x^{\mathrm{true},(d)}_{si}(t_0,t_1)
=\log f^{(d)}_{si}(t_1)-\log f^{(d)}_{si}(t_0).
$$

Its population variance across the full subject reference universe, of size $J_s=|\mathcal{I}_s|$, was

$$
v^{(d)}_s(t_0,t_1)
=\frac{1}{J_s}\sum_{i\in\mathcal{I}_s}
\left[
\Delta x^{\mathrm{true},(d)}_{si}(t_0,t_1)
-\overline{\Delta x}^{\mathrm{true},(d)}_s(t_0,t_1)
\right]^2.
$$

This corresponds to `ddof=0`. Let $\mathcal{P}_s(h)$ be the eligible intervals of subject $s$ with lag $h=t_1-t_0$, and let $\mathcal{S}_h$ contain subjects with at least one such interval. The oracle profile was

$$
V_d(h)
=\frac{1}{|\mathcal{S}_h|}
\sum_{s\in\mathcal{S}_h}
\frac{1}{|\mathcal{P}_s(h)|}
\sum_{(t_0,t_1)\in\mathcal{P}_s(h)}
v^{(d)}_s(t_0,t_1).
$$

Thus interval variances were first averaged within subject and lag, then averaged equally across subjects available at that lag. The oracle used analysis-eligible visits but all reference clonotypes; it was not restricted to `common4` transitions, inferred-abundance bins or the nine-subject complete-case recovery cohort.

Positive doses targeted the relative increase

$$
R_d=\frac{V_d(5)-V_d(1)}{V_d(1)}
$$

at values of 0.25, 0.50 and 1.00. For each target, the calibrator bracketed $q_d$, expanded the upper bound when required, and applied up to 22 bisection iterations with an absolute target tolerance of $5\times10^{-4}$. The candidate with the smallest absolute target error was retained. Common random fields were fixed throughout this search. Synthetic counts, inferred states, empirical cross-replicate covariance and recovered temporal slopes were not used to select the dose.

For the zero-accumulation scenario, $q_0$ was fixed at zero without tuning its finite-sample oracle profile to have an exactly zero slope or exactly zero realized $R$. This plateau retained the non-accumulating fast component; it was distinct from the generator's strict `null` preset, which sets both components to zero.

| Scenario | Target $R$ | $\sigma_{\mathrm{fast}}$ | Calibrated $q_d$ |
|:---|---:|---:|---:|
| `R0p00` | 0.00 | 0.20 | 0 |
| `R0p25` | 0.25 | 0.20 | 0.005333333333333334 |
| `R0p50` | 0.50 | 0.20 | 0.01142857142857143 |
| `R1p00` | 1.00 | 0.20 | 0.026666666666666672 |

These are the parameter values recorded for the reported series. A new reproduction reads the values generated in `oracle/01_scenario_table.csv`, rather than inserting rounded values from this table. Target and realized $R$, oracle profiles and the search trace are exported separately. The oracle temporal slope, $\beta_d^{\mathrm{oracle}}$, was an unweighted least-squares linear fit of $V_d(h)$ across $h=1,\ldots,5$.

## 5. Read-count generation conditional on empirical support

For each eligible subject, visit and replicate, let

$$
A_{str}=\{i:C^{\mathrm{emp}}_{sir}(t)>0\}
$$

be the frozen empirical-positive clonotype set and let $K_{str}=|A_{str}|$. Given simulated frequencies, independent Gamma factors were drawn for clonotypes in this set:

$$
G^{(d)}_{sir}(t)
\sim\operatorname{Gamma}
\left(\kappa_{st},\ \mathrm{scale}=1/\kappa_{st}\right),
\qquad i\in A_{str}.
$$

Frequency-weighted allocation probabilities were then defined by

$$
p^{(d)}_{sir}(t)
=\frac{f^{(d)}_{si}(t)G^{(d)}_{sir}(t)}
{\sum_{j\in A_{str}}f^{(d)}_{sj}(t)G^{(d)}_{sjr}(t)}.
$$

Each empirical-positive clonotype received one mandatory read. The remaining reads were allocated as

$$
\begin{aligned}
\left(M^{(d)}_{sir}(t)\right)_{i\in A_{str}}
&\sim\operatorname{Multinomial}
\left(N_{sr}(t)-K_{str},
\left(p^{(d)}_{sir}(t)\right)_{i\in A_{str}}\right),\\[4pt]
C^{(d)}_{sir}(t)
&=\begin{cases}
1+M^{(d)}_{sir}(t), & i\in A_{str},\\
0, & i\notin A_{str}.
\end{cases}
\end{aligned}
$$

The generator required $0<K_{str}\leq N_{sr}(t)$ and a finite positive $\kappa_{st}$. If the sum of Gamma-weighted frequencies was non-finite or non-positive, the numerical safeguard used the unperturbed simulated frequencies within the active set before normalization.

This construction preserved the exact sequencing depth, observed richness and clonotype-specific positivity mask of each empirical replicate. It was not an unconditioned draw of independent negative-binomial counts. Original empirical count vectors were not reused as synthetic observations, although their frequencies, fitted dispersion and positive sets informed the calibration.

Each scenario comprised 112 synthetic repertoires from the 56 eligible replicate pairs. Count outputs were written to `repertoires/` with `aaSeqCDR3`, `readCount` and recalculated `readFraction` columns. Known log-frequency trajectories were written separately to `ground_truth/` and never supplied to production state inference.

### 5.1 Support verification

For each subject and clonotype, let $k_{si}$ be the number of eligible visits positive in at least one replicate, and $b_{si}$ the number positive in both replicates. The generator checked the implied all-interval and four-measure support totals using

$$
N_{\mathrm{transitions}}
=\sum_s\sum_i\binom{k_{si}}{2},
\qquad
N_{\mathrm{common4}}
=\sum_s\sum_i\binom{b_{si}}{2}.
$$

For the six-timepoint design, all visit pairs have lags within 1–5 weeks. The frozen totals were 2,822,922 transitions and 1,182,069 `common4` transitions. Replicate-positive, union-positive and both-positive state counts were also checked against the calibration. These generator-level checks were kept distinct from checks on the support actually retained after synthetic state inference and downstream analysis.

## 6. Recovery through the longitudinal-analysis pipeline

Synthetic count repertoires were processed using the production state-inference, trajectory-assembly and transition-assembly workflow. State inference estimated abundances from the synthetic technical-replicate counts; it did not read the oracle trajectories. The downstream transition table retained replicate-specific counts, depths and observed log-frequencies alongside the model-based state representation.

For positive counts, the observed log-frequency and displacement were

$$
x^{\mathrm{obs},(d)}_{sir}(t)
=\log\left[\frac{C^{(d)}_{sir}(t)}{N_{sr}(t)}\right],
\qquad
\Delta x^{\mathrm{obs},(d)}_{sir}(t_0,t_1)
=x^{\mathrm{obs},(d)}_{sir}(t_1)-x^{\mathrm{obs},(d)}_{sir}(t_0).
$$

Fluctuation recovery used `common4`: both replicates had to contain positive counts at both interval endpoints. The conditioning coordinate was the midpoint of the inferred endpoint log-frequencies,

$$
x^{(d)}_{\mathrm{mid},si}(t_0,t_1)
=\frac{\widehat{x}^{(d)}_{si}(t_0)
+\widehat{x}^{(d)}_{si}(t_1)}{2}.
$$

The empirical 20-bin abundance grid was retained across scenarios. The reproduction configuration defines this grid with equal-width binning and the empirical 0.01 and 0.99 quantile limits; synthetic Step 11 reads the exported empirical edges instead of estimating new ones.

Within subject $s$, lag $h$ and bin $b$, let $\mathcal{T}^{(d)}_{shb}$ contain the eligible transitions and let $n^{(d)}_{shb}$ denote their count. Each transition is one clonotype–interval record. The subject-level cross-replicate covariance was

$$
\widehat{C}^{(d)}_{shb}
=\frac{1}{n^{(d)}_{shb}-1}
\sum_{u\in\mathcal{T}^{(d)}_{shb}}
\left(\Delta x^{(d)}_{uA}-\overline{\Delta x}^{(d)}_{A,shb}\right)
\left(\Delta x^{(d)}_{uB}-\overline{\Delta x}^{(d)}_{B,shb}\right).
$$

Both means were computed on the same transition set. Covariance was calculated across transitions pooled within subject, lag and bin, rather than by averaging interval-specific covariances. It remained signed and was not truncated at zero. No pseudo-longitudinal covariance was subtracted, and operational TT/TF/FT/FF classifications did not define the primary recovery population.

## 7. Fixed abundance core and temporal summaries

The temporal core was selected from the empirical grid before analysis of the synthetic scenarios. At each lag, a subject–lag–bin cell was valid when its cross-replicate covariance was finite and at least two transitions contributed. A bin qualified if it contained valid cells for at least

$$
\max\left\{2,\left\lceil0.90\,n_{\mathrm{available}}(h)\right\rceil\right\}
$$

subjects at every requested lag. The largest contiguous qualifying bin run defined the core; ties were resolved in favour of the run starting at lower abundance. The complete-case cohort comprised subjects with a valid primary covariance in every selected bin at every lag.

The reported analysis used 17 bins, IDs 3–19, spanning $x_{\mathrm{mid}}=-12.4129$ to $-8.5494$, and nine complete-case subjects: 1, 2, 3, 4, 5, 7, 8, 9 and 10. These same bins were supplied to synthetic Step 13. Complete-case subjects were nevertheless determined from each scenario's outputs and checked for agreement with the empirical cohort; they were not assumed to be unchanged merely because the detection masks were fixed. Membership of individual transitions in inferred-abundance bins was allowed to differ across doses.

The minimum of two transitions for a subject-level covariance cell is distinct from the Step-11 pooled-bin reporting threshold of 50 transitions. The latter was not imposed as a minimum of 50 transitions in every subject–lag–bin cell.

Let $\mathcal{B}$ be the fixed core and $\mathcal{S}$ the complete-case subject set. The primary temporal summary gave equal weight to bins within subject, then equal weight to subjects:

$$
m^{(d)}_s(h)
=\frac{1}{|\mathcal{B}|}\sum_{b\in\mathcal{B}}\widehat{C}^{(d)}_{shb},
\qquad
M_d(h)=\frac{1}{|\mathcal{S}|}\sum_{s\in\mathcal{S}}m^{(d)}_s(h).
$$

An unweighted least-squares fit across the five lag values described the temporal profile as

$$
M_d(h)=a_d+\beta_d^{\mathrm{recovered}}(h-1).
$$

The fitted signed slope, in squared natural-log-frequency units per week, was the primary recovery descriptor. The displacement covariance was not divided by lag before fitting. A positive effect was described as resolved when the lower endpoint of its 95% subject-bootstrap interval exceeded zero.

Step 13 additionally exported a sensitivity estimate using bin weights $n^{(d)}_{shb}-1$ within subject while retaining equal subject weights. It also evaluated mean within-replicate displacement variance and replicate-specific excess on the same fixed cohort. Subject slopes, sign-flip tests, binwise slopes and finite-lag model descriptors were retained as ancillary outputs; the dose-response summaries selected `equal_bin_equal_subject` and `cross_cov` explicitly.

## 8. Joint subject bootstrap and incremental recovery

Uncertainty was estimated using 2,000 subject-bootstrap draws with seed `123`. Subjects were ordered by identifier. For each bootstrap draw, nine subject indices were sampled with replacement and the same draw was propagated across all lags, metrics and estimators. Scenario analyses used the same ordered subject set and seed. Their matching was checked by reconstructing the subject draws and the primary bootstrap curves and slopes from the subject-level outputs.

For draw $k$, with resampled subjects $s^{*(k)}_1,\ldots,s^{*(k)}_9$, the primary curve was

$$
M_d^{*(k)}(h)
=\frac{1}{9}\sum_{j=1}^{9}m^{(d)}_{s^{*(k)}_j}(h).
$$

The temporal slope was refitted to each resampled curve. The 2.5th and 97.5th percentiles supplied the interval endpoints. These resamples did not regenerate frequencies or counts and did not refit the state-inference model.

Absolute recovered slopes were retained for all four scenarios. Secondary recovery summaries compared each positive dose with the plateau:

$$
\Delta\beta_d^{\mathrm{recovered}}
=\beta_d^{\mathrm{recovered}}-\beta_0^{\mathrm{recovered}},
\qquad
\Delta\beta_d^{\mathrm{oracle}}
=\beta_d^{\mathrm{oracle}}-\beta_0^{\mathrm{oracle}}.
$$

Paired bootstrap differences were formed within matching draws,

$$
\Delta\beta_d^{*(k)}
=\beta_d^{\mathrm{recovered},*(k)}
-\beta_0^{\mathrm{recovered},*(k)}.
$$

Intervals for these contrasts were the corresponding percentiles of the paired difference distribution, not differences between marginal interval endpoints. The summarizer matches records by bootstrap ID; the separate verification step establishes that these IDs identify the same subject resamples across scenarios.

For nonzero oracle increments, incremental recovery was

$$
\rho_d
=\frac{\Delta\beta_d^{\mathrm{recovered}}}
{\Delta\beta_d^{\mathrm{oracle}}}.
$$

Its bootstrap distribution divided each paired recovered increment by the fixed oracle increment. A descriptive through-origin fit across positive doses gave

$$
\widehat{\lambda}
=\frac{\sum_{d>0}\Delta\beta_d^{\mathrm{oracle}}
\Delta\beta_d^{\mathrm{recovered}}}
{\sum_{d>0}\left(\Delta\beta_d^{\mathrm{oracle}}\right)^2}.
$$

An unconstrained descriptive linear fit of absolute recovered against oracle slopes across all four scenarios was also exported. Neither fit was used to recalibrate the doses or correct the empirical results. The oracle and recovered quantities have different support and weighting definitions; their ratio is therefore a conditional recovery benchmark, not a universal measure of estimator accuracy.

This experiment comprises one coupled realization per dose under the specified reference, random fields and count-generation procedure. Subject-bootstrap draws quantify resampling uncertainty within these realizations. They do not constitute independent simulations, establish a power curve or define a general detection limit for temporal accumulation. The validation does not assess recovery outside the inherited detection-support structure.

---

# Part II. Current GitHub computational workflow

## 9. Scope and repository integration

The current public workflow **does not rebuild the finalized empirical ClonoDynamics analysis from scratch inside the positive-control branch**. Instead, it consumes the already completed empirical outputs from the `core` and `validation` analysis blocks, freezes the required geometry and support contracts, generates the synthetic scenarios, and then reruns only the production components required for recovery.

This differs from the historical manuscript-reconstruction runbook, which contained a separate from-zero empirical reconstruction phase. The scientific positive-control model in Part I is unchanged; only the public orchestration layer has been simplified and aligned with the final repository.

The relevant repository layout is:

```text
code/
|-- orchestrator_clonodynamics.py
|-- 01_core/
|   |-- 2-multirepresentation_clonotype_state_inference.py
|   |-- noiseK_latent.py
|   |-- 4-multirepresentation_trajectory_assembly.py
|   `-- 5-longitudinal_transition_assembly.py
|
`-- 03_validation/
    |-- 11-cross_replicate_fluctuation_dynamics.py
    |-- 13-temporal_fluctuation_scaling.py
    `-- positive_controls/
        |-- 01_freeze_empirical_calibration.py
        |-- 02_calibrate_oracle_doses.py
        |-- 03_generate_support_conditioned_repertoires.py
        |-- 04_summarize_temporal_recovery.py
        |-- 05_plot_synthetic_validation.py
        `-- verify_workflow.py
```

The master launcher embeds the validated positive-control orchestration engine. Separate runtime `orchestrator_positive_controls.py` files are therefore not required in the final public layout.

## 10. Required inputs

The interactive positive-control workflow requests four inputs:

1. the corrected longitudinal repertoire dataset;
2. the completed `01_core` results root;
3. the completed `03_validation` results root (by default this may be the same shared longitudinal results root);
4. the destination directory for positive-control outputs.

With the standard ClonoDynamics layout, the destination defaults to:

```text
<validation-results-root>/positive_controls/
```

The empirical results must already contain the finalized outputs required by the freeze step, including:

```text
Step 1
    repertoire_characteristics.csv

Step 2
    latent_pair_qc.csv
    latent_params_by_pair.csv

Step 5
    longitudinal_transitions.parquet
    transition_assembly_report.md

Step 11
    00_run_config.json
    01_fluctuation_bin_edges.csv
    06_support_by_dt.csv

Step 13
    00_run_config.json
    02_complete_case_subjects.csv
    05_temporal_slope_summary.csv
```

The current positive-control workflow verifies source-code versions/hashes and input contracts before executing downstream phases.

## 11. Interactive execution

Launch ClonoDynamics from the repository root:

```bash
python3 code/orchestrator_clonodynamics.py
```

Select:

```text
Analysis
  -> positive_controls
```

The positive-control engine then requests the paths listed above and displays its execution plan.

Existing outputs are handled conservatively at each stage. A completed compatible stage may be skipped. Rewriting a stage deletes only the outputs owned by that stage and regenerates them. If an upstream calibration, oracle, generated scenario, or synthetic analysis changes during the current run, dependent downstream outputs are marked stale and cannot be silently reused.

## 12. Current ordered stages

The public workflow contains seven stages.

### 12.1 Stage 1 — freeze finalized empirical calibration

`code/03_validation/positive_controls/01_freeze_empirical_calibration.py` reads the corrected empirical repertoire tables together with the finalized empirical Step-1, Step-2, Step-5, Step-11 and Step-13 outputs.

It creates:

```text
positive_controls/calibration/
|-- 00_manifest.json
|-- 01_raw_sampling_design.csv
|-- 02_analysis_sampling_design.csv
|-- 03_step2_pair_parameters.csv
|-- 04_observation_structure_targets.json
|-- 05_cohort_and_analysis_summary.json
|-- 06_step11_bin_edges.csv
|-- 07_step11_support_by_dt.csv
|-- 08_step13_run_config.json
|-- 08_step13_core_selection.csv
|-- 09_step13_complete_case_subjects.csv
|-- 10_step13_temporal_slope_summary.csv
|-- 11_subject_union_references.csv
|-- 12_source_repertoire_audit.csv
|-- 13_validation_contract.json
`-- subject_union_repertoires/
```

For the reported study, the frozen cohort-specific guards are 10 subjects, 114 repertoires, 57 raw pairs, 56 analysis-eligible pairs, and 130 within-subject interval pairs. These values are reproducibility guards for the reported dataset, not generic ClonoDynamics requirements.

### 12.2 Stage 2 — oracle dose calibration

`code/03_validation/positive_controls/02_calibrate_oracle_doses.py` reads the frozen calibration and constructs the four reported scenarios:

```text
R0p00
R0p25
R0p50
R1p00
```

The finalized orchestration settings are:

```text
reference operational alpha = 0.05
sigma_fast                 = 0.20
oracle/generator seed      = 20260918
oracle targets R           = 0, 0.25, 0.50, 1.00
maximum calibration steps  = 22
target tolerance           = 0.0005
```

Outputs are written to:

```text
positive_controls/oracle/
|-- 00_qbio_oracle_calibration.json
|-- 01_scenario_table.csv
|-- 02_oracle_temporal_profiles.csv
|-- 03_search_trace.csv
`-- 04_manifest.json
```

Later stages read the calibrated values directly from `01_scenario_table.csv`; they do not re-enter rounded values from documentation.

### 12.3 Stage 3 — generate support-conditioned scenarios

`code/03_validation/positive_controls/03_generate_support_conditioned_repertoires.py` is called once for each oracle scenario, using the scenario-specific `sigma_fast`, `q_bio`, and seed recorded by Stage 2.

Each scenario is written under:

```text
positive_controls/scenarios/<scenario>/
|-- repertoires/
|-- ground_truth/
|-- 00_manifest.json
|-- 01_generated_sampling_summary.csv
|-- 02_oracle_pair_summary.csv
`-- 03_support_audit_by_subject.csv
```

For the reported support-conditioned design, each scenario contains 112 synthetic repertoire files corresponding to the 56 frozen analysis-eligible replicate pairs. The known trajectories in `ground_truth/` are never supplied to the production inference pipeline.

### 12.4 Stage 4 — production recovery pipeline

Each generated scenario is analyzed with the same production implementations used by ClonoDynamics:

```text
Step 2 -> Step 4 -> Step 5 -> Step 11 -> Step 13
```

The scenario-specific result tree is:

```text
positive_controls/scenarios/<scenario>/
|-- 2-clonotype_state_inference/
|-- 4-multirepresentation_trajectory_assembly/
|-- 5-longitudinal_transition_assembly/
|-- 11-cross_replicate_fluctuation_dynamics/
`-- 13-temporal_fluctuation_scaling/
```

The recovery analysis differs from an unconstrained new ClonoDynamics analysis only where the validation design explicitly freezes empirical geometry:

- Step 11 receives the empirical bin edges from `calibration/06_step11_bin_edges.csv`;
- Step 13 receives the empirical selected core-bin list from the frozen calibration;
- the operational threshold remains `alpha = 0.05`;
- Step 11/13 use 2,000 subject-bootstrap draws with seed `123`;
- the generated ground truth is not an analysis input.

Immediately after each synthetic Step 2, `verify_workflow.py fitted-pairs` checks agreement with the frozen eligible-pair design.

### 12.5 Stage 5 — strict workflow verification

`code/03_validation/positive_controls/verify_workflow.py completed` verifies the completed validation workspace before recovery summaries are accepted.

The verification layer checks, among other contracts:

- scenario parameters and seeds against the oracle table;
- fitted-pair identity against the frozen empirical design;
- transition/support consistency;
- Step-11 abundance edges;
- Step-13 core and complete-case cohort;
- Step-11/Step-13 signature linkage;
- reconstruction of the primary bootstrap subject draws, temporal curves, and slopes;
- matching bootstrap IDs across scenarios before paired differences are calculated.

Outputs are written to:

```text
positive_controls/verification/
|-- verification_checks.csv
|-- verification_summary.json
|-- bootstrap_subject_draws.csv
`-- verified_paired_increments.csv
```

Verification tests computational consistency. It does not require an estimated slope to have a desired sign.

### 12.6 Stage 6 — summarize temporal recovery

`code/03_validation/positive_controls/04_summarize_temporal_recovery.py` reads the oracle and scenario-specific Step-13 outputs and selects the primary `equal_bin_equal_subject/cross_cov` recovery estimates.

Outputs are:

```text
positive_controls/summary/
|-- 00_validation_summary.csv
|-- 01_temporal_profiles.csv
|-- 02_incremental_recovery_summary.csv
|-- 03_recovery_fit_summary.csv
`-- 04_validation_metadata.json
```

Absolute recovered slopes and baseline-adjusted increments remain separate. Paired contrast intervals use matched bootstrap records.

### 12.7 Stage 7 — publication plots and source data

`code/03_validation/positive_controls/05_plot_synthetic_validation.py` reads the numerical summary tables without rerunning simulation or inference.

The standard output directory is:

```text
positive_controls/figures/
```

and includes PDF figure panels together with:

```text
Table_validation_dose_response_publication.csv
SourceData_Figure_validation_A_temporal_profiles.csv
SourceData_Figure_validation_B_incremental_recovery.csv
```

The global ClonoDynamics figure workflow can later regenerate the positive-control publication figures from the finalized `summary/` directory without rerunning the positive-control analysis.

## 13. Current scientific settings reproduced by the orchestration engine

The finalized public positive-control workflow freezes the following reported-study settings:

```text
alpha                 = 0.05
sigma_fast            = 0.20
oracle/generator seed = 20260918
bootstrap seed        = 123
bootstrap replicates  = 2000
scenarios             = R0p00, R0p25, R0p50, R1p00
```

For every synthetic scenario, the production recovery chain uses:

```text
Step 2
    pair-specific latent inference from the synthetic technical counts

Step 4
    multirepresentation trajectory assembly

Step 5
    all finite-time transitions at lags 1-5

Step 11
    signed cross-replicate covariance on common4
    frozen empirical abundance-bin edges

Step 13
    frozen empirical abundance core
    scenario-verified complete-case cohort
    equal-bin/equal-subject primary estimator
    2000 joint subject-bootstrap draws
```

The hidden simulated trajectories remain isolated from this production chain.

## 14. Reproducibility and stopping rules

The master records the positive-control configuration, source-code inventory/hashes, effective commands, logs, and stage manifest under the selected positive-control workspace.

The workflow is fail-fast for attempted stages. Missing or incompatible required empirical inputs stop execution before downstream generation. Completed stage outputs can be reused only when their declared contracts are still current. When an upstream stage is regenerated, dependent downstream outputs are treated as stale and must be regenerated rather than silently mixed with the new upstream state.

This policy is intentionally stricter than treating the positive-control directory as a collection of independent files: the reported recovery summaries are valid only as a linked calibration → oracle → generation → recovery → verification chain.

## 15. Study-specific versus general quantities

Several numerical values in Part I are deliberately retained because this document also records the reported validation experiment. They include cohort size, repertoire count, eligible-pair count, selected abundance-core bins, complete-case subjects, calibrated `q_bio` values, and observed support totals.

These quantities should be read as **reported-study configuration**, not as software defaults that define ClonoDynamics generally. The general methodological contract is instead:

1. freeze an empirical design and analysis geometry;
2. calibrate known cumulative signals without using downstream recovery;
3. generate counts conditional on the frozen empirical support;
4. run the unchanged production inference/recovery chain;
5. verify support and bootstrap pairing;
6. compare recovered and oracle temporal effects without using the oracle to correct the empirical result.
