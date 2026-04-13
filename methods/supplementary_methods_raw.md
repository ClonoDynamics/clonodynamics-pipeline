---
title: "Supplementary Methods"
toc: true
toc-depth: 2
numbersections: true
---
\newpage

# Quality Control of RepSeq

This document describes in detail the logic, mathematical assumptions, command-line arguments, internal workflow, and outputs of the script `1-powerlaw_qc_analysis.py`. 

## General Purpose of the Script

The script performs **analysis only**. It does not generate figures. Its goal is to implement two conceptually distinct but related analyses on longitudinal immune repertoire replicate files:

1. **Step A1: empirical heavy-tail validation and tail-model comparison**
2. **Step A2: technical quality control based on replicate concordance of the power-law exponent gamma**

The script therefore serves as an upstream analytical stage whose outputs are tabular summaries that can later support figure generation, downstream filtering, or methodological justification.

In conceptual terms, the script answers two questions:

- Are the clone-size distributions in each repertoire consistent with a heavy-tailed structure, and is a power law favored over a truncated log-normal in the high-frequency tail?
- Are the two technical replicates of the same sample mutually concordant in their inferred tail exponent gamma?

These two tasks are deliberately separated:

- **A1** addresses the **biological-statistical structure** of the repertoire distribution.
- **A2** addresses **technical reproducibility** across replicate measurements.

---

## Expected Input Structure

The script expects a directory containing repertoire files named according to a regular-expression pattern. By default, the filename pattern is:

```text
^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])$
```

This means that filenames are expected to encode:

- `subject`: subject identifier
- `time`: timepoint identifier
- `replica`: replicate number, either 1 or 2

Examples of valid filenames under the default pattern are:

```text
1_1-1
1_1-2
3_4-1
3_4-2
```

If the files carry an extension such as `.tsv`, the filename pattern must be adjusted at launch using the `--pattern` argument.

Each repertoire file must contain at least the columns:

- `aaSeqCDR3`
- `readCount`

Only these two columns are used by the script. Any additional columns are ignored. The script reads the files using `pandas.read_csv` with a separator defined by `--file_sep`.

---

## High-Level Workflow

The script proceeds through the following stages:

1. Parse command-line arguments.
2. Index all files in the input directory and identify complete replicate pairs.
3. For each replicate pair:
   - read both files
   - compute repertoire depth
   - convert counts to frequencies
   - fit a power law separately to each replicate
   - store replicate-level fit statistics
   - compute replicate-pair discordance metrics
   - perform threshold-based tail-model comparison for each replicate
4. Aggregate replicate-level and pair-level results.
5. Apply gamma-only quality control to each replicate pair.
6. Write analysis tables to disk.

The script never edits input files and never produces plots.

---

## Utility Functions

### `ensure_dir`

This function creates the output directory if it does not already exist.

If `p` is a path, then:

$$
\text{ensure\_dir}(p) \to p
$$

with side effect:

$$
\text{mkdir}(p, \text{parents}=\text{True}, \text{exist\_ok}=\text{True})
$$

Its role is purely infrastructural.

### `robust_mad`

This function computes the **median absolute deviation** (MAD), a robust estimator of scale.

Given a vector $x = (x_1, \dots, x_n)$, define the median:

$$
\tilde{x} = \operatorname{median}(x)
$$

Then the MAD is:

$$
\operatorname{MAD}(x) = \operatorname{median}(|x_i - \tilde{x}|)
$$

Only finite values are retained. If no finite values are available, the function returns `NaN`.

This quantity is used later to define a robust outlier score for replicate discordance in gamma.

### `robust_z`

This function defines a robust z-score:

$$
z_i = \frac{x_i - \operatorname{median}(x)}{\operatorname{MAD}(x)}
$$

If the MAD is not finite or is numerically too small, the function returns `NaN` for all entries.

This robust standardization is used because gamma discordance may not be well modeled by a Gaussian distribution and may contain outliers.

### `depth_ratio`

Given two replicate sequencing depths $d_1$ and $d_2$, the function returns:

$$
\text{depth\_ratio}(d_1, d_2) = \frac{\max(d_1, d_2)}{\min(d_1, d_2)}
$$

provided that both are strictly positive. Otherwise, it returns `NaN`.

This metric is only diagnostic in this script. It does **not** enter the final quality-control rule.

---

## Mathematical Basis of the Power-Law Fit

The script uses a **continuous Clauset-style power-law fitting procedure** on clonotype frequencies. Frequencies are treated as positive real values.

### Input Frequencies

For each repertoire, the script converts raw counts into relative frequencies:

$$
f_i = \frac{c_i}{\sum_j c_j}
$$

where:

- $c_i$ is the read count of clonotype $i$
- $\sum_j c_j$ is the total repertoire depth

Only strictly positive finite frequencies are kept.

### Validation of Positive Values

The function `validate_positive_1d`:

- converts input to floating point
- removes non-finite values
- removes non-positive values
- sorts the remaining values in increasing order

Thus the analysis domain is:

$$
x_i > 0
$$

for all retained clonotype frequencies.

### Continuous Pareto Model

The fitted power-law density on the tail $x \ge x_{\min}$ is:

$$
p(x \mid \alpha, x_{\min}) = (\alpha - 1) x_{\min}^{\alpha - 1} x^{-\alpha}
$$

valid for:

$$
\alpha > 1, \qquad x_{\min} > 0
$$

The script denotes the fitted exponent as both `alpha_hat` and `gamma_hat`. In practice, for the purpose of this workflow, gamma is simply the fitted power-law exponent.

### Maximum-Likelihood Estimator of the Exponent

Given the tail sample $x_1, \dots, x_n$ with all $x_i \ge x_{\min}$, the continuous maximum-likelihood estimator is:

$$
\hat{\alpha} = 1 + \frac{n}{\sum_{i=1}^{n} \log(x_i / x_{\min})}
$$

This is implemented in `continuous_alpha_mle`.

If the denominator is non-positive, the estimate is considered invalid and the function returns `NaN`.

### Model CDF

The theoretical cumulative distribution function of the continuous Pareto law is:

$$
F(x) = 1 - \left(\frac{x}{x_{\min}}\right)^{1 - \alpha}
$$

for $x \ge x_{\min}$.

This is implemented by `powerlaw_cdf_continuous`.

### Kolmogorov-Smirnov Distance

For a sorted tail sample $x_{(1)} \le \dots \le x_{(n)}$, the empirical CDF is:

$$
F_{\text{emp}}(x_{(k)}) = \frac{k}{n}
$$

The model CDF is $F_{\text{model}}(x_{(k)})$. The KS distance is:

$$
D_{\text{KS}} = \max_k \left| F_{\text{emp}}(x_{(k)}) - F_{\text{model}}(x_{(k)}) \right|
$$

This statistic is used to choose the best tail threshold $x_{\min}$.

### KS-Based Estimation of the Tail Threshold

The function `estimate_xmin_ks_continuous` evaluates candidate values of $x_{\min}$ drawn from the sorted unique observed frequencies.

For each candidate threshold:

1. Keep the tail sample $x_i \ge x_{\min}$.
2. Require at least `min_tail_size` observations.
3. Estimate $\hat{\alpha}$ by MLE.
4. Compute the KS distance.
5. Retain the threshold that minimizes the KS distance.

Formally, the selected threshold is:

$$
\hat{x}_{\min} = \arg\min_{x_{\min}} D_{\text{KS}}(x_{\min})
$$

subject to the tail-size constraint.

The output consists of:

- `xmin_hat`
- `alpha_hat`
- `ks_distance`
- `n_tail`

### Final Per-Repertoire Power-Law Fit

The function `fit_powerlaw_on_frequencies` wraps the above machinery and returns a `PowerLawFitResult` object containing:

- fit success or failure
- diagnostic message
- number of total positive observations
- number of tail observations
- estimated exponent
- estimated threshold
- KS distance
- tail fraction

The tail fraction is:

$$
\text{tail\_fraction} = \frac{n_{\text{tail}}}{n_{\text{total}}}
$$

---

## Tail-Model Comparison: Power Law Versus Truncated Log-Normal

The second major component of Step A1 is a model comparison restricted to the high-frequency tail.

Importantly, the script does **not** compare models using model-specific thresholds. Instead, it defines tail regions by **frequency quantiles**. This is methodologically important because it enforces comparison on the **same support** for both candidate models.

### Threshold Quantiles

Let the vector of positive frequencies in a repertoire be $x_1, \dots, x_n$. For each quantile level $q$ in the user-provided list, the script defines:

$$
x_{\min}^{(q)} = Q_q(x)
$$

where $Q_q(x)$ is the empirical $q$-quantile.

The tail is then:

$$
\{x_i : x_i \ge x_{\min}^{(q)}\}
$$

A tail is analyzed only if it contains at least `scan_min_tail_size` observations.

### Power-Law Log-Likelihood on the Tail

For a tail sample $x_1, \dots, x_n$ and fitted exponent $\alpha$, the continuous Pareto log-likelihood is:

$$
\ell_{\text{PL}} = n \log(\alpha - 1) + n(\alpha - 1) \log(x_{\min}) - \alpha \sum_{i=1}^{n} \log(x_i)
$$

This is implemented in `powerlaw_tail_loglik`.

### Truncated Log-Normal Model

The competing model is a log-normal distribution truncated below at $x_{\min}$.

Define:

$$
z_i = \log x_i
$$

The plug-in parameter estimates are:

$$
\hat{\mu} = \frac{1}{n} \sum_{i=1}^{n} z_i
$$

and

$$
\hat{\sigma} = \sqrt{\frac{1}{n-1} \sum_{i=1}^{n} (z_i - \hat{\mu})^2}
$$

The untruncated normal density in log-space is combined with:

- the Jacobian factor $1/x$
- the survival-function normalization induced by truncation at $x_{\min}$

The resulting truncated log-normal log-likelihood is:

$$
\ell_{\text{TLN}} = \sum_{i=1}^{n} \left[ \log \phi(z_i \mid \mu, \sigma) - \log x_i \right] - n \log S(x_{\min})
$$

where:

- $\phi$ is the Gaussian density in log-space
- $S(x_{\min}) = 1 - \Phi\left( \frac{\log x_{\min} - \mu}{\sigma} \right)$ is the survival probability above the threshold

This is implemented in `truncated_lognormal_tail_loglik`.

### Delta Log-Likelihood

The model comparison statistic is:

$$
\Delta \ell = \ell_{\text{PL}} - \ell_{\text{TLN}}
$$

Interpretation:

- $\Delta \ell > 0$: the power law is preferred
- $\Delta \ell < 0$: the truncated log-normal is preferred

The script also defines a per-clonotype average:

$$
\overline{\Delta \ell} = \frac{\Delta \ell}{n_{\text{tail}}}
$$

This normalization is useful when tail sizes differ across thresholds or repertoires.

### Threshold Scan Output

For each repertoire and each threshold quantile, the script records:

- threshold quantile
- threshold value `xmin_threshold`
- total number of positive clonotypes
- tail size and tail fraction
- power-law parameters and KS statistic
- truncated log-normal parameters
- total and mean delta log-likelihood
- Boolean indicator `pl_better`

This is the raw material for evaluating whether model preference is stable across increasingly stringent high-frequency tails.

---

## File Reading, Pair Indexing, and Replicate Processing

### `read_rep_file`

Each file is read and restricted to the columns `aaSeqCDR3` and `readCount`.

The function:

- coerces `readCount` to numeric
- removes non-finite entries
- removes non-positive entries

Thus, every retained clonotype satisfies:

$$
\text{readCount}_i > 0
$$

### `index_pairs`

This function scans the input directory, matches filenames to the user-provided regex pattern, and stores files indexed by `(subject, time, replica)`.

A pair is considered complete only if **both** replicate 1 and replicate 2 exist for the same `(subject, time)`.

The output is a list of tuples:

```text
(subject, time, file_rep1, file_rep2)
```

### `fit_one_replicate`

This function encapsulates the per-replicate analysis:

1. Read the repertoire file.
2. Compute total depth:

$$
\text{depth} = \sum_i \text{readCount}_i
$$

3. Convert counts to frequencies:

$$
f_i = \frac{\text{readCount}_i}{\text{depth}}
$$

4. Fit the continuous power law on frequencies.

The function returns:

- the fitted `PowerLawFitResult`
- the total depth
- the vector of frequencies

---

## Gamma-Only Quality Control

The second major output of the script is a QC classification for each replicate pair based **only** on gamma concordance.

This is intentionally restrictive. Other quantities such as `xmin`, KS distance, depth ratio, tail size, or tail fraction are retained as diagnostics but do not affect PASS or FAIL.

### Pairwise Discordance in Gamma

For each sample pair, the script extracts:

- $\gamma_1$: fitted exponent in replicate 1
- $\gamma_2$: fitted exponent in replicate 2

Then it computes:

$$
\Delta \gamma = |\gamma_1 - \gamma_2|
$$

This is the central QC variable.

### Pairwise Mean Gamma

A descriptive quantity also recorded is:

$$
\gamma_{\text{mean}} = \operatorname{mean}(\gamma_1, \gamma_2)
$$

computed as a NaN-aware mean.

### Robust Standardization of Gamma Discordance

Across all replicate pairs, the script computes a robust z-score for $\Delta \gamma$:

$$
z_{\Delta \gamma} = \frac{\Delta \gamma - \operatorname{median}(\Delta \gamma)}{\operatorname{MAD}(\Delta \gamma)}
$$

This score is not a classical Gaussian z-score. It is a median-and-MAD standardized outlier score.

### QC Threshold

The user provides a cutoff `gamma_mad_cutoff`, with default:

$$
\text{gamma\_mad\_cutoff} = 3.5
$$

A pair is marked as a gamma outlier if:

$$
|z_{\Delta \gamma}| > \text{gamma\_mad\_cutoff}
$$

### PASS/FAIL Rule

The QC rule is intentionally simple:

A pair is **FAIL** if either of the following conditions holds:

1. one replicate fit failed
2. gamma discordance is an outlier

Formally:

$$
\text{FAIL if } (\neg \text{fit\_success\_rep1}) \lor (\neg \text{fit\_success\_rep2}) \lor \left( |z_{\Delta \gamma}| > c \right)
$$

where $c$ is `gamma_mad_cutoff`.

Otherwise:

$$
\text{PASS}
$$

There is **no WARN class** in this version of the script.

### QC Reason Field

For each failing pair, the script composes a reason string such as:

- `fail:fit_failed`
- `fail:delta_gamma`
- `fail:fit_failed,delta_gamma`

This makes the failure mechanism explicit in the output table.

---

## Aggregate Summaries

After computing per-replicate and per-pair results, the script generates summary tables.

### Overall QC Summary

The overall summary includes:

- number of pairs
- number of PASS
- number of FAIL
- median `delta_gamma`
- median `delta_xmin`
- median `depth_ratio`
- median `gamma_mean`

These are descriptive summaries of the full dataset.

### Subject-Level QC Summary

The script also groups by subject and computes:

- number of pairs
- number of PASS
- number of FAIL
- median `delta_gamma`
- median `delta_xmin`
- median `depth_ratio`
- median `gamma_mean`

This is useful for identifying donors with consistently unstable replicate concordance.

---

## Detailed Explanation of Command-Line Arguments

The script exposes the following CLI arguments.

### Input and Output Arguments

#### `--data_dir`

Path to the directory containing replicate files.

This argument is required.

#### `--out_dir`

Path to the output directory where all analysis tables are written.

This argument is required.

#### `--file_sep`

Field separator used when reading repertoire files.

Default:

```text
\t
```

This corresponds to tab-separated files.

#### `--pattern`

Regular-expression pattern used to parse filenames.

Default:

```text
^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])$
```

This pattern is critical because it determines how files are grouped into replicate pairs.

### Power-Law Fit Arguments

#### `--min_tail_size`

Minimum number of observations required in the tail during KS-based power-law fitting.

Default:

```text
50
```

If a repertoire contains fewer than `min_tail_size` positive frequencies, the fit fails.

If a candidate threshold leaves fewer than `min_tail_size` tail observations, that threshold is ignored.

#### `--candidate_step`

Step size used to thin candidate $x_{\min}$ values.

Default:

```text
1
```

If set to 1, every unique observed value is considered as a candidate threshold. Larger values subsample the candidate set and accelerate computation, at the possible cost of precision in threshold estimation.

### Threshold-Scan Arguments

#### `--threshold_quantiles`

Comma-separated list of quantile levels used to define tail thresholds for model comparison.

Default:

```text
0.70,0.75,0.80,0.85,0.90,0.95
```

Each value must satisfy:

$$
0 < q < 1
$$

These quantiles determine increasingly strict definitions of the high-frequency tail.

#### `--scan_min_tail_size`

Minimum number of observations required for the threshold-scan model comparison.

Default:

```text
50
```

This is conceptually separate from `--min_tail_size`, although they may be set equal in practice.

### QC Argument

#### `--gamma_mad_cutoff`

Robust outlier cutoff applied to $z_{\Delta \gamma}$.

Default:

```text
3.5
```

This controls the stringency of replicate concordance QC.

Larger values make the QC more permissive; smaller values make it stricter.

---

## Detailed Description of Output Files

The script writes nine output tables.

### `A1_repertoire_tail_fits.csv`

This file contains one row per replicate and summarizes the KS-based power-law fit.

It includes:

- subject, time, replicate, pair ID
- file name
- repertoire depth
- fit success and message
- total positive clonotypes
- tail size
- tail fraction
- fitted gamma
- fitted xmin
- KS distance

This table answers: how well can each replicate be described by a fitted power-law tail?

### `A1_repertoire_threshold_scan.csv`

This file contains one row per replicate per threshold quantile.

It includes the full tail-model comparison on quantile-defined tails:

- threshold quantile
- threshold value
- tail size and tail fraction
- power-law exponent and log-likelihood
- truncated log-normal log-likelihood and parameters
- total and mean delta log-likelihood
- indicator of whether power law is preferred

This is the central table for evaluating robustness of tail-model preference.

### `A1_repertoire_model_compare.csv`

This file aggregates the threshold scan at the level of each replicate.

It summarizes, for each repertoire:

- number of threshold levels analyzed
- mean and median normalized delta log-likelihood
- fraction of thresholds where power law is preferred
- median tail fraction

This answers whether model preference is stable across thresholds within a repertoire.

### `A1_pair_tail_preference_summary.csv`

This file aggregates the threshold-scan comparison at the replicate-pair level.

It summarizes:

- mean of the replicate-level median delta log-likelihood
- mean fraction of thresholds favoring power law
- minimum fraction of thresholds favoring power law
- number of replicates contributing to the pair

This table is useful when model preference must be summarized at the sample level rather than at the individual replicate level.

### `A2_replicate_powerlaw_fits.csv`

This file stores the raw replicate-level power-law fit outputs used by the QC step.

It is essentially the pre-QC per-replicate fit table.

### `A2_pair_gamma_qc.csv`

This is the main QC output.

Each row corresponds to one replicate pair and includes:

- replicate filenames
- fit success flags
- depths and depth ratio
- gamma values and delta gamma
- xmin values and delta xmin
- KS distances
- tail sizes and tail fractions
- robust z-score of delta gamma
- outlier flag
- PASS/FAIL label
- QC reason

This is the key file that can later be used to retain or exclude samples in downstream analysis.

### `A2_gamma_qc_thresholds.csv`

This file records the QC rule and the numerical threshold used.

It documents the cutoff and the precise PASS/FAIL logic, improving reproducibility.

### `A2_summary_overall.csv`

This is a single-row summary of the full dataset after QC classification.

### `A2_summary_by_subject.csv`

This file provides the same type of descriptive summary stratified by subject.

---

## Conceptual Interpretation of Step A1

Step A1 is not merely a fitting exercise. It is a methodological justification stage.

Its purpose is to determine whether the high-frequency clone-size distribution is compatible with a heavy tail and whether a power-law tail provides a better account than a truncated log-normal on matched support regions.

Two design choices are especially important.

### Quantile-Based Tail Thresholds

By scanning quantiles rather than fitting separate model-specific thresholds, the script ensures that both candidate models are evaluated on the same tail subset. This avoids a common methodological ambiguity in tail-model comparison.

### Per-Clonotype Normalized Delta Log-Likelihood

The quantity:

$$
\overline{\Delta \ell} = \frac{\ell_{\text{PL}} - \ell_{\text{TLN}}}{n_{\text{tail}}}
$$

allows comparison across repertoires and thresholds even when the number of tail clonotypes differs.

This is particularly important when the tail becomes progressively smaller at high quantiles.

---

## Conceptual Interpretation of Step A2

Step A2 uses replicate concordance in gamma as a technical reproducibility metric.

The rationale is that if two technical replicates are measuring the same underlying repertoire, then the inferred power-law exponent of the high-frequency tail should be relatively stable across replicates. Large discrepancies in gamma may indicate technical instability, insufficient sampling, fitting failure, or other sources of unreliability.

The use of a robust z-score rather than a standard z-score makes the QC rule resistant to skewness and extreme values in the empirical distribution of $\Delta \gamma$.

At the same time, the choice to ignore other diagnostics in the final PASS/FAIL rule makes the procedure deliberately transparent and conservative in scope.

---

## What the Script Does Not Do

It is important to state explicitly what is outside the scope of this script.

The script does **not**:

- generate plots
- fit discrete power laws
- compare against exponential tails
- use `xmin`, KS, depth ratio, or tail fraction in the PASS/FAIL decision
- create a WARN category
- denoise repertoires
- merge replicates
- construct trajectories
- perform downstream dynamical inference

It is therefore an **analysis-and-diagnostics stage only**.

---

## Computational Logic of the Main Loop

Inside the main loop over replicate pairs, the script performs the following operations in sequence.

For each pair `(subject, time)`:

1. Fit replicate 1.
2. Fit replicate 2.
3. Store replicate-level fit outputs in `rep_rows`.
4. Compute pair-level diagnostics and store them in `pair_rows`.
5. For each replicate separately, run the threshold scan if enough positive frequencies are available.
6. Append threshold-scan tables to a list.

After the loop:

1. Convert accumulated rows into DataFrames.
2. Apply gamma-only QC.
3. Compute overall and subject-level summaries.
4. Aggregate threshold-scan tables at repertoire and pair level.
5. Save all outputs.

This structure is linear, explicit, and reproducible.

---

## Possible Failure Modes

Several types of failure are handled explicitly.

### Missing Required Columns

If a repertoire file lacks either `aaSeqCDR3` or `readCount`, the script raises an error.

### Zero or Non-Positive Depth

If the sum of counts is zero, the replicate fit fails.

### Too Few Positive Frequencies

If the number of positive frequencies is smaller than `min_tail_size`, the power-law fit fails.

### No Valid KS-Based Threshold

If none of the candidate thresholds yields a valid power-law tail with sufficient size and valid exponent, the fit fails.

### Too Few Tail Observations in Threshold Scan

If a quantile-defined tail contains fewer than `scan_min_tail_size` observations, that threshold is skipped.

These failure modes are propagated transparently into the output tables.

---

## Reproducibility and Transparency

A strength of the script is that it records both:

- the raw fit outputs
- the exact QC decision rule

This makes it possible to audit all decisions post hoc.

In particular, because the QC rule depends only on fit success and robust gamma discordance, one can later assess the influence of excluded samples without ambiguity.

---

## Recommended Interpretation in a Methods Section

In manuscript terms, the script supports the following methodological narrative.

First, clone-size distributions are examined for heavy-tailed structure and model preference in the high-frequency regime, using quantile-defined tails and matched-support comparison between power-law and truncated log-normal models.

Second, technical reproducibility is assessed by fitting the power-law exponent separately in each technical replicate and quantifying replicate discordance in gamma. Samples are classified as PASS or FAIL using a robust outlier rule based exclusively on delta gamma, while additional fit diagnostics are retained for descriptive purposes.

This separation between structural validation and technical QC is a central conceptual feature of the script.

\newpage

# Noise Model Pipeline for RepSeq Data


## Overview

The script `noise_pipeline.py` is an adaptation of the NoisET Framework of M. B. Koraichi, et al. J. Phys. Chem. A 126, 7407–7414 (2022).It implements a streamlined end-to-end pipeline for technical noise modelling, empirical denoising, and generation of trajectory-ready observability outputs from longitudinal immune repertoire sequencing data.

This pipeline is an adaptation of the NoisET framework originally developed by Mora and collaborators, reformulated here for a longitudinal setting with explicit technical replicates, external quality-control filtering, empirical null pooling, and downstream preparation for clonotype trajectory analysis.

In its current form, the pipeline performs the following operations:

1. selection of replicate pairs passing an external QC step  
2. pairwise fitting of a Negative Binomial noise model with power-law prior  
3. extraction of per-clonotype log-probabilities under the fitted null model  
4. construction of pooled null distributions at clone, subject, or global level  
5. empirical denoising through tail probabilities derived from those pooled nulls  
6. generation of per-timepoint observability calls and trajectory-ready summaries  
7. optional sensitivity analyses over different significance thresholds  

The pipeline therefore provides the bridge between raw replicate counts and a probabilistic definition of clonotype observability suitable for downstream dynamical inference.

---

## General Statistical Rationale

Immune repertoire sequencing data are affected by multiple sources of variability, including sampling noise, library preparation effects, amplification noise, and sequencing stochasticity. As a consequence, the observed count of a clonotype is not a direct measurement of its latent biological abundance, but a noisy realization of an underlying frequency.

The central idea of the pipeline is that technical replicate pairs contain direct information about this experimental noise. By fitting a probabilistic model to replicate counts, one can estimate the expected variability under a null hypothesis of purely technical fluctuation. This null model can then be used to assign to each clonotype a probability of being compatible with technical noise, and therefore to distinguish detectable biological signal from stochastic background.

This logic follows the core philosophy of NoisET: infer a probabilistic noise model from technical replicates, then use that model to calibrate the detectability of clonotypes.

---

## Inputs

The pipeline expects a directory containing technical replicate files named in the form

`<subject>_<time>-1`  
`<subject>_<time>-2`

or any equivalent pattern specified through a regular expression.

Each file must contain at least the columns:

- `aaSeqCDR3`
- `readCount`

An external QC table is also required. This table must contain at least:

- `subject`
- `time`
- `QC_pass`

Only replicate pairs with `QC_pass == True` are processed.

---

## External QC Filtering

This version of the pipeline does not perform QC internally. Instead, QC is assumed to have been computed upstream, and the pipeline uses an external table to decide which replicate pairs are admissible for modelling.

Formally, let

$$
\mathcal{P} = \{(s,t)\}
$$

be the set of all subject-time pairs present in the repertoire directory, and let

$$
\mathcal{P}_{\mathrm{QC}} = \{(s,t) \in \mathcal{P} : QC\_pass(s,t) = \mathrm{True}\}
$$

be the subset passing QC. The pipeline restricts all subsequent inference to the set

$$
\mathcal{P}_{\mathrm{QC}}.
$$

Pairs failing QC are excluded from:

- noise-model fitting
- null-pool construction
- denoising
- trajectory-ready outputs
- parameter aggregation

This design ensures that the statistical model is fitted only to replicate pairs judged reliable by the external QC procedure.

---

## Pairwise Data Representation

For each QC-pass pair, the two replicate files are merged by clonotype identity, using `aaSeqCDR3` as matching key. Missing counts are imputed as zero.

If the two replicate counts for clonotype $i$ are denoted by

$$
c_{i1}, \quad c_{i2},
$$

and the total sequencing depths of the two replicates are

$$
N_1, \quad N_2,
$$

then the corresponding observed frequencies are

$$
f_{i1} = \frac{c_{i1}}{N_1}, \qquad f_{i2} = \frac{c_{i2}}{N_2}.
$$

The pipeline also computes an auxiliary geometric-mean frequency:

$$
f_{i,\mathrm{geo}} =
\begin{cases}
\sqrt{f_{i1}f_{i2}} & \text{if } c_{i1} > 0 \text{ and } c_{i2} > 0, \\
f_{i1} & \text{if } c_{i1} > 0 \text{ and } c_{i2} = 0, \\
f_{i2} & \text{if } c_{i1} = 0 \text{ and } c_{i2} > 0, \\
0 & \text{otherwise.}
\end{cases}
$$

This quantity is not the inferential core of the model, but provides a convenient summary of clonotype abundance across the replicate pair.

---

## Core Noise Model: `noiseK_nb.py`

The inferential core of the pipeline is the function

`fit_noiseK_nb_powerlaw(...)`

implemented in `noiseK_nb.py`.

This function fits a probabilistic model in which:

1. each clonotype has an unknown latent frequency $f$  
2. latent frequencies are drawn from a power-law prior  
3. observed replicate counts arise from a Negative Binomial observation process conditional on $f$  

This is the essential NoisET-style component that transforms replicate counts into a likelihood-based model of technical noise.

---

## Power-Law Prior on Latent Frequencies

The latent clonotype frequency $f$ is assumed to lie in the interval

$$
f \in [f_{\min}, 1].
$$

The lower bound is defined from the mean sequencing depth of the replicate pair:

$$
N_{\mathrm{total}} = \frac{N_1 + N_2}{2},
$$

and

$$
f_{\min} = \frac{1}{N_{\mathrm{total}}}.
$$

The latent-frequency domain is discretized on a logarithmic grid:

$$
f^{(1)}, f^{(2)}, \dots, f^{(G)},
$$

where $G$ is the user-defined `grid_size`.

The prior over latent frequencies is a discrete approximation of a power-law distribution:

$$
\pi(f) \propto f^{-\gamma},
$$

where $\gamma > 0$ is an inferential parameter.

On the discretized grid, the log-prior is implemented as

$$
\log \pi_g = -\gamma \log f^{(g)} - \log Z,
$$

where $Z$ is the normalizing constant computed numerically through log-sum-exp stabilization.

This choice encodes the biologically motivated assumption that clonotype frequency distributions are heavy-tailed, with many rare clonotypes and comparatively few expanded ones.

---

## Negative Binomial Observation Model

Conditional on latent frequency $f$, the expected count in replicate $r$ is

$$
\mu_r = f N_r.
$$

Observed counts are modelled using a Negative Binomial distribution with mean $\mu_r$ and dispersion parameter $\kappa$ (called `k` in the implementation).

For a count $k_{\mathrm{obs}}$, the log-probability mass function is

$$
\log P\bigl(k_{\mathrm{obs}} \mid \mu, \kappa \bigr)
=
\log \Gamma(k_{\mathrm{obs}}+\kappa)
-
\log \Gamma(\kappa)
-
\log \Gamma(k_{\mathrm{obs}}+1)
+
\kappa \log \left( \frac{\kappa}{\kappa+\mu} \right)
+
k_{\mathrm{obs}} \log \left( \frac{\mu}{\kappa+\mu} \right).
$$

Equivalently, this corresponds to a Negative Binomial law parameterized by mean and inverse-noise scale, with variance

$$
\mathrm{Var}(K \mid f)
=
\mu + \frac{\mu^2}{\kappa}.
$$

Thus, when $\kappa$ is large, the model approaches Poisson-like noise, whereas smaller $\kappa$ implies stronger overdispersion.

For each clonotype $i$, the two replicate counts are treated as conditionally independent given the latent frequency:

$$
P(c_{i1}, c_{i2} \mid f)
=
P(c_{i1} \mid f, N_1, \kappa)\;
P(c_{i2} \mid f, N_2, \kappa).
$$

---

## Marginal Likelihood of a Clonotype

Because the latent frequency $f$ is not directly observed, it is integrated out with respect to the power-law prior.

For clonotype $i$, the marginal probability is

$$
P(c_{i1}, c_{i2} \mid \gamma, \kappa)
=
\int_{f_{\min}}^1
P(c_{i1} \mid f, N_1, \kappa)\;
P(c_{i2} \mid f, N_2, \kappa)\;
\pi(f \mid \gamma)\;
df.
$$

In practice, this integral is approximated on the logarithmic grid:

$$
P(c_{i1}, c_{i2} \mid \gamma, \kappa)
\approx
\sum_{g=1}^{G}
P(c_{i1} \mid f^{(g)}, N_1, \kappa)\;
P(c_{i2} \mid f^{(g)}, N_2, \kappa)\;
\pi_g.
$$

The implementation works in log-space. For each clonotype and each grid point, it computes

$$
\log P_{ig}^{(1)} = \log P(c_{i1} \mid f^{(g)}, N_1, \kappa),
$$

$$
\log P_{ig}^{(2)} = \log P(c_{i2} \mid f^{(g)}, N_2, \kappa),
$$

then forms the joint log-density

$$
\log J_{ig}
=
\log P_{ig}^{(1)}
+
\log P_{ig}^{(2)}
+
\log \pi_g.
$$

The clonotype-specific marginal log-probability is then

$$
\log P_i
=
\log \sum_{g=1}^{G} \exp(\log J_{ig}),
$$

computed numerically via `logsumexp`.

This quantity is precisely what the pipeline stores as the per-clonotype `logP`.

---

## Pairwise Log-Likelihood and Parameter Estimation

For a replicate pair containing $M$ clonotypes, the total log-likelihood is

$$
\mathcal{L}(\gamma, \kappa)
=
\sum_{i=1}^{M}
\log P(c_{i1}, c_{i2} \mid \gamma, \kappa).
$$

The function `fit_noiseK_nb_powerlaw(...)` estimates $\gamma$ and $\kappa$ by numerical maximization of this likelihood, or equivalently by minimization of the negative log-likelihood:

$$
(\hat{\gamma}, \hat{\kappa})
=
\arg \min_{\gamma > 0, \kappa > 0}
\left[
-\mathcal{L}(\gamma,\kappa)
\right].
$$

Optimization is performed with the L-BFGS-B algorithm under positivity constraints:

$$
\gamma > 0, \qquad \kappa > 0.
$$

The user provides initial values through:

- `gamma_init`
- `k_init`

and the grid resolution through:

- `grid_size`

The fitted model returns:

- success flag
- optimizer message
- maximized log-likelihood
- estimated parameters
- per-clonotype marginal log-probabilities

The parameter dictionary contains at least:

- `gamma`
- `k`
- `fmin`
- `N_total_mean`
- `N_total_min`
- `N_total_max`

---

## Interpretation of the Fitted Parameters

The parameter $\gamma$ controls the heaviness of the prior tail. Larger $\gamma$ implies stronger concentration toward low-frequency clonotypes, whereas smaller $\gamma$ corresponds to a flatter tail with relatively more mass in expanded clonotypes.

The parameter $\kappa$ controls experimental overdispersion. In the variance relation

$$
\mathrm{Var}(K \mid f)
=
\mu + \frac{\mu^2}{\kappa},
$$

smaller $\kappa$ implies larger noise beyond Poisson sampling, reflecting stronger technical variability.

The fitted pairwise log-likelihood summarizes how well the model explains the joint replicate counts under the inferred parameters.

The per-clonotype `logP` values quantify how compatible each clonotype is with the fitted technical-noise model.

---

## Per-Clonotype Log-Probability Tables

After model fitting, the pipeline stores for each QC-pass pair a per-clonotype table containing at least:

- clonotype identity
- replicate counts
- replicate frequencies
- geometric-mean frequency
- pair metadata
- per-clonotype `logP`

These files are written under the `per_clone_long_logP` directory and serve as the basis for all downstream empirical denoising steps.

---

## Null Model Choices in the Pipeline

Once the pairwise `logP` values have been computed, the pipeline constructs an empirical null distribution by pooling them according to one of three choices.

### Clone-level null

Each pair defines its own pool:

$$
\mathrm{pool\_id} = (subject, time).
$$

This is the most local null and preserves pair-specific noise structure.

### Subject-level null

All QC-pass pairs from the same subject are pooled:

$$
\mathrm{pool\_id} = subject.
$$

However, this is used only if the subject has at least a specified number of QC-pass pairs. If not, the pipeline falls back to the global pool.

### Global null

All QC-pass pairs across all subjects and times are pooled together:

$$
\mathrm{pool\_id} = \mathrm{GLOBAL}.
$$

This provides the largest empirical null but may average over subject-specific noise features.

---

## Quantile-Based Observability Threshold

For each pool, the pipeline computes a quantile cutoff of the empirical `logP` distribution.

If the chosen significance level is $\alpha$, then the cutoff is

$$
\log P_{\mathrm{cutoff}} = Q_{\alpha}(\log P),
$$

where $Q_{\alpha}$ denotes the empirical $\alpha$-quantile.

A clonotype-specific cutoff is assigned according to the pool to which that clonotype belongs. If a subject-level cutoff is missing, the pipeline uses the global cutoff as fallback.

This step yields an initial cutoff-based observability indicator:

$$
\mathrm{observable\_cutoff} =
\mathbf{1}\{\log P_i \ge \log P_{\mathrm{cutoff}}\}.
$$

This indicator is kept as an internal diagnostic, but the final observability status is defined later through empirical tail probabilities.

---

## Empirical Denoising Through Tail Probabilities

The pipeline then computes, for each clonotype, its empirical cumulative probability within the chosen null pool.

Let the sorted pool-specific `logP` values be

$$
\ell_{(1)} \le \ell_{(2)} \le \cdots \le \ell_{(n)}.
$$

For a clonotype with value $\ell$, the empirical lower-tail probability is

$$
p_{\mathrm{emp,low}}(\ell)
=
\frac{1}{n}
\sum_{j=1}^{n}
\mathbf{1}\{\ell_{(j)} \le \ell\}.
$$

In the implementation this is computed with a right-sided empirical CDF.

The upper-tail probability is then

$$
p_{\mathrm{emp,high}}(\ell)
=
1 - p_{\mathrm{emp,low}}(\ell).
$$

Depending on the user-selected option `tail`, the pipeline defines

$$
p\_\mathrm{value}
=
\begin{cases}
p_{\mathrm{emp,low}} & \text{if tail = low}, \\
p_{\mathrm{emp,high}} & \text{if tail = high}.
\end{cases}
$$

Final observability is defined as

$$
\mathrm{observable}
=
\mathbf{1}\{p\_\mathrm{value} < \alpha\}.
$$

Thus, observability is not based directly on raw counts, but on the position of a clonotype within the empirical null distribution calibrated from technical replicates.

---

## Alpha-Sensitivity Analysis

The pipeline can optionally evaluate the robustness of observability calls over a range of significance thresholds, by default

$$
\alpha \in \{10^{-2}, 5 \times 10^{-3}, 10^{-3}\}.
$$

For each $\alpha$, it computes:

- the fraction of observable rows
- the fraction of clonotypes ever observable
- the distribution of the number of observable timepoints per clonotype
- pairwise stability between different $\alpha$ values

This allows one to assess whether downstream conclusions are sensitive to the chosen significance threshold.

---

## Aggregation of Noise Parameters

For reporting purposes, the fitted noise parameters can be aggregated across QC-pass successful pairs.

If the null choice is `subject`, aggregation is performed per subject using medians across pairs.

Otherwise, a global median summary is produced.

This yields compact descriptors of the experimental noise landscape after QC filtering.

---

## Trajectory-Ready Outputs

After denoising, the pipeline aggregates the per-clonotype observability information into trajectory-ready tables.

For each subject, clonotype, and timepoint, it computes:

- whether the clonotype is observable at that time
- the minimum lower-tail empirical probability across entries

This yields a binary observability time series

$$
o_i(t) \in \{0,1\}
$$

for each clonotype $i$.

A second table records, for each clonotype, the number of timepoints at which it is observable:

$$
n_{i,\mathrm{obs}} = \sum_t o_i(t).
$$

These outputs are specifically designed to support downstream construction of longitudinal clonotype trajectories and transition datasets.

---

## Relation to NoisET

This pipeline should be understood as an adaptation of the NoisET statistical framework of Mora et al., not as an unrelated de novo method.

The core NoisET logic retained here is:

1. infer technical noise from replicate pairs  
2. describe latent frequencies through a heavy-tailed prior  
3. use a Negative Binomial observation model  
4. integrate over latent frequencies to obtain clonotype-level probabilities  

The present implementation extends that framework in several directions relevant for longitudinal immune-repertoire analysis:

- external QC-based inclusion of replicate pairs
- automated handling of longitudinal file structures
- explicit null pooling at clone, subject, or global level
- empirical denoising through pool-wise ECDFs
- trajectory-ready observability summaries
- optional alpha-sensitivity diagnostics

Thus, the script preserves the fundamental probabilistic rationale of NoisET while embedding it in a broader longitudinal analysis workflow.

---

## Conceptual Role in the Full Analysis Pipeline

The complete transformation implemented by the pipeline can be summarized as

$$
\text{raw replicate counts}
\;\rightarrow\;
\text{pairwise noise-model fit}
\;\rightarrow\;
\text{per-clonotype } \log P
\;\rightarrow\;
\text{empirical denoising}
\;\rightarrow\;
\text{trajectory observability}.
$$

This transformation is essential because downstream dynamical analyses should not be performed directly on raw counts. Instead, they should be based on clonotypes whose detectability has been calibrated against an experimentally inferred null model.

---

## Interpretation

The central output of the pipeline is not merely a filtered repertoire, but a probabilistically calibrated representation of clonotype observability.

Each clonotype is evaluated with respect to a noise model inferred from technical replicates. This makes it possible to separate two distinct components of observed variability:

- fluctuations that are statistically compatible with technical noise
- fluctuations that remain detectable after calibration against that noise

In this sense, the pipeline provides the statistical foundation for subsequent inference on clonotype turnover, persistence, and stochastic dynamics over time.

---

## Conclusion

The script `noise_pipeline.py`, together with the core inferential module `noiseK_nb.py`, implements a complete noise-aware framework for longitudinal RepSeq preprocessing.

Its essential components are:

- QC-based selection of reliable replicate pairs
- likelihood-based inference of technical noise using a Negative Binomial model with power-law prior
- extraction of per-clonotype marginal log-probabilities
- empirical denoising through pooled null distributions
- generation of trajectory-ready observability outputs

Because the model is explicitly adapted from the NoisET framework of Mora et al., the pipeline retains a principled probabilistic interpretation while being directly usable for longitudinal clonotype-dynamics studies.# Supplementary Methods

\newpage

#  Construction of longitudinal clonotype trajectories

## Overview

The script `trajectory_builder.py` builds a standardized long-format trajectory table from a denoised per-clonotype input table. Its role in the workflow is purely structural: it converts one or more rows per clonotype, subject, and time point into a single longitudinal representation that can be used downstream for transition construction and dynamical inference.

It performs only formatting and aggregation:

- selection of one frequency column via `--freq-col`
- retention of the required identifiers `subject`, `time`, and `aaSeqCDR3`
- aggregation of duplicate rows at the same subject-clonotype-time combination
- computation of a log-transformed abundance variable
- preservation of the `observable` column if present, or creation of that column with default `False` if absent

The script does **not** perform biological filtering, dropout correction, transition construction, detectability thresholding, or model fitting. This documentation is based on the uploaded file `3-trajectory_builder.py`.

## Input structure

The script expects a CSV input such as `per_clone_denoised_*.csv` with at least:

- `subject`
- `time`
- `aaSeqCDR3`

and one frequency column selected at runtime with `--freq-col` (default: `freq_geo_x`).

An optional column may also be present:

- `observable`

The script stops with an error if any required identifier is missing or if the requested frequency column is not found. This ensures that trajectory building is deterministic and reproducible.

## Command-line arguments

### `--per-clone`

Path to the input denoised per-clone table.

### `--outdir`

Directory in which outputs are written.

### `--freq-col`

Name of the frequency column to use as the abundance variable.

### `--epsilon`

Small positive constant added before log transformation. Default: `1e-15`.

The log-transformed frequency is defined as

$$
\mathrm{log\_freq}_i = \ln(f_i + \varepsilon),
$$

where $f_i$ is the selected frequency value for row $i$, and $\varepsilon > 0$ is the user-defined constant.

## Processing steps

### Input reading and validation

The input CSV is loaded into a pandas DataFrame. The script checks for the presence of the required identifiers

$$
\{\texttt{subject},\ \texttt{time},\ \texttt{aaSeqCDR3}\}.
$$

It also checks that the column specified by `--freq-col` exists. If not, execution stops.

### Type coercion of identifiers

The script standardizes the key identifiers as follows:

- `subject` is coerced to numeric and cast to integer
- `time` is coerced to numeric
- `aaSeqCDR3` is cast to string

This prevents mixed internal representations such as `"1"` and `1` from being treated as different subjects.

### Frequency validity filtering

The selected frequency column is converted to numeric with coercion of invalid entries to missing values. Rows are retained only when the selected frequency is finite and non-negative.

If $f_i$ denotes the selected raw frequency value for row $i$, the row is retained only if

$$
f_i \in \mathbb{R}, \qquad f_i \ge 0, \qquad |f_i| < \infty.
$$

Rows with missing, non-numeric, infinite, or negative frequency values are discarded. No additional abundance filter is applied. Zero values are allowed.

This is the only direct row-level filtering step. It is a mathematical validity check, not a biological selection rule.

### Standardization of the abundance column

After filtering, the selected frequency values are copied into a standardized output column named `freq`. This makes the downstream schema independent of the original frequency column name.

### Treatment of `observable`

If `observable` exists in the input, the script converts it robustly to Boolean. String values such as `true`, `t`, `1`, `yes`, and `y` are interpreted as `True`; all other values are interpreted as `False`. If the input column is already Boolean, it is preserved.

If `observable` is absent, the script creates it and assigns `False` to all rows.

This guarantees that the output trajectory table always contains an `observable` column, even if observability has not yet been propagated from upstream steps.

### Aggregation over duplicate subject-clone-time entries

For each unique triplet

$$
(\texttt{subject},\ \texttt{aaSeqCDR3},\ \texttt{time}),
$$

the script keeps one row only. If duplicate rows exist, it aggregates them as:

- `freq` = arithmetic mean
- `observable` = maximum

Suppose a given subject-clone-time triplet has $m$ rows with frequencies $f_1, \dots, f_m$ and Boolean observability indicators $o_1, \dots, o_m$, where each $o_j \in \{0,1\}$. The aggregated values are

$$
\bar f = \frac{1}{m}\sum_{j=1}^{m} f_j,
$$

and

$$
o^* = \max_{1 \le j \le m} o_j.
$$

Thus the aggregated row is considered observable if any duplicated row is observable.

This aggregation is structural rather than inferential. The script is not fitting a model to reconcile duplicates; it is simply enforcing a unique row per subject, clonotype, and time point.

### Log transformation

After aggregation, the script computes

$$
\texttt{log\_freq} = \ln(\texttt{freq} + \varepsilon).
$$

This step is important because downstream dynamical analyses are typically carried out in log-abundance space. The addition of $\varepsilon$ prevents undefined values when `freq = 0`:

$$
\ln(0 + \varepsilon) = \ln(\varepsilon).
$$

The script uses the natural logarithm. Therefore downstream displacements are naturally written as

$$
\Delta x = \ln f(t_1) - \ln f(t_0),
$$

which, when $\varepsilon$ is negligible relative to the frequency scale, approximates the log fold-change

$$
\Delta x \approx \ln\left(\frac{f(t_1)}{f(t_0)}\right).
$$

This representation is especially useful when clone frequencies span multiple orders of magnitude.

## Statistical role of the script

The script itself is not a statistical inference engine, but its design choices affect the structure of the downstream dataset.

### No abundance thresholding

No detectability threshold, abundance cutoff, or clone-selection rule is applied at this step. All mathematically valid denoised observations are retained.

### No model-based smoothing

The script does not shrink or smooth frequencies. It passes them through directly after duplicate aggregation.

### No dropout correction

Missing time points are not imputed, and absent rows are not interpreted biologically at this stage.

### No transition definition

The output is a state table, not a transition table. Pairwise temporal transitions must be built later from this long-format representation.

This separation is methodologically useful because it keeps structural preprocessing distinct from dynamical assumptions.

## Output files

### `trajectories_long.csv`

This is the main output. It contains one row per unique `(subject, aaSeqCDR3, time)` combination after aggregation, with columns:

- `subject`
- `aaSeqCDR3`
- `time`
- `freq`
- `observable`
- `log_freq`

Rows are sorted by `subject`, `aaSeqCDR3`, and `time`.

### `included_clones.csv`

This file contains the unique set of `(subject, aaSeqCDR3)` pairs represented after aggregation.

If $N_{\mathrm{clone}}$ denotes the number of unique subject-clone pairs, then

$$
N_{\mathrm{clone}} =
\#\{(s,c): \exists t \text{ such that } (s,c,t) \text{ is present in } \texttt{trajectories\_long.csv}\}.
$$

### `trajectory_build_report.md`

This report contains run metadata and basic output sizes, including:

- timestamp
- input file path
- selected frequency column
- epsilon value
- number of rows after frequency validity filtering
- number of rows after aggregation
- number of subjects
- number of unique subject-clone pairs
- number of rows in the final trajectory table

This file is useful for provenance and auditability.

## Interpretation of counts

Let

- $N_{\mathrm{in}}$ = number of rows in the raw input
- $N_{\mathrm{valid}}$ = number of rows after removing invalid frequencies
- $N_{\mathrm{agg}}$ = number of unique `(subject, aaSeqCDR3, time)` combinations after aggregation
- $N_{\mathrm{clone}}$ = number of unique `(subject, aaSeqCDR3)` pairs

Then in general,

$$
N_{\mathrm{agg}} \le N_{\mathrm{valid}} \le N_{\mathrm{in}},
$$

with strict inequality whenever duplicate rows or invalid frequencies are present.

The final `trajectories_long.csv` table contains exactly $N_{\mathrm{agg}}$ rows.

## Why log-frequency trajectories are useful downstream

If $f(t)$ denotes clonotype frequency at time $t$, defining

$$
x(t) = \ln(f(t) + \varepsilon)
$$

allows temporal changes to be written as

$$
\Delta x = x(t_1) - x(t_0).
$$

When $\varepsilon$ is small relative to $f$, this is approximately the log fold-change. This is advantageous because repertoire frequencies usually span orders of magnitude and many downstream dynamical models are naturally formulated in log-abundance space.

## Limitations

This script is intentionally minimal.

First, duplicate rows are collapsed using the arithmetic mean of frequency values. This is convenient and deterministic, but it is not the only possible summary.

Second, `observable` is propagated by the maximum. This preserves any evidence of observability but does not quantify uncertainty.

Third, zero frequencies are permitted and mapped to $\ln(\varepsilon)$. This is numerically convenient, but the biological interpretation depends on upstream preprocessing.

Fourth, missing time points are not imputed. Structural absence from the trajectory table should therefore not automatically be interpreted as biological extinction.

## Reproducibility

The script is fully deterministic. Given the same input file and the same arguments, it produces the same outputs. No random seed is required because no stochastic step is involved.

This makes it an appropriate preprocessing stage before downstream procedures that may involve thresholding, bootstrap resampling, or model fitting.

## Practical summary

In summary, `3-trajectory_builder.py`:

1. validates required identifiers and the selected frequency column
2. standardizes identifier types
3. removes rows with invalid frequency values
4. standardizes the selected abundance variable as `freq`
5. preserves or creates `observable`
6. collapses duplicate subject-clone-time entries by mean frequency and max observability
7. computes `log_freq = ln(freq + epsilon)`
8. exports a long-format trajectory table, an included-clone list, and a provenance report

It is therefore a foundational structural step for transition building and downstream dynamical analysis, while deliberately avoiding any inferential decision beyond basic mathematical consistency.

\newpage

# Construction of Clonotype Transition Patterns

The script `3-build_transition_patterns.py` constructs a dataset of clonotype frequency transitions over time from longitudinal repertoire data.

Its purpose is to convert a trajectory-based representation of clonotype dynamics into a transition-based representation suitable for stochastic inference, including the estimation of drift and diffusion functions.

---

## Input Data Structure

The input is a long-format table where each row represents a clonotype at a given timepoint. The required columns are:

- subject: subject identifier  
- aaSeqCDR3: clonotype identifier  
- time: discrete time index  
- freq: clonotype frequency  
- observable: detectability indicator (0/1 or True/False)  

An optional column:

- log_freq: precomputed log-frequency  

Column names can be remapped via command-line arguments.

---

## Definition of the Dynamical Variable

The analysis is performed in log-frequency space. The state variable is:

$$
x(t) = \log f(t)
$$

If the column `log_freq` is present and numeric, it is used directly as the state variable.

Otherwise, the logarithm of frequency is computed with a lower bound:

$$
x(t) = \log\big(\max(f(t), f_{\min})\big)
$$

### Minimum frequency parameter

The parameter `min_freq` defines the lower bound:

$$
f(t) \leftarrow \max(f(t), f_{\min})
$$

Default value:

$$
f_{\min} = 10^{-12}
$$

This prevents numerical instability in the low-frequency regime.

---

## Transition Definition

A transition is defined between two timepoints:

$$
(t_0, t_1), \quad t_1 > t_0
$$

with time interval:

$$
\Delta t = t_1 - t_0
$$

and state variables:

$$
x_0 = x(t_0), \quad x_1 = x(t_1)
$$

The displacement is:

$$
\Delta x = x_1 - x_0
$$

Each transition is represented as:

$$
(x_0, \Delta x, \Delta t)
$$

---

## Temporal Constraint

The parameter `max_dt` limits the maximum time lag:

$$
\Delta t \leq \Delta t_{\max}
$$

- Default: 6  
- If `max_dt <= 0`, no filtering is applied  

This parameter controls the temporal scale of the transitions included in the analysis.

---

## Transition Patterns

By default, only the ALL transition pattern is constructed. The ADJ pattern is generated only if explicitly requested at runtime using the `--build_adj` flag.

### All pairs (all)

All ordered pairs of timepoints are included:

$$
\forall i < j
$$

The number of transitions per trajectory with $ n $ timepoints is:

$$
N = \binom{n}{2}
$$

This pattern captures both short- and long-range temporal dynamics.

---

### Consecutive observed pairs (adj, optional)

Only transitions between consecutive observed timepoints are included:

$$
(t_i, t_{i+1})
$$

In general:

$$
\Delta t \neq 1
$$

because intermediate timepoints may be missing.

This pattern focuses on locally observed dynamics.

---

## Observation Classes

Each timepoint has a detectability state:

$$
o(t) \in \{0,1\}
$$

Missing values are treated conservatively as:

$$
o(t) = 0
$$

Transitions are classified according to endpoint observability:

- TT: observed → observed  
- TF: observed → not observed  
- FT: not observed → observed  
- FF: not observed → not observed  

---

## Weighting System

If enabled, each transition is assigned a weight:

$$
w = w_{\text{clone}} \cdot w_{\Delta t} \cdot w_{\text{class}}
$$

This weighting scheme corrects structural imbalances in the dataset.


-Clone Weight

The clone-level weight is:

$$
w_{\text{clone}} = \frac{1}{N_{\text{transitions per clone}}}
$$

This prevents trajectories with many timepoints from dominating the dataset.


-Time-Lag Weight

If the option `balance_dt` is enabled:

$$
w_{\Delta t} = \frac{1}{N(\Delta t)}
$$

Otherwise:

$$
w_{\Delta t} = 1
$$

This compensates for uneven sampling across time intervals.


- Class Weights

The class-dependent weight is:

$$
w_{\text{class}} =
\begin{cases}
1 & \text{TT} \\
\gamma_{\text{cross}} & \text{TF, FT} \\
\gamma_{FF} & \text{FF}
\end{cases}
$$

Default values:

$$
\gamma_{\text{cross}} = 0.25, \quad \gamma_{FF} = 0.1
$$

Interpretation:

- TT transitions are fully observed and most reliable  
- TF and FT transitions are partially observed (censored)  
- FF transitions correspond to unobserved states  

-Weight Normalization

If enabled, each weight component is normalized:

$$
w_i \leftarrow \frac{w_i}{\mathbb{E}[w_i]}
$$

ensuring that:

$$
\mathbb{E}[w_i] \approx 1
$$

This stabilizes numerical behavior while preserving relative weights.

---

## Pattern-Specific Behavior

- The ALL pattern supports full weighting, including time-lag balancing  
- The ADJ pattern typically does not apply time-lag weighting  

---


**Default behavior**

- The ALL pattern is always generated  
- The ADJ pattern is generated only if `--build_adj` is specified  

This design reflects the primary use of ALL transitions for dynamical inference, while ADJ transitions are used for validation or sensitivity analyses.


## Output Files

The script produces by default:

- transitions_all.csv  

If `--build_adj` is specified, it additionally produces:

- transitions_adj.csv  

Each row contains:

$$
(t_0, t_1, \Delta t, x_0, x_1, \Delta x)
$$

as well as metadata such as observation class, pattern type, and optional weights.



- Count Tables

The script also outputs aggregated counts (for each generated pattern):

$$
N(\text{subject}, \text{obs\_class}, \Delta t)
$$

These are used for diagnostic and stratification purposes.

---

## Stochastic Interpretation

The transitions can be interpreted as samples from a stochastic process:

$$
x(t + \Delta t) = x(t) + b(x(t)) \, \Delta t + \sqrt{2 D(x(t)) \, \Delta t} \, \eta
$$

where $\eta \sim \mathcal{N}(0,1)$.


This implies:

$$
b(x) \approx \frac{\mathbb{E}[\Delta x \mid x]}{\Delta t}
$$

$$
D(x) \approx \frac{\mathrm{Var}[\Delta x \mid x]}{2\Delta t}
$$

---

## Statistical Dependencies and Bias Control

The construction of transition datasets introduces statistical dependencies and sampling biases.

### Dependence in the ALL pattern

Multiple transitions may share the same underlying observations, for example:

$$
(t_0 \to t_2), \quad (t_0 \to t_3), \quad (t_1 \to t_3)
$$

This induces correlation among transitions and violates independence assumptions.

---

### ADJ pattern

Restricting to consecutive transitions reduces overlap but does not eliminate dependence between transitions from the same clonotype.

---

### Sampling imbalance

The dataset is intrinsically heterogeneous:

$$
P(\Delta t), \quad P(\text{clone}), \quad P(\text{obs\_class})
$$

are not uniform.

---

### Role of weighting

The weighting scheme:

$$
w = w_{\text{clone}} \cdot w_{\Delta t} \cdot w_{\text{class}}
$$

reduces bias but does not restore statistical independence.

---

### Implications for inference

Estimated quantities such as drift and diffusion should be interpreted as effective estimators under dependence rather than independent samples.

---

### Cluster bootstrap

Reliable uncertainty estimation requires resampling at the level of:

$$
(\text{subject}, \text{clonotype})
$$

to preserve the dependence structure.

---

### Practical consequence

Valid inference requires explicit handling of dependence and appropriate resampling strategies.

---

## Role in the Analysis Pipeline

The script performs the transformation:

$$
\text{trajectories} \rightarrow \text{transitions} \rightarrow \text{stochastic inference}
$$

It is a central step connecting denoised trajectories to dynamic modeling.

---

## Conclusion

The script transforms longitudinal clonotype data into a transition-based representation:

$$
\{x(t)\} \;\longrightarrow\; \{(x_0, \Delta x, \Delta t, w)\}
$$

This provides the minimal sufficient structure for quantitative modeling of clonotype dynamics in a stochastic framework.

\newpage

# Model-agnostic local diagnostics and finite-time inference

The script `5-finite_time_diagnostics.py` performs **model-agnostic characterization of clonotype dynamics** from a transition dataset. It implements two complementary analysis blocks:

1. **Local diagnostics**: estimation of conditional displacement statistics and variance structure as a function of initial abundance and time lag.
2. **Finite-time inference**: extraction of effective drift and diffusion parameters from transition statistics across multiple time intervals.

All computations are performed in **natural logarithm space (ln)**, regardless of the input representation.

The script operates on transition-level data generated from clonotype trajectories and provides the statistical foundation for subsequent dynamical modeling.

---

## Input Data Structure

The input is a transition table with at least the following columns:

- `subject`: subject identifier  
- `aaSeqCDR3`: clonotype identifier  
- `dt`: time lag between observations  
- `x0`: log-frequency at initial time  
- `x1`: log-frequency at final time  
- `dx`: displacement $x_1 - x_0$  
- `obs_class`: detectability class  

Optional:

- `w`: transition weight  

---

## Log-space normalization

The script accepts input in either:

- natural log (`ln`)
- base-10 log (`log10`)

If input is in base-10:

$$
x^{(\ln)} = x^{(\log_{10})} \cdot \ln(10)
$$

All outputs are expressed in **natural log space**:

- $\ln f_0$
- $\ln f_1$
- $\Delta x$

Additionally, the geometric mean abundance is defined as:

$$
\ln f_{\mathrm{geo}} = \frac{1}{2}(\ln f_0 + \ln f_1)
$$

---

## Filtering and selection

Transitions are filtered according to:

### Time constraints

$$
dt > dt_{\min}, \quad dt \le dt_{\max}
$$

Optional exact selection:

$$
dt \in \{dt_1, dt_2, \ldots\}
$$

### Detectability

- `--only_TT`: restrict to TT transitions  
- or specify classes explicitly  

---

## Weighting

If a weight column is provided:

$$
w_i > 0
$$

weights are incorporated in all estimators.

If `weight_col = "1"`:

$$
w_i = 1
$$

---

## Binning strategy

Two variables are used:

- $\ln f_0$ (initial abundance)
- $\ln f_{\mathrm{geo}}$ (mean abundance)

Bins are defined using:

### Quantile binning

$$
\text{edges} = \mathrm{quantile}(x, q)
$$

ensuring approximately equal sample sizes per bin.

### Fixed binning

User-supplied bin edges can be used for consistency across datasets.

---

## Local diagnostics

### Drift estimation

For each bin:

$$
b(x_0, dt) \approx \mathrm{median}(\Delta x \mid x_0, dt)
$$

Additionally:

$$
P(\Delta x > 0 \mid x_0, dt)
$$

is computed.

Bootstrap resampling provides confidence intervals.

---

### Mean squared displacement (MSD)

Global:

$$
\mathrm{MSD}(dt) = \mathbb{E}[\Delta x^2 \mid dt]
$$

Binned:

$$
\mathrm{MSD}(x_0, dt), \quad \mathrm{MSD}(f_{\mathrm{geo}}, dt)
$$

---

### Variance estimation via MAD

Robust variance is estimated using:

$$
\mathrm{MAD} = \mathrm{median}(|\Delta x - \mathrm{median}(\Delta x)|)
$$

Converted to standard deviation:

$$
\sigma \approx 1.4826 \cdot \mathrm{MAD}
$$

Variance:

$$
\mathrm{Var} \approx \sigma^2
$$

---

## Finite-time inference

### Zero-crossing point

The drift curve satisfies:

$$
b(x^*) = 0
$$

The zero-crossing $x^*$ is identified by sign change.

Interpolation:

$$
x^* = x_i - \frac{b(x_i)}{b(x_{i+1}) - b(x_i)} (x_{i+1} - x_i)
$$

---

### Local slope estimation

A local linear fit is performed:

$$
b(x) \approx a + s(x - x^*)
$$

Slope:

$$
s = \frac{db}{dx}
$$

---

### Characteristic timescale

$$
\tau = \frac{dt}{-s}
$$

Valid only when:

$$
s < 0
$$

---

### Diffusion coefficient

$$
D = \frac{\mathrm{Var}(\Delta x)}{2dt}
$$

---

## Bootstrap procedure

Cluster bootstrap is used with clustering over:

$$
(\text{subject}, \text{clonotype})
$$

This preserves temporal dependencies.

Confidence intervals are computed as:

$$
[\alpha/2, 1 - \alpha/2]
$$

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

$$
x(t + dt) = x(t) + b(x)dt + \sqrt{2D(x)dt}\,\eta
$$

with:

$$
\eta \sim \mathcal{N}(0,1)
$$

Estimators:

$$
b(x) \approx \frac{\mathbb{E}[\Delta x]}{dt}
$$

$$
D(x) \approx \frac{\mathrm{Var}(\Delta x)}{2dt}
$$

---

## Bias and dependence

### Non-independence

Transitions share observations:

$$
(t_0 \to t_2), (t_0 \to t_3)
$$

### Heterogeneity

$$
P(dt), P(x), P(\text{class})
$$

are non-uniform.

### Implication

Bootstrap at cluster level is required.

---

## Role in pipeline

This script performs:

$$
\text{transitions} \rightarrow \text{moments} \rightarrow \text{dynamical inference}
$$

\newpage

# finite-time validation

The script `6-finite_time_fit.py` performs **finite-time validation of an Ornstein-Uhlenbeck (OU) model** using empirical estimates obtained from transition data.

It does not estimate drift/diffusion from raw data; instead, it:

1. Takes as input summary statistics
2. Calibrates an OU process at a reference timescale
3. Predicts variance across time
4. Compares predictions to empirical data

---

## OU model

$$
dx = -\lambda x \, dt + \sqrt{2D} \, dW_t
$$

---

## Calibration

$$
\lambda = -\text{slope}(dt_{\mathrm{ref}})
$$

$$
D = \lambda \cdot \frac{\mathrm{Var}(dt_{\mathrm{ref}})}{1 - e^{-2\lambda dt_{\mathrm{ref}}}}
$$

---

## Predictions

$$
\mathrm{Var}(dt) = \frac{D}{\lambda}(1 - e^{-2\lambda dt})
$$

$$
D_{\mathrm{app}} = \frac{\mathrm{Var}(dt)}{2dt}
$$

---

## Outputs

- `ou_pred_vs_empirical.csv`
- `ou_fit_meta.json`

---

## Interpretation

Agreement with OU:

- variance saturates
- apparent diffusion decreases

---

## Role

Validation step linking empirical dynamics to OU theory.

\newpage

# Drift model fitting and diffusion-from-residuals

The script `7-fit_drift_models_analysis_and_residuals.py` performs **parametric fitting of drift functions** on transition-level data at a fixed lag $dt$, and derives **diffusion estimates from model residuals**. It is **analysis-only** (no plots), producing tabular outputs for downstream diagnostics and figures.

Core steps:
1. Filter transitions (by $dt$, class, weights)
2. Fit multiple drift models $\hat b(x_0)$ to displacements $\Delta x$
3. Compare models via AIC/BIC
4. Compute per-transition residuals and squared residuals
5. Estimate **conditional diffusion** $D(x_0)$ from residuals (binned)
6. (Optional) Cluster bootstrap for parameter and model-comparison uncertainty

---

## Input

A transitions CSV with at least:
- `dt`, `x0`, `dx`

Optional:
- `obs_class` (for filtering)
- weight column (e.g. `w_clone`)
- identifiers (`subject`, `aaSeqCDR3`, `time0`, `time1`)

Filtering:

$$
dt = \text{--dt}, \qquad \text{obs\_class} = \text{--filter\_class (e.g. TT)}
$$

Finite rows only are retained.

---

## Models

The script fits a **set of drift models** $\hat b(x) \approx \mathbb{E}[\Delta x \mid x_0=x]$.

### OU-linear

$$
\hat b(x) = a + b x
$$

Parameters: $a,b$.  
Equilibrium: $x^* = -a/b$.  
Timescale: $\tau = -1/b$ (for $b<0$).

### Tanh-saturating

$$
\hat b(x) = A \tanh\!\left(\frac{x^* - x}{\Delta}\right)
$$

### Rational-saturating

$$
\hat b(x) = A \frac{z}{\sqrt{1+z^2}}, \quad z=\frac{x^*-x}{\Delta}
$$

### Hinge-plateau (smooth)

Smooth interpolation between linear regime and plateau:

$$
\hat b(x) = (a+bx)\,s(x) + (-A)\,[1-s(x)], \quad s(x)=1-\sigma\!\left(\frac{x-x_c}{\Delta}\right)
$$

---

## Estimation

### Loss function

All models are fit by minimizing weighted SSE:

$$
\mathrm{SSE} = \sum_i w_i (\Delta x_i - \hat b(x_{0,i}))^2
$$

- Unweighted: $w_i=1$
- Weighted: sanitized and clipped weights (quantile cap `--weight_cap_q`), then rescaled to mean 1.

Robust optimization uses `least_squares(..., loss="soft_l1")` for nonlinear models.

### Standardization (OU-linear)

For numerical stability:

$$
z = \frac{x-\mu}{\sigma}
$$

Fit in $z$-space, then map back to $x$-space.

---

## Model comparison

Information criteria (Gaussian residual approximation):

$$
\mathrm{AIC} = n\ln(\mathrm{SSE}/n) + 2k,\qquad
\mathrm{BIC} = n\ln(\mathrm{SSE}/n) + k\ln n
$$

Outputs:
- `fit_summary.csv` (per-model parameters, SSE, AIC/BIC)
- `model_compare.csv` ($\Delta$AIC/$\Delta$BIC vs baseline)

---

## Fitted values and residuals

For each model:

$$
\hat{\Delta x}_i = \hat b(x_{0,i}),\qquad
r_i = \Delta x_i - \hat{\Delta x}_i,\qquad
r_i^2
$$

Saved in:
- `fitted_values.csv`

---

## Diffusion from residuals

Assuming:

$$
\Delta x = b(x_0)\,dt + \epsilon, \quad \mathbb{E}[\epsilon]=0,\ \mathrm{Var}(\epsilon)=2D(x_0)dt
$$

Estimate diffusion via residuals.

### Mean-based

$$
\mathbb{E}[r^2 \mid x_0] \approx 2D(x_0)dt
$$

### Robust (MAD)

$$
\mathrm{MAD} = \mathrm{median}(|r - \mathrm{median}(r)|),\quad
\sigma \approx 1.4826\,\mathrm{MAD},\quad
\mathrm{Var}_{\mathrm{MAD}} = \sigma^2
$$

---

## Binning

Quantile bins over $x_0$ with constraints:
- `nbins_diff`
- `min_per_bin_diff`

Per-bin outputs:
- median residual
- MAD
- robust variance (`var_mad`)
- mean $r^2$

Saved in:
- `diffusion_binned.csv`

Long-format residuals:
- `diffusion_residuals_long.csv`

---

## Drift grid

A regular grid $x\in[q_{lo},q_{hi}]$ is built:

$$
x \in [Q_{q_{lo}}(x_0),\, Q_{q_{hi}}(x_0)]
$$

Predictions $\hat b(x)$ for each model are saved in:
- `drift_function_grid.csv`

---

## Bootstrap (optional)

Cluster bootstrap over:

$$
(\text{subject}, \text{aaSeqCDR3})
$$

Outputs:
- `bootstrap_reps.csv`
- `bootstrap_summary.csv`

Statistics:
- parameter medians and CI
- $\Delta$AIC/$\Delta$BIC distributions

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
- **Diffusion estimates** reveal heteroscedasticity $D(x)$

Key diagnostic:
- If OU is insufficient, nonlinear models reduce SSE and AIC
- Residual variance decreasing with $x$ suggests abundance-dependent noise

---

## Limitations

- Fixed $dt$ analysis
- Gaussian residual assumption in IC
- Residual-based diffusion assumes correct drift model
- Weighting may affect variance estimates

---

## Role in pipeline

$$
\text{transitions} \rightarrow \text{drift fit} \rightarrow \text{residuals} \rightarrow \text{diffusion}
$$

\newpage

# Fokker–Planck closed-vs-open diagnostics in log-frequency space


The script `8-fp_diagnostics_build.py` constructs diagnostic grids for comparing empirical clonotype dynamics with the stationary solution of a one-dimensional Fokker–Planck equation in log-frequency space,

$$
x = \ln f.
$$

Its purpose is not to fit new dynamical models, but to combine previously estimated drift and diffusion functions with the empirical distribution of clonotype abundances and thereby evaluate whether the observed dynamics are compatible with a closed stationary process or instead show signatures of open, non-conservative behavior.

It reads:

- a drift grid $b(x)$ from a fitted drift model,
- a binned diffusion estimate $D(x)$,
- an empirical trajectory table used to reconstruct the observed density $p_{\mathrm{emp}}(x)$,

and returns a unified grid containing:

- the model drift $b(x)$,
- the interpolated diffusion $D(x)$,
- the empirical density $p_{\mathrm{emp}}(x)$,
- the closed stationary density $p_{\mathrm{closed}}(x)$,
- the density ratio $p_{\mathrm{emp}}/p_{\mathrm{closed}}$,
- the empirical probability current $J(x)$,
- the source term $S(x)=dJ/dx$.

This documentation is based on the uploaded script `8-fp_diagnostics_build.py`.

## Conceptual role in the workflow

Upstream steps provide:

1. a fitted deterministic drift function $b(x)$,
2. a diffusion estimate $D(x)$ from residual fluctuations,
3. a trajectory table giving empirical occupancy of state space.

This script combines these ingredients to test two related questions:

1. **Closed stationary consistency**  
   If the repertoire behaved as a stationary one-dimensional diffusion with drift $b(x)$ and diffusion $D(x)$, what stationary density would it imply?

2. **Open-process diagnostics**  
   Given the empirical density $p_{\mathrm{emp}}(x)$, what probability current $J(x)$ and source term $S(x)$ would be required to sustain it under the same drift–diffusion field?

In this sense, the script bridges

$$
\text{drift fit} + \text{diffusion estimate} + \text{empirical occupancy}
\;\longrightarrow\;
\text{stationary FP diagnostics}.
$$

## Inputs

The script requires three CSV files.

###  `drift_function_grid.csv`

This file contains a grid of $x$ values and one or more fitted drift functions. The relevant model-specific drift column must be one of:

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

### `diffusion_binned.csv`

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

### `trajectories_long.csv`

This file provides empirical occupancy of state space. The script extracts the log-frequency coordinate from:

1. `log_freq`
2. `x`
3. `log(freq)` if only `freq` is available

Optionally, it can restrict the empirical density to rows with pointwise `observable == True`.

## State-space representation

All calculations are performed in log-frequency space:

$$
x = \ln f.
$$

This choice is natural because:

- clone frequencies span several orders of magnitude,
- upstream drift and diffusion fits are already expressed in $x$-space,
- stationary Fokker–Planck equations are more interpretable when multiplicative abundance changes become additive.

The diagnostic grid is therefore defined over a one-dimensional axis $x$, inherited from the fitted drift grid.

## Filtering logic

The script identifies one consistent dynamical slice across drift and diffusion inputs by matching:

- a single $dt$,
- a single `obs_class`,
- a single drift model,
- and, if specified, exact weight metadata.

Formally, the diagnostics are computed for one selected combination

$$
(dt,\ \mathrm{obs\_class},\ \mathrm{model},\ \mathrm{weighting}).
$$

This is important because the drift and diffusion estimates may differ substantially across time lag, detectability class, weighting scheme, and fitted drift model. The script therefore avoids mixing incompatible slices.

## Drift extraction

After filtering, the script extracts:

- a grid $x_1,\dots,x_n$,
- a model-specific drift function $b(x)$.

The drift grid is sorted by $x$, and a minimum number of points is required in order to compute stable numerical derivatives later.

If the requested drift model is $m$, then the script uses the corresponding column:

$$
b(x) = b_m(x).
$$

No further smoothing of the drift is applied inside this script; it assumes that the upstream fit has already produced a suitable grid representation.

## Diffusion construction

The diffusion input is not necessarily provided directly as $D(x)$, but rather as a variance proxy for displacements over lag $dt$. The script converts this to diffusion using

$$
D(x) = \frac{\mathrm{Var}(\Delta x \mid x, dt)}{2\,dt}.
$$

Depending on the available column, the variance proxy is defined as:

### Option 1 — `var_mad`

A robust variance estimate already computed upstream:

$$
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathrm{var\_mad}(x).
$$

### Option 2 — `mean_resid2`

The conditional mean squared residual:

$$
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathbb{E}[r^2 \mid x].
$$

### Option 3 — `mad`

If only the median absolute deviation is available, the script uses

$$
\mathrm{Var}_{\mathrm{proxy}}(x) = \mathrm{mad}(x)^2.
$$

After constructing $D(x)$ on the diffusion bins, the script interpolates it onto the drift grid $x_{\mathrm{grid}}$. If interpolation leaves missing values at the edges, it fills them by one-dimensional interpolation across the available finite points.

Thus the final diffusion field is

$$
D_{\mathrm{on}}(x_i), \qquad i=1,\dots,n,
$$

defined on the same grid as the drift.

## Empirical density $p_{\mathrm{emp}}(x)$

The empirical density is constructed from the trajectory table by histogramming state-space occupancy onto the drift grid.

### Choice of empirical variable

If not overridden by `--x_col_traj`, the script uses:

- `log_freq`, if present,
- else `x`,
- else $\ln(\mathrm{freq})$.

Only finite values are retained.

### Optional observability filter

If `--traj_observable_only` is set, the empirical density is built only from rows satisfying

$$
\mathrm{observable} = \mathrm{True}.
$$

This is a pointwise filter on state occupancy and should not be confused with transition classes such as TT, TF, FT, FF.

### Histogram construction on a predefined grid

The script defines histogram bin edges from the drift grid midpoints. Let the grid be

$$
x_1 < x_2 < \cdots < x_n.
$$

Then the histogram edges are built as midpoint boundaries between adjacent grid points, with extrapolated outer edges at the ends. Counts are converted to a density by dividing by total counts and bin widths.

Thus the empirical density satisfies approximately

$$
\int p_{\mathrm{emp}}(x)\,dx = 1.
$$

### Optional smoothing

The empirical density can be smoothed in index-space by Gaussian convolution with standard deviation `smooth_p_sigma_pts` measured in grid points. After smoothing, the density is re-normalized by numerical integration.

## Closed stationary density $p_{\mathrm{closed}}(x)$

For a one-dimensional Fokker–Planck equation with drift $b(x)$ and diffusion $D(x)$, the stationary zero-current solution satisfies

$$
J(x)=0.
$$

Under the convention used in this script, the resulting closed stationary density is

$$
p_{\mathrm{closed}}(x)
\propto
\frac{1}{D(x)}
\exp\!\left(\int^x \frac{b(y)}{D(y)}\,dy\right).
$$

### Numerical implementation

The script computes:

1. a safe diffusion field

$$
D_{\mathrm{safe}}(x) = \max(D(x), \varepsilon),
$$

with $\varepsilon = 10^{-12}$,

2. the integrand

$$
I'(x) = \frac{b(x)}{D_{\mathrm{safe}}(x)},
$$

3. its cumulative trapezoidal integral

$$
I(x) = \int^x \frac{b(y)}{D_{\mathrm{safe}}(y)}\,dy,
$$

4. the unnormalized density

$$
p_{\mathrm{unn}}(x)
=
\frac{1}{D_{\mathrm{safe}}(x)}
\exp(I(x)).
$$

For numerical stability, the maximum of $I(x)$ is subtracted before exponentiation. Finally, the density is normalized using trapezoidal integration:

$$
p_{\mathrm{closed}}(x)
=
\frac{p_{\mathrm{unn}}(x)}
{\int p_{\mathrm{unn}}(x)\,dx}.
$$

## Empirical current $J(x)$

The script next computes the current implied by the empirical density $p_{\mathrm{emp}}(x)$ under the same drift and diffusion fields:

$$
J(x) = b(x)\,p_{\mathrm{emp}}(x) - \frac{d}{dx}\bigl[D(x)\,p_{\mathrm{emp}}(x)\bigr].
$$

This quantity measures local imbalance between deterministic transport and diffusive redistribution.

The derivative is computed numerically by `np.gradient` on the drift grid. Because numerical derivatives amplify noise, the raw current is smoothed by Gaussian convolution with width `smooth_J_sigma_pts` measured in grid points:

$$
J_s(x) = \mathcal{G}_\sigma * J(x).
$$

The smoothed current is what is written to output as `J_emp`.

## Source term $S(x)$

The source term is defined as the derivative of the smoothed current:

$$
S(x) = \frac{dJ_s(x)}{dx}.
$$

In a strictly closed stationary system with zero current, one expects

$$
J(x) = 0, \qquad S(x)=0.
$$

Nonzero $J(x)$ or $S(x)$ indicates that the empirical density cannot be explained purely as a closed stationary balance under the supplied drift and diffusion fields.

As with the current, the raw derivative is smoothed again before saving:

$$
S_s(x) = \mathcal{G}_\sigma * \frac{dJ_s(x)}{dx}.
$$

The saved column `S_emp` corresponds to this smoothed source profile.

## Density ratio

The script also computes the pointwise ratio

$$
\frac{p_{\mathrm{emp}}(x)}{p_{\mathrm{closed}}(x)}.
$$

This ratio directly identifies abundance regions where the empirical occupancy is enriched relative to the closed stationary prediction, depleted relative to it, or broadly consistent with it.

## Summary metrics: JS and KS

To compare the empirical and closed stationary densities globally, the script computes two summary metrics.

### Jensen–Shannon divergence

The densities are converted to discrete probability masses on the grid by multiplying by local grid spacing $dx$, then normalized. The Jensen–Shannon divergence is

$$
\mathrm{JS}(P,Q)
=
\frac{1}{2}\mathrm{KL}(P\|M)
+
\frac{1}{2}\mathrm{KL}(Q\|M),
\qquad
M=\frac{P+Q}{2}.
$$

### Kolmogorov–Smirnov distance

The script also computes

$$
\mathrm{KS}(P,Q)
=
\max_x
\left|
F_P(x) - F_Q(x)
\right|.
$$

Both metrics are stored in the summary output.

## Output files

The script writes three files.

### `fp_<tag>.grid.csv`

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

### `fp_<tag>.summary.json`

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

###  `fp_<tag>.summary.txt`

This is a plain-text key=value version of the same summary metadata.

## Tag construction and reproducibility

The output tag is built from the requested parameter combination, typically including:

- $dt$
- `obs_class`
- drift model
- optional weight metadata
- optional observability restriction

This ensures that diagnostic grids from different configurations do not overwrite one another and remain traceable to the exact analysis settings used.

## Numerical conventions and safeguards

Several safeguards are built into the implementation.

### Safe diffusion floor

Diffusion is clipped below by a small positive constant:

$$
D(x) \ge 10^{-12}.
$$

### Interpolation gap filling

If diffusion interpolation leaves missing values on the drift grid, the script fills them by interpolation from nearby finite points, provided enough finite support exists.

### Smoothing in grid-index space

Gaussian smoothing is performed in index space, not in physical $x$-space. Therefore the smoothing width is measured in grid points, not in units of $\ln f$.

### Normalization by trapezoidal integration

All densities are normalized by numerical integration using the trapezoidal rule.

## Statistical interpretation

This script does not fit a new Fokker–Planck model. Instead, it asks whether the previously fitted drift and diffusion functions are jointly consistent with the observed abundance distribution.

### Closed-process interpretation

If

$$
p_{\mathrm{emp}}(x) \approx p_{\mathrm{closed}}(x),
\qquad
J(x) \approx 0,
\qquad
S(x) \approx 0,
$$

then the observed density is broadly compatible with a closed stationary diffusion.

### Open-process interpretation

If instead one finds substantial deviations such as

$$
p_{\mathrm{emp}}(x) \neq p_{\mathrm{closed}}(x),
\qquad
J(x)\neq 0,
\qquad
S(x)\neq 0,
$$

then the empirical density behaves as if additional mechanisms are present, such as influx or efflux of probability mass, nonstationarity, misspecified drift, misspecified diffusion, or other violations of the closed one-dimensional approximation.

## Limitations

Several limitations should be kept in mind.

- The diagnostics depend directly on upstream drift and diffusion estimates.
- The method assumes a one-dimensional diffusion in $x=\ln f$.
- The empirical density depends on histogram construction, smoothing, and observability filtering.
- Current and source estimates involve numerical derivatives and are noise-sensitive.
- This script produces point estimates only and does not propagate bootstrap uncertainty.

## Practical summary

In practical terms, `8-fp_diagnostics_build.py` performs the following sequence:

1. select one drift model and one dynamical slice,
2. reconstruct $D(x)$ from a residual-based variance proxy,
3. build the empirical state-density $p_{\mathrm{emp}}(x)$,
4. compute the closed stationary prediction $p_{\mathrm{closed}}(x)$,
5. derive the empirical current $J(x)$,
6. derive the source term $S(x)$,
7. summarize agreement between empirical and closed stationary densities by JS and KS.

This makes the script a central step for evaluating whether the inferred clonotype dynamics are compatible with a closed stationary Fokker–Planck description.

## Conclusion

This script synthesizes fitted drift, residual-based diffusion, and empirical state occupancy into a unified Fokker–Planck diagnostic framework. It provides the quantitative ingredients needed to compare closed stationary predictions with empirical repertoire organization and to identify abundance regions in which the observed system behaves as if it were open, nonstationary, or otherwise inconsistent with a closed one-dimensional diffusion process.

# Full-pipeline bootstrap of FP flux and source 

The script `9-boot_full_pipeline_flux.py` implements a **full uncertainty propagation framework** for Fokker–Planck diagnostics by combining:


- drift refitting
- diffusion re-estimation
- empirical density reconstruction
- flux and source computation

within a **cluster bootstrap scheme**.

Unlike previous steps (point estimates), this script quantifies **confidence intervals for $J(x)$ and $S(x)$** by resampling the entire pipeline.

---

## Conceptual pipeline

Each bootstrap replicate performs:

$$
\text{transitions}
\rightarrow b(x)
\rightarrow D(x)
\rightarrow p(x)
\rightarrow J(x), S(x)
$$

This ensures that **all sources of uncertainty are propagated jointly**.

---

## Bootstrap scheme

### Cluster resampling

Data are resampled at the clonotype level:

$$
(\text{subject}, \text{aaSeqCDR3})
$$

This preserves temporal structure and avoids pseudoreplication.

---

## Step-by-step algorithm

For each bootstrap replicate $r$:

### 1. Resample transitions

$$
\{(x_0, \Delta x)\}_r
$$

### 2. Refit drift

Estimate:

$$
b_r(x)
$$

using:
- `OU_linear`
- `hinge_plateau_smooth`

---

### 3. Estimate diffusion

Two modes:

#### Default (recommended)

$$
r_i = \Delta x_i - b_r(x_{0,i})
$$

$$
D_r(x) = \frac{\mathrm{Var}(r \mid x)}{2dt}
$$

#### Optional

Use precomputed residuals.

---

### 4. Resample trajectories to reconstruct empirical density

$$
p_r(x)
$$

Supports:

- ignore
- filter (observable only)
- weight (weighted density)

---

### 5. Compute flux

$$
J_r(x) = b_r(x)p_r(x) - \frac{d}{dx}[D_r(x)p_r(x)]
$$

---

### 6. Compute source

$$
S_r(x) = \frac{dJ_r(x)}{dx}
$$

---

## Diffusion estimation

Residual-based.

### Robust estimator

$$
\mathrm{Var} \approx (1.4826 \cdot \mathrm{MAD})^2
$$

### Alternative

$$
\mathrm{Var} = \mathbb{E}[r^2]
$$

---

## Empirical density

Estimated via histogram on a fixed grid:

$$
p(x) = \frac{\text{counts}}{N \cdot \Delta x}
$$

Optional Gaussian smoothing is applied.

---

## Bulk region

To avoid edge artifacts, metrics are computed on:

- trimmed grid (remove extremes)
- density-threshold filter

---

## Global metrics

Computed on each bootstrap:

- $\int |J(x)|\,dx$
- $\int J(x)^2\,dx$
- $\int J(x)\,dx$
- $\int |S(x)|\,dx$
- $\int S(x)^2\,dx$

---

## Confidence intervals

For each quantity:

$$
\text{mean},\quad \text{CI}_{\alpha}
$$

computed across bootstrap replicates.

---

## Outputs

### Curve summaries

- `J_boot_summary.csv`
- `S_boot_summary.csv`

### Global metrics

- `global_boot_summary.csv`

### Drift and diffusion summaries

- `drift_boot_summary.csv`
- `diffusion_boot_summary.csv`

### Raw bootstrap arrays

- `J_boot_reps.npy`
- `S_boot_reps.npy`

### Figures

- `flux_J_with_CI.png`
- `source_S_with_CI.png`

### Report

- `report.txt`

---

## Observable handling

Three modes:

- ignore -> use all data
- filter -> keep observable only
- weight -> weighted density

---

## Interpretation

This step answers:

> Are flux and source significantly different from zero?

Key diagnostic:

- If the CI of $\int J(x)\,dx$ includes 0, the result is consistent with a closed system.
- Otherwise, it supports evidence for open dynamics.

---

## Advantages

- full uncertainty propagation
- model comparison under noise
- robustness to sampling variability

---

## Limitations

- computationally intensive
- depends on upstream model choices
- bootstrap assumes cluster independence

---

## Role in workflow

Final validation step:

$$
\text{Drift + Diffusion}
\rightarrow \text{FP diagnostics}
\rightarrow \text{Bootstrap validation}
$$

---

## Conclusion

This script provides a **rigorous statistical validation layer** for Fokker-Planck diagnostics, allowing quantitative assessment of whether observed repertoire dynamics are compatible with a closed stochastic system or require additional open-process mechanisms.