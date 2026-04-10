# Detailed Description of `1-powerlaw_qc_analysis.py`

This document describes in detail the logic, mathematical assumptions, command-line arguments, internal workflow, and outputs of the script `1-powerlaw_qc_analysis.py`. 

## 1. General Purpose of the Script

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

## 2. Expected Input Structure

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

## 3. High-Level Workflow

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

## 4. Utility Functions

### 4.1 `ensure_dir`

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

### 4.2 `robust_mad`

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

### 4.3 `robust_z`

This function defines a robust z-score:

$$
z_i = \frac{x_i - \operatorname{median}(x)}{\operatorname{MAD}(x)}
$$

If the MAD is not finite or is numerically too small, the function returns `NaN` for all entries.

This robust standardization is used because gamma discordance may not be well modeled by a Gaussian distribution and may contain outliers.

### 4.4 `depth_ratio`

Given two replicate sequencing depths $d_1$ and $d_2$, the function returns:

$$
\text{depth\_ratio}(d_1, d_2) = \frac{\max(d_1, d_2)}{\min(d_1, d_2)}
$$

provided that both are strictly positive. Otherwise, it returns `NaN`.

This metric is only diagnostic in this script. It does **not** enter the final quality-control rule.

---

## 5. Mathematical Basis of the Power-Law Fit

The script uses a **continuous Clauset-style power-law fitting procedure** on clonotype frequencies. Frequencies are treated as positive real values.

### 5.1 Input Frequencies

For each repertoire, the script converts raw counts into relative frequencies:

$$
f_i = \frac{c_i}{\sum_j c_j}
$$

where:

- $c_i$ is the read count of clonotype $i$
- $\sum_j c_j$ is the total repertoire depth

Only strictly positive finite frequencies are kept.

### 5.2 Validation of Positive Values

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

### 5.3 Continuous Pareto Model

The fitted power-law density on the tail $x \ge x_{\min}$ is:

$$
p(x \mid \alpha, x_{\min}) = (\alpha - 1) x_{\min}^{\alpha - 1} x^{-\alpha}
$$

valid for:

$$
\alpha > 1, \qquad x_{\min} > 0
$$

The script denotes the fitted exponent as both `alpha_hat` and `gamma_hat`. In practice, for the purpose of this workflow, gamma is simply the fitted power-law exponent.

### 5.4 Maximum-Likelihood Estimator of the Exponent

Given the tail sample $x_1, \dots, x_n$ with all $x_i \ge x_{\min}$, the continuous maximum-likelihood estimator is:

$$
\hat{\alpha} = 1 + \frac{n}{\sum_{i=1}^{n} \log(x_i / x_{\min})}
$$

This is implemented in `continuous_alpha_mle`.

If the denominator is non-positive, the estimate is considered invalid and the function returns `NaN`.

### 5.5 Model CDF

The theoretical cumulative distribution function of the continuous Pareto law is:

$$
F(x) = 1 - \left(\frac{x}{x_{\min}}\right)^{1 - \alpha}
$$

for $x \ge x_{\min}$.

This is implemented by `powerlaw_cdf_continuous`.

### 5.6 Kolmogorov-Smirnov Distance

For a sorted tail sample $x_{(1)} \le \dots \le x_{(n)}$, the empirical CDF is:

$$
F_{\text{emp}}(x_{(k)}) = \frac{k}{n}
$$

The model CDF is $F_{\text{model}}(x_{(k)})$. The KS distance is:

$$
D_{\text{KS}} = \max_k \left| F_{\text{emp}}(x_{(k)}) - F_{\text{model}}(x_{(k)}) \right|
$$

This statistic is used to choose the best tail threshold $x_{\min}$.

### 5.7 KS-Based Estimation of the Tail Threshold

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

### 5.8 Final Per-Repertoire Power-Law Fit

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

## 6. Tail-Model Comparison: Power Law Versus Truncated Log-Normal

The second major component of Step A1 is a model comparison restricted to the high-frequency tail.

Importantly, the script does **not** compare models using model-specific thresholds. Instead, it defines tail regions by **frequency quantiles**. This is methodologically important because it enforces comparison on the **same support** for both candidate models.

### 6.1 Threshold Quantiles

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

### 6.2 Power-Law Log-Likelihood on the Tail

For a tail sample $x_1, \dots, x_n$ and fitted exponent $\alpha$, the continuous Pareto log-likelihood is:

$$
\ell_{\text{PL}} = n \log(\alpha - 1) + n(\alpha - 1) \log(x_{\min}) - \alpha \sum_{i=1}^{n} \log(x_i)
$$

This is implemented in `powerlaw_tail_loglik`.

### 6.3 Truncated Log-Normal Model

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

### 6.4 Delta Log-Likelihood

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

### 6.5 Threshold Scan Output

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

## 7. File Reading, Pair Indexing, and Replicate Processing

### 7.1 `read_rep_file`

Each file is read and restricted to the columns `aaSeqCDR3` and `readCount`.

The function:

- coerces `readCount` to numeric
- removes non-finite entries
- removes non-positive entries

Thus, every retained clonotype satisfies:

$$
\text{readCount}_i > 0
$$

### 7.2 `index_pairs`

This function scans the input directory, matches filenames to the user-provided regex pattern, and stores files indexed by `(subject, time, replica)`.

A pair is considered complete only if **both** replicate 1 and replicate 2 exist for the same `(subject, time)`.

The output is a list of tuples:

```text
(subject, time, file_rep1, file_rep2)
```

### 7.3 `fit_one_replicate`

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

## 8. Step A2: Gamma-Only Quality Control

The second major output of the script is a QC classification for each replicate pair based **only** on gamma concordance.

This is intentionally restrictive. Other quantities such as `xmin`, KS distance, depth ratio, tail size, or tail fraction are retained as diagnostics but do not affect PASS or FAIL.

### 8.1 Pairwise Discordance in Gamma

For each sample pair, the script extracts:

- $\gamma_1$: fitted exponent in replicate 1
- $\gamma_2$: fitted exponent in replicate 2

Then it computes:

$$
\Delta \gamma = |\gamma_1 - \gamma_2|
$$

This is the central QC variable.

### 8.2 Pairwise Mean Gamma

A descriptive quantity also recorded is:

$$
\gamma_{\text{mean}} = \operatorname{mean}(\gamma_1, \gamma_2)
$$

computed as a NaN-aware mean.

### 8.3 Robust Standardization of Gamma Discordance

Across all replicate pairs, the script computes a robust z-score for $\Delta \gamma$:

$$
z_{\Delta \gamma} = \frac{\Delta \gamma - \operatorname{median}(\Delta \gamma)}{\operatorname{MAD}(\Delta \gamma)}
$$

This score is not a classical Gaussian z-score. It is a median-and-MAD standardized outlier score.

### 8.4 QC Threshold

The user provides a cutoff `gamma_mad_cutoff`, with default:

$$
\text{gamma\_mad\_cutoff} = 3.5
$$

A pair is marked as a gamma outlier if:

$$
|z_{\Delta \gamma}| > \text{gamma\_mad\_cutoff}
$$

### 8.5 PASS/FAIL Rule

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

### 8.6 QC Reason Field

For each failing pair, the script composes a reason string such as:

- `fail:fit_failed`
- `fail:delta_gamma`
- `fail:fit_failed,delta_gamma`

This makes the failure mechanism explicit in the output table.

---

## 9. Aggregate Summaries

After computing per-replicate and per-pair results, the script generates summary tables.

### 9.1 Overall QC Summary

The overall summary includes:

- number of pairs
- number of PASS
- number of FAIL
- median `delta_gamma`
- median `delta_xmin`
- median `depth_ratio`
- median `gamma_mean`

These are descriptive summaries of the full dataset.

### 9.2 Subject-Level QC Summary

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

## 10. Detailed Explanation of Command-Line Arguments

The script exposes the following CLI arguments.

### 10.1 Input and Output Arguments

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

### 10.2 Power-Law Fit Arguments

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

### 10.3 Threshold-Scan Arguments

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

### 10.4 QC Argument

#### `--gamma_mad_cutoff`

Robust outlier cutoff applied to $z_{\Delta \gamma}$.

Default:

```text
3.5
```

This controls the stringency of replicate concordance QC.

Larger values make the QC more permissive; smaller values make it stricter.

---

## 11. Detailed Description of Output Files

The script writes nine output tables.

### 11.1 `A1_repertoire_tail_fits.csv`

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

### 11.2 `A1_repertoire_threshold_scan.csv`

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

### 11.3 `A1_repertoire_model_compare.csv`

This file aggregates the threshold scan at the level of each replicate.

It summarizes, for each repertoire:

- number of threshold levels analyzed
- mean and median normalized delta log-likelihood
- fraction of thresholds where power law is preferred
- median tail fraction

This answers whether model preference is stable across thresholds within a repertoire.

### 11.4 `A1_pair_tail_preference_summary.csv`

This file aggregates the threshold-scan comparison at the replicate-pair level.

It summarizes:

- mean of the replicate-level median delta log-likelihood
- mean fraction of thresholds favoring power law
- minimum fraction of thresholds favoring power law
- number of replicates contributing to the pair

This table is useful when model preference must be summarized at the sample level rather than at the individual replicate level.

### 11.5 `A2_replicate_powerlaw_fits.csv`

This file stores the raw replicate-level power-law fit outputs used by the QC step.

It is essentially the pre-QC per-replicate fit table.

### 11.6 `A2_pair_gamma_qc.csv`

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

### 11.7 `A2_gamma_qc_thresholds.csv`

This file records the QC rule and the numerical threshold used.

It documents the cutoff and the precise PASS/FAIL logic, improving reproducibility.

### 11.8 `A2_summary_overall.csv`

This is a single-row summary of the full dataset after QC classification.

### 11.9 `A2_summary_by_subject.csv`

This file provides the same type of descriptive summary stratified by subject.

---

## 12. Conceptual Interpretation of Step A1

Step A1 is not merely a fitting exercise. It is a methodological justification stage.

Its purpose is to determine whether the high-frequency clone-size distribution is compatible with a heavy tail and whether a power-law tail provides a better account than a truncated log-normal on matched support regions.

Two design choices are especially important.

### 12.1 Quantile-Based Tail Thresholds

By scanning quantiles rather than fitting separate model-specific thresholds, the script ensures that both candidate models are evaluated on the same tail subset. This avoids a common methodological ambiguity in tail-model comparison.

### 12.2 Per-Clonotype Normalized Delta Log-Likelihood

The quantity:

$$
\overline{\Delta \ell} = \frac{\ell_{\text{PL}} - \ell_{\text{TLN}}}{n_{\text{tail}}}
$$

allows comparison across repertoires and thresholds even when the number of tail clonotypes differs.

This is particularly important when the tail becomes progressively smaller at high quantiles.

---

## 13. Conceptual Interpretation of Step A2

Step A2 uses replicate concordance in gamma as a technical reproducibility metric.

The rationale is that if two technical replicates are measuring the same underlying repertoire, then the inferred power-law exponent of the high-frequency tail should be relatively stable across replicates. Large discrepancies in gamma may indicate technical instability, insufficient sampling, fitting failure, or other sources of unreliability.

The use of a robust z-score rather than a standard z-score makes the QC rule resistant to skewness and extreme values in the empirical distribution of $\Delta \gamma$.

At the same time, the choice to ignore other diagnostics in the final PASS/FAIL rule makes the procedure deliberately transparent and conservative in scope.

---

## 14. What the Script Does Not Do

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

## 15. Computational Logic of the Main Loop

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

## 16. Possible Failure Modes

Several types of failure are handled explicitly.

### 16.1 Missing Required Columns

If a repertoire file lacks either `aaSeqCDR3` or `readCount`, the script raises an error.

### 16.2 Zero or Non-Positive Depth

If the sum of counts is zero, the replicate fit fails.

### 16.3 Too Few Positive Frequencies

If the number of positive frequencies is smaller than `min_tail_size`, the power-law fit fails.

### 16.4 No Valid KS-Based Threshold

If none of the candidate thresholds yields a valid power-law tail with sufficient size and valid exponent, the fit fails.

### 16.5 Too Few Tail Observations in Threshold Scan

If a quantile-defined tail contains fewer than `scan_min_tail_size` observations, that threshold is skipped.

These failure modes are propagated transparently into the output tables.

---

## 17. Reproducibility and Transparency

A strength of the script is that it records both:

- the raw fit outputs
- the exact QC decision rule

This makes it possible to audit all decisions post hoc.

In particular, because the QC rule depends only on fit success and robust gamma discordance, one can later assess the influence of excluded samples without ambiguity.

---

## 18. Recommended Interpretation in a Methods Section

In manuscript terms, the script supports the following methodological narrative.

First, clone-size distributions are examined for heavy-tailed structure and model preference in the high-frequency regime, using quantile-defined tails and matched-support comparison between power-law and truncated log-normal models.

Second, technical reproducibility is assessed by fitting the power-law exponent separately in each technical replicate and quantifying replicate discordance in gamma. Samples are classified as PASS or FAIL using a robust outlier rule based exclusively on delta gamma, while additional fit diagnostics are retained for descriptive purposes.

This separation between structural validation and technical QC is a central conceptual feature of the script.

---

## 19. Summary

In summary, `1-powerlaw_qc_analysis.py` performs two tightly connected but conceptually distinct tasks. Step A1 evaluates whether repertoire clone-size distributions exhibit a robust heavy-tailed structure and whether a power-law tail is favored over a truncated log-normal across quantile-defined tail thresholds. Step A2 then uses replicate concordance of the inferred power-law exponent gamma as a transparent and robust technical quality-control criterion.

The script is therefore best understood as a diagnostic and analytical entry point for downstream repertoire analysis. It validates the use of a heavy-tail representation, documents model preference across tail thresholds, and identifies replicate pairs whose tail exponent estimates are insufficiently concordant for confident downstream use.
