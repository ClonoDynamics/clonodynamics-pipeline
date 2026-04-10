# Noise Model Pipeline for RepSeq Data
## Adaptation of the NoisET Framework of Mora et al.

## 1. Overview

The script `noise_pipeline.py` implements a streamlined end-to-end pipeline for technical noise modelling, empirical denoising, and generation of trajectory-ready observability outputs from longitudinal immune repertoire sequencing data.

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

## 2. General Statistical Rationale

Immune repertoire sequencing data are affected by multiple sources of variability, including sampling noise, library preparation effects, amplification noise, and sequencing stochasticity. As a consequence, the observed count of a clonotype is not a direct measurement of its latent biological abundance, but a noisy realization of an underlying frequency.

The central idea of the pipeline is that technical replicate pairs contain direct information about this experimental noise. By fitting a probabilistic model to replicate counts, one can estimate the expected variability under a null hypothesis of purely technical fluctuation. This null model can then be used to assign to each clonotype a probability of being compatible with technical noise, and therefore to distinguish detectable biological signal from stochastic background.

This logic follows the core philosophy of NoisET: infer a probabilistic noise model from technical replicates, then use that model to calibrate the detectability of clonotypes.

---

## 3. Inputs

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

## 4. External QC Filtering

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

## 5. Pairwise Data Representation

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

## 6. Core Noise Model: `noiseK_nb.py`

The inferential core of the pipeline is the function

`fit_noiseK_nb_powerlaw(...)`

implemented in `noiseK_nb.py`.

This function fits a probabilistic model in which:

1. each clonotype has an unknown latent frequency $f$  
2. latent frequencies are drawn from a power-law prior  
3. observed replicate counts arise from a Negative Binomial observation process conditional on $f$  

This is the essential NoisET-style component that transforms replicate counts into a likelihood-based model of technical noise.

---

## 7. Power-Law Prior on Latent Frequencies

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

## 8. Negative Binomial Observation Model

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

## 9. Marginal Likelihood of a Clonotype

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

## 10. Pairwise Log-Likelihood and Parameter Estimation

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

## 11. Interpretation of the Fitted Parameters

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

## 12. Per-Clonotype Log-Probability Tables

After model fitting, the pipeline stores for each QC-pass pair a per-clonotype table containing at least:

- clonotype identity
- replicate counts
- replicate frequencies
- geometric-mean frequency
- pair metadata
- per-clonotype `logP`

These files are written under the `per_clone_long_logP` directory and serve as the basis for all downstream empirical denoising steps.

---

## 13. Null Model Choices in the Pipeline

Once the pairwise `logP` values have been computed, the pipeline constructs an empirical null distribution by pooling them according to one of three choices.

### 13.1 Clone-level null

Each pair defines its own pool:

$$
\mathrm{pool\_id} = (subject, time).
$$

This is the most local null and preserves pair-specific noise structure.

### 13.2 Subject-level null

All QC-pass pairs from the same subject are pooled:

$$
\mathrm{pool\_id} = subject.
$$

However, this is used only if the subject has at least a specified number of QC-pass pairs. If not, the pipeline falls back to the global pool.

### 13.3 Global null

All QC-pass pairs across all subjects and times are pooled together:

$$
\mathrm{pool\_id} = \mathrm{GLOBAL}.
$$

This provides the largest empirical null but may average over subject-specific noise features.

---

## 14. Quantile-Based Observability Threshold

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

## 15. Empirical Denoising Through Tail Probabilities

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

## 16. Alpha-Sensitivity Analysis

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

## 17. Aggregation of Noise Parameters

For reporting purposes, the fitted noise parameters can be aggregated across QC-pass successful pairs.

If the null choice is `subject`, aggregation is performed per subject using medians across pairs.

Otherwise, a global median summary is produced.

This yields compact descriptors of the experimental noise landscape after QC filtering.

---

## 18. Trajectory-Ready Outputs

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

## 19. Relation to NoisET

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

## 20. Conceptual Role in the Full Analysis Pipeline

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

## 21. Interpretation

The central output of the pipeline is not merely a filtered repertoire, but a probabilistically calibrated representation of clonotype observability.

Each clonotype is evaluated with respect to a noise model inferred from technical replicates. This makes it possible to separate two distinct components of observed variability:

- fluctuations that are statistically compatible with technical noise
- fluctuations that remain detectable after calibration against that noise

In this sense, the pipeline provides the statistical foundation for subsequent inference on clonotype turnover, persistence, and stochastic dynamics over time.

---

## 22. Conclusion

The script `noise_pipeline.py`, together with the core inferential module `noiseK_nb.py`, implements a complete noise-aware framework for longitudinal RepSeq preprocessing.

Its essential components are:

- QC-based selection of reliable replicate pairs
- likelihood-based inference of technical noise using a Negative Binomial model with power-law prior
- extraction of per-clonotype marginal log-probabilities
- empirical denoising through pooled null distributions
- generation of trajectory-ready observability outputs

Because the model is explicitly adapted from the NoisET framework of Mora et al., the pipeline retains a principled probabilistic interpretation while being directly usable for longitudinal clonotype-dynamics studies.