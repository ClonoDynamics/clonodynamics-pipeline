# Model-agnostic diagnostics of clonotype frequency dynamics

To characterize clonotype frequency dynamics without imposing an explicit parametric model, we analyzed a table of lag-sampled transitions derived from longitudinal repertoire sequencing data. Each transition corresponds to the change in clonotype frequency between two observation times separated by a lag $\Delta t$. The analysis is explicitly model-agnostic: it does not fit a stochastic model but instead computes empirical summary statistics that describe directional tendencies and fluctuation amplitudes across temporal and abundance scales.

The transition table must contain at minimum the variables

subject, aaSeqCDR3, dt, x0, x1, dx, obs_class

where `subject` identifies the individual, `aaSeqCDR3` the clonotype sequence, `dt` the time lag between observations, `x0` and `x1` the log-frequencies at the initial and final time points, `dx = x1 - x0` the log-frequency increment, and `obs_class` the observation class used for optional filtering. Additional columns may include a transition weight and metadata describing observation paths. Only rows with finite numeric values for the variables required in each calculation were retained.

---

# Log-frequency convention

Log-frequencies may be provided either in natural logarithm or base-10 logarithm. The convention is controlled by the parameter `x_base`.

If log-frequencies are expressed in natural logarithm,

$$
\log_{10} f_0 = \frac{x_0}{\ln 10},
\qquad
\log_{10} f_1 = \frac{x_1}{\ln 10}.
$$

If they are already in base-10 logarithm,

$$
\log_{10} f_0 = x_0,
\qquad
\log_{10} f_1 = x_1.
$$

A geometric abundance coordinate was defined as the midpoint between the two log-frequencies,

$$
\log_{10} f_{\mathrm{geo}}
=
\frac{1}{2}
(\log_{10} f_0 + \log_{10} f_1).
$$

The variable $\log_{10} f_0$ was used for drift diagnostics, whereas $\log_{10} f_{\mathrm{geo}}$ was used for fluctuation diagnostics.

---

# Filtering

Transitions may be filtered by observation class and by time lag.

Observation-class filtering retains only rows satisfying

$$
\mathrm{obs\_class} = c
$$

where $c$ is specified by the user.

Time lag filtering retains rows satisfying

$$
dt > dt_{\min}, \qquad dt \le dt_{\max}.
$$

After filtering, numeric columns were coerced to finite values and `dt` was converted to integer.

---

# Optional subsampling

To control computational cost and partially balance sample sizes across lags, the script allows optional subsampling. When enabled, up to a fixed number of rows are sampled without replacement independently within each lag $\Delta t$. Subsampling affects the analyzed dataset itself rather than only downstream summaries.

---

# Transition weights

The analysis optionally incorporates per-transition weights.

If the parameter `weight_col` is set to `"1"`, all observations receive equal weight. Otherwise the selected column is converted into analysis weights

$$
w_i^\ast =
\begin{cases}
w_i & \text{if } w_i>0 \text{ and finite} \\
0 & \text{otherwise}
\end{cases}
$$

Rows with zero weight remain in the dataset but do not contribute numerically to weighted summaries.

---

# Weighted summary statistics

## Weighted mean

For a variable $z_i$ the weighted mean is

$$
\bar z_w =
\frac{\sum_i w_i^\ast z_i}{\sum_i w_i^\ast}.
$$

## Weighted probability of positive increment

Directional asymmetry of the displacement distribution is quantified as

$$
P_w(\Delta x > 0)
=
\frac{\sum_i w_i^\ast \mathbf{1}(\Delta x_i > 0)}{\sum_i w_i^\ast}.
$$

where $\mathbf{1}(\cdot)$ is the indicator function.

## Weighted median

For values $z_i$ with weights $w_i^\ast$, the weighted median $\tilde z_w$ is defined as the smallest value such that the cumulative ordered weight reaches at least half of the total weight,

$$
\tilde z_w = Q_w(0.5).
$$

More generally, weighted quantiles are obtained by sorting observations by value and identifying the point where the cumulative weight exceeds the target quantile level.

---

# Quantile binning

To obtain abundance-resolved diagnostics, transitions were stratified using quantile bins.

For an abundance variable $a_i$, quantile edges are defined by

$$
q_j = \frac{j}{K}, \quad j=0,\dots,K
$$

where $K$ is the number of bins. Empirical quantiles at these probabilities define the bin boundaries.

Duplicate edges are removed to maintain valid intervals even when the abundance distribution contains repeated values.

Each bin is described by

- lower bound
- upper bound
- bin center
- bin identifier

For initial-frequency bins

$$
x_{\mathrm{lo}},\; x_{\mathrm{hi}},\; x_{\mathrm{center}}
$$

and for geometric-frequency bins

$$
f_{\mathrm{lo}},\; f_{\mathrm{hi}},\; f_{\mathrm{center}}.
$$

The bin center corresponds to the midpoint between interval edges.

---

# Counts table

Before computing drift and fluctuation diagnostics, the script generates a counts table summarizing sample support. Counts are reported

- per lag $\Delta t$
- per lag and initial-frequency bin
- per lag and geometric-frequency bin

This table is purely descriptive and does not use weights or bootstrap estimation.

---

# Drift diagnostics

Transitions are grouped jointly by lag $\Delta t$ and by quantile bin of $\log_{10} f_0$.

Within each group two quantities are computed.

## Median increment

$$
m(\Delta t,x\text{-bin})
=
\mathrm{median}_w
\left(
\Delta x
\mid
\Delta t,\log_{10} f_0 \in x\text{-bin}
\right)
$$

## Probability of positive displacement

$$
p_+(\Delta t,x\text{-bin})
=
P_w
\left(
\Delta x > 0
\mid
\Delta t,\log_{10} f_0 \in x\text{-bin}
\right)
$$

Bins containing fewer than `min_rows_per_bin` transitions are retained in the output but assigned missing summary values.

---

# Mean squared displacement

Fluctuation amplitude across time lags is summarized by the weighted mean squared displacement

$$
\mathrm{MSD}(\Delta t)
=
E_w[(\Delta x)^2 \mid \Delta t]
=
\frac{\sum_i w_i^\ast (\Delta x_i)^2}{\sum_i w_i^\ast}.
$$

---

# Frequency-resolved MSD

To test whether fluctuation amplitude depends on abundance, MSD is also computed after stratifying transitions by $\log_{10} f_{\mathrm{geo}}$.

$$
\mathrm{MSD}(\Delta t,f\text{-bin})
=
E_w[(\Delta x)^2
\mid
\Delta t,\log_{10} f_{\mathrm{geo}}\in f\text{-bin}]
$$

Groups with insufficient observations are retained but assigned missing MSD estimates.

---

# Bootstrap confidence intervals

Uncertainty is estimated using nonparametric bootstrap resampling performed independently within each group.

For a group containing $n$ observations, bootstrap samples are generated by sampling $n$ indices with replacement.

Let $\hat\theta$ denote the observed statistic and

$$
\hat\theta^{(1)},\dots,\hat\theta^{(B)}
$$

the bootstrap replicates. Percentile confidence intervals are computed as

$$
\theta_{\mathrm{lo}} =
Q_{0.025}(\hat\theta^{(1)},\dots,\hat\theta^{(B)})
$$

$$
\theta_{\mathrm{hi}} =
Q_{0.975}(\hat\theta^{(1)},\dots,\hat\theta^{(B)})
$$

Bootstrap intervals are computed only when the group contains at least 50 observations.

---

# Output tables

The script generates four main output tables.

## dt_freqbin_counts.csv

Counts of transitions by lag and abundance bins.

## drift_by_dt_xbins.csv

Drift diagnostics including

- lag
- abundance bin
- sample size
- weighted median increment
- bootstrap confidence interval
- probability of positive increment
- bootstrap confidence interval

## msd_by_dt.csv

Lag-level mean squared displacement estimates and bootstrap intervals.

## msd_by_dt_freqbins.csv

Frequency-resolved MSD estimates by lag and geometric-frequency bin.

---

# Methodological rationale

This diagnostic framework provides a purely empirical characterization of clonotype frequency dynamics. Drift statistics summarize the central tendency and directional asymmetry of log-frequency increments as functions of abundance and time lag, without assuming a particular functional form of the drift term. Mean squared displacement quantifies fluctuation amplitude across temporal scales and allows testing whether stochastic variability depends on abundance regime. Weighted estimators allow transition-specific reliability to be incorporated when available, while bootstrap intervals provide a nonparametric measure of uncertainty. Together these diagnostics provide a model-independent empirical description of repertoire dynamics that can guide the formulation and validation of downstream stochastic models.