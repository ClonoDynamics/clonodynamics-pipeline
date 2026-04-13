# Construction of Clonotype Transition Patterns

## 1. Purpose of the Script

The script `3-build_transition_patterns.py` constructs a dataset of clonotype frequency transitions over time from longitudinal repertoire data.

Its purpose is to convert a trajectory-based representation of clonotype dynamics into a transition-based representation suitable for stochastic inference, including the estimation of drift and diffusion functions.

---

## 2. Input Data Structure

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

## 3. Definition of the Dynamical Variable

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

## 4. Transition Definition

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

## 5. Temporal Constraint

The parameter `max_dt` limits the maximum time lag:

$$
\Delta t \leq \Delta t_{\max}
$$

- Default: 6  
- If `max_dt <= 0`, no filtering is applied  

This parameter controls the temporal scale of the transitions included in the analysis.

---

## 6. Transition Patterns

By default, only the ALL transition pattern is constructed. The ADJ pattern is generated only if explicitly requested at runtime using the `--build_adj` flag.

### 6.1 All pairs (all)

All ordered pairs of timepoints are included:

$$
\forall i < j
$$

The number of transitions per trajectory with \( n \) timepoints is:

$$
N = \binom{n}{2}
$$

This pattern captures both short- and long-range temporal dynamics.

---

### 6.2 Consecutive observed pairs (adj, optional)

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

## 7. Observation Classes

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

## 8. Weighting System

If enabled, each transition is assigned a weight:

$$
w = w_{\text{clone}} \cdot w_{\Delta t} \cdot w_{\text{class}}
$$

This weighting scheme corrects structural imbalances in the dataset.

---

## 9. Clone Weight

The clone-level weight is:

$$
w_{\text{clone}} = \frac{1}{N_{\text{transitions per clone}}}
$$

This prevents trajectories with many timepoints from dominating the dataset.

---

## 10. Time-Lag Weight

If the option `balance_dt` is enabled:

$$
w_{\Delta t} = \frac{1}{N(\Delta t)}
$$

Otherwise:

$$
w_{\Delta t} = 1
$$

This compensates for uneven sampling across time intervals.

---

## 11. Class Weights

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

---

## 12. Weight Normalization

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

## 13. Pattern-Specific Behavior

- The ALL pattern supports full weighting, including time-lag balancing  
- The ADJ pattern typically does not apply time-lag weighting  

---


**Default behavior**

- The ALL pattern is always generated  
- The ADJ pattern is generated only if `--build_adj` is specified  

This design reflects the primary use of ALL transitions for dynamical inference, while ADJ transitions are used for validation or sensitivity analyses.


## 14. Output Files

The script produces by default:

- transitions_all.csv  

If `--build_adj` is specified, it additionally produces:

- transitions_adj.csv  

Each row contains:

$$
(t_0, t_1, \Delta t, x_0, x_1, \Delta x)
$$

as well as metadata such as observation class, pattern type, and optional weights.

---

## 15. Count Tables

The script also outputs aggregated counts (for each generated pattern):

$$
N(\text{subject}, \text{obs\_class}, \Delta t)
$$

These are used for diagnostic and stratification purposes.

---

## 16. Stochastic Interpretation

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

## 17. Statistical Dependencies and Bias Control

The construction of transition datasets introduces statistical dependencies and sampling biases.

### 17.1 Dependence in the ALL pattern

Multiple transitions may share the same underlying observations, for example:

$$
(t_0 \to t_2), \quad (t_0 \to t_3), \quad (t_1 \to t_3)
$$

This induces correlation among transitions and violates independence assumptions.

---

### 17.2 ADJ pattern

Restricting to consecutive transitions reduces overlap but does not eliminate dependence between transitions from the same clonotype.

---

### 17.3 Sampling imbalance

The dataset is intrinsically heterogeneous:

$$
P(\Delta t), \quad P(\text{clone}), \quad P(\text{obs\_class})
$$

are not uniform.

---

### 17.4 Role of weighting

The weighting scheme:

$$
w = w_{\text{clone}} \cdot w_{\Delta t} \cdot w_{\text{class}}
$$

reduces bias but does not restore statistical independence.

---

### 17.5 Implications for inference

Estimated quantities such as drift and diffusion should be interpreted as effective estimators under dependence rather than independent samples.

---

### 17.6 Cluster bootstrap

Reliable uncertainty estimation requires resampling at the level of:

$$
(\text{subject}, \text{clonotype})
$$

to preserve the dependence structure.

---

### 17.7 Practical consequence

Valid inference requires explicit handling of dependence and appropriate resampling strategies.

---

## 18. Role in the Analysis Pipeline

The script performs the transformation:

$$
\text{trajectories} \rightarrow \text{transitions} \rightarrow \text{stochastic inference}
$$

It is a central step connecting denoised trajectories to dynamic modeling.

---

## 19. Conclusion

The script transforms longitudinal clonotype data into a transition-based representation:

$$
\{x(t)\} \;\longrightarrow\; \{(x_0, \Delta x, \Delta t, w)\}
$$

This provides the minimal sufficient structure for quantitative modeling of clonotype dynamics in a stochastic framework.