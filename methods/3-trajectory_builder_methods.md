# Supplementary Methods 

— Construction of longitudinal clonotype trajectories (`3-trajectory_builder.py`)

## Overview

This script builds a standardized long-format trajectory table from a denoised per-clonotype input table. Its role in the workflow is purely structural: it converts one or more rows per clonotype, subject, and time point into a single longitudinal representation that can be used downstream for transition construction and dynamical inference.

It performs only formatting and aggregation:
- selection of one frequency column via `--freq-col`
- retention of the required identifiers `subject`, `time`, and `aaSeqCDR3`
- aggregation of duplicate rows at the same subject–clonotype–time combination
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

\[
\mathrm{log\_freq}_i = \ln(f_i + \varepsilon),
\]

where \(f_i\) is the selected frequency value for row \(i\), and \(\varepsilon > 0\) is the user-defined constant.

## Processing steps

### 1. Input reading and validation

The input CSV is loaded into a pandas DataFrame. The script checks for the presence of the required identifiers

\[
\{\texttt{subject},\ \texttt{time},\ \texttt{aaSeqCDR3}\}.
\]

It also checks that the column specified by `--freq-col` exists. If not, execution stops.

### 2. Type coercion of identifiers

The script standardizes the key identifiers as follows:

- `subject` is coerced to numeric and cast to integer
- `time` is coerced to numeric
- `aaSeqCDR3` is cast to string

This prevents mixed internal representations such as `"1"` and `1` from being treated as different subjects.

### 3. Frequency validity filtering

The selected frequency column is converted to numeric with coercion of invalid entries to missing values. Rows are retained only when the selected frequency is finite and non-negative.

If \(f_i\) denotes the selected raw frequency value for row \(i\), the row is retained only if

\[
f_i \in \mathbb{R}, \qquad f_i \ge 0, \qquad |f_i| < \infty.
\]

Rows with missing, non-numeric, infinite, or negative frequency values are discarded. No additional abundance filter is applied. Zero values are allowed.

This is the only direct row-level filtering step. It is a mathematical validity check, not a biological selection rule.

### 4. Standardization of the abundance column

After filtering, the selected frequency values are copied into a standardized output column named `freq`. This makes the downstream schema independent of the original frequency column name.

### 5. Treatment of `observable`

If `observable` exists in the input, the script converts it robustly to Boolean. String values such as `true`, `t`, `1`, `yes`, and `y` are interpreted as `True`; all other values are interpreted as `False`. If the input column is already Boolean, it is preserved.

If `observable` is absent, the script creates it and assigns `False` to all rows.

This guarantees that the output trajectory table always contains an `observable` column, even if observability has not yet been propagated from upstream steps.

### 6. Aggregation over duplicate subject–clone–time entries

For each unique triplet

\[
(subject,\ aaSeqCDR3,\ time),
\]

the script keeps one row only. If duplicate rows exist, it aggregates them as:

- `freq` = arithmetic mean
- `observable` = maximum

Suppose a given subject–clone–time triplet has \(m\) rows with frequencies \(f_1, \dots, f_m\) and Boolean observability indicators \(o_1, \dots, o_m\), where each \(o_j \in \{0,1\}\). The aggregated values are

\[
\bar f = \frac{1}{m}\sum_{j=1}^{m} f_j,
\]

and

\[
o^* = \max_{1 \le j \le m} o_j.
\]

Thus the aggregated row is considered observable if any duplicated row is observable.

This aggregation is structural rather than inferential. The script is not fitting a model to reconcile duplicates; it is simply enforcing a unique row per subject, clonotype, and time point.

### 7. Log transformation

After aggregation, the script computes

\[
\texttt{log\_freq} = \ln(\texttt{freq} + \varepsilon).
\]

This step is important because downstream dynamical analyses are typically carried out in log-abundance space. The addition of \(\varepsilon\) prevents undefined values when `freq = 0`:

\[
\ln(0 + \varepsilon) = \ln(\varepsilon).
\]

The script uses the natural logarithm. Therefore downstream displacements are naturally written as

\[
\Delta x = \ln f(t_1) - \ln f(t_0),
\]

which, when \(\varepsilon\) is negligible relative to the frequency scale, approximates the log fold-change

\[
\Delta x \approx \ln\left(\frac{f(t_1)}{f(t_0)}\right).
\]

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

If \(N_{\mathrm{clone}}\) denotes the number of unique subject–clone pairs, then

\[
N_{\mathrm{clone}} = \#\{(s,c): \exists t \text{ such that } (s,c,t) \text{ is present in } \texttt{trajectories\_long.csv}\}.
\]

### `trajectory_build_report.md`

This report contains run metadata and basic output sizes, including:

- timestamp
- input file path
- selected frequency column
- epsilon value
- number of rows after frequency validity filtering
- number of rows after aggregation
- number of subjects
- number of unique subject–clone pairs
- number of rows in the final trajectory table

This file is useful for provenance and auditability.

## Interpretation of counts

Let

- \(N_{\mathrm{in}}\) = number of rows in the raw input
- \(N_{\mathrm{valid}}\) = number of rows after removing invalid frequencies
- \(N_{\mathrm{agg}}\) = number of unique `(subject, aaSeqCDR3, time)` combinations after aggregation
- \(N_{\mathrm{clone}}\) = number of unique `(subject, aaSeqCDR3)` pairs

Then in general,

\[
N_{\mathrm{agg}} \le N_{\mathrm{valid}} \le N_{\mathrm{in}},
\]

with strict inequality whenever duplicate rows or invalid frequencies are present.

The final `trajectories_long.csv` table contains exactly \(N_{\mathrm{agg}}\) rows.

## Why log-frequency trajectories are useful downstream

If \(f(t)\) denotes clonotype frequency at time \(t\), defining

\[
x(t) = \ln(f(t) + \varepsilon)
\]

allows temporal changes to be written as

\[
\Delta x = x(t_1) - x(t_0).
\]

When \(\varepsilon\) is small relative to \(f\), this is approximately the log fold-change. This is advantageous because repertoire frequencies usually span orders of magnitude and many downstream dynamical models are naturally formulated in log-abundance space.

## Limitations

This script is intentionally minimal.

First, duplicate rows are collapsed using the arithmetic mean of frequency values. This is convenient and deterministic, but it is not the only possible summary.

Second, `observable` is propagated by the maximum. This preserves any evidence of observability but does not quantify uncertainty.

Third, zero frequencies are permitted and mapped to \(\ln(\varepsilon)\). This is numerically convenient, but the biological interpretation depends on upstream preprocessing.

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
6. collapses duplicate subject–clone–time entries by mean frequency and max observability
7. computes `log_freq = ln(freq + epsilon)`
8. exports a long-format trajectory table, an included-clone list, and a provenance report

It is therefore a foundational structural step for transition building and downstream dynamical analysis, while deliberately avoiding any inferential decision beyond basic mathematical consistency.

## Source

This documentation was written from the uploaded script `3-trajectory_builder.py`. 
