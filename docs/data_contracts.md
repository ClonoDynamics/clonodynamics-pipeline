# ClonoDynamics data contracts

## Purpose

This document defines the canonical data contracts used by the current ClonoDynamics pipeline.

It complements the scientific methods in [`methods/`](methods/README.md) by specifying:

- what the pipeline expects as input;
- how technical replicates, biological time points, clonotypes, trajectories, and transitions are represented;
- which fields define estimator support;
- which fields are operational annotations rather than biological states;
- which quantities may be missing by construction;
- which contracts are inherited by downstream ClonoDynamics modules.

The contracts below describe the current production workflow. They are intended to be stable interfaces between pipeline blocks.

---

# 1. Global conventions

## 1.1 Biological and technical units

ClonoDynamics distinguishes the following units:

| Unit | Meaning |
|---|---|
| biological subject | independent biological individual |
| biological time point | sampling time within a biological subject |
| technical replicate | independently measured repertoire from the same biological subject-time state |
| clonotype | amino-acid CDR3 sequence identified by `aaSeqCDR3` |
| transition | one subject × clonotype × endpoint-pair record |
| pseudo configuration | randomized reassembly of technical measurements; **not** a biological subject |

Technical replicates must never be counted as independent biological subjects.

Pseudo configurations must never be counted as biological replicates.

## 1.2 Logarithms

All log-frequency variables used by the current pipeline use the natural logarithm.

## 1.3 Missing observed log-frequency

Observed log-frequency is defined only for a positive observed count.

If a technical-replicate count is zero:

```text
observed relative frequency = 0
observed log-frequency      = undefined / missing
```

The pipeline does not introduce a finite pseudocount-based log-frequency for a zero-count replicate.

This missingness is expected and is part of the measurement-support contract.

## 1.4 Operational observation labels

Operational observation labels are defined from the endpoint observation statistic and an operational threshold `alpha`.

At threshold `alpha`:

```text
T = p_value < alpha
F = p_value >= alpha
```

For a transition:

```text
TT = T at t0 and T at t1
TF = T at t0 and F at t1
FT = F at t0 and T at t1
FF = F at t0 and F at t1
```

These are **operational observation states**.

They must not be interpreted as:

- biological presence or absence;
- clonotype birth or extinction;
- biological persistence or disappearance.

The primary reference threshold is:

```text
alpha = 0.05
```

Alternative thresholds are evaluated downstream without refitting the latent-state model.

---

# 2. Raw repertoire-table contract

ClonoDynamics starts from processed clonotype-level repertoire tables.

Raw FASTQ/BAM processing is outside the current ClonoDynamics analysis contract.

## 2.1 Required columns

Each repertoire table must contain at least:

| Column | Required | Meaning |
|---|---:|---|
| `aaSeqCDR3` | yes | amino-acid CDR3 clonotype identifier |
| `readCount` | yes | observed sequencing read count for the clonotype |
| `readFraction` | no | optional input field; the pipeline can reconstruct relative frequency from `readCount` |

Within a repertoire, repeated rows with the same `aaSeqCDR3` are collapsed by summing `readCount`.

`readFraction`, when present, is not treated as authoritative if it disagrees with counts.

## 2.2 Count requirements

`readCount` must be numeric and non-negative.

For the observed repertoire characterization, non-positive rows are excluded from the positive-count repertoire representation.

For paired state inference, zero counts arise naturally after outer-joining the two technical replicates on the union of detected clonotypes.

## 2.3 Sequencing depth

For one repertoire:

```text
sequencing depth = sum(readCount)
```

Observed relative frequency is reconstructed as:

```text
observed frequency = readCount / sequencing depth
```

---

# 3. Genuine longitudinal filename contract

The canonical longitudinal filename pattern is:

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

Interpretation:

```text
1_2-1.tsv
| | `-- technical replicate 1
| `---- biological time point 2
`------ biological subject 1
```

The two technical replicates for one state are therefore identified by the same subject and time and by replicate labels `1` and `2`.

The current production readers support the table formats used by the pipeline, including:

```text
.tsv
.txt
.csv
.parquet
.pq
.feather
.arrow
```

Compressed text formats may be accepted where supported by the corresponding reader.

## 3.1 Complete replicate pairs

Step 2 performs latent-state inference only on complete technical-replicate pairs.

A complete pair has:

```text
same biological subject
same biological time
replicate 1 present
replicate 2 present
```

Incomplete states may be documented by upstream characterization but cannot produce a paired Step-2 latent-state fit.

---

# 4. Pseudo-reference source-data contract

The current pseudo-longitudinal reference uses one biological source represented by 12 independent technical repertoire measurements.

The canonical current naming convention is:

```text
1_1-1.tsv
1_1-2.tsv
...
1_1-12.tsv
```

These 12 measurements are **not** 12 biological subjects and are **not** 12 biological time points.

They originate from the same biological source and are reassembled computationally into:

```text
6 pseudo-time positions × 2 technical replicates
```

The current production ensemble contains:

```text
2,000 pseudo configurations
```

Twelve technical measurements define:

```text
66 unique unordered measurement pairs
```

The reusable pair bank fits each unique pair once and reuses the result across pseudo configurations.

Pseudo time and pseudo lag are technical randomization coordinates and must not be interpreted as biological elapsed time.

---

# 5. Step-2 clonotype-state contract

Step 2 constructs one paired-replicate state table for each successfully fitted biological subject-time pair.

## 5.1 Clonotype universe

The two replicate tables are outer-joined by `aaSeqCDR3`.

Therefore the state table contains the union of clonotypes detected in either replicate.

A clonotype detected in only one replicate receives:

```text
positive count in one replicate
zero count in the other replicate
```

A synthetic `(0,0)` clonotype row is never created.

## 5.2 Canonical state identity

For the genuine longitudinal workflow, the state identity is:

```text
subject × time × aaSeqCDR3
```

The production state table must contain at most one row per state identity.

Duplicate state keys are an upstream contract violation.

## 5.3 Representation layers

The Step-2/Step-4 state contract preserves conceptually distinct layers.

### Model-based latent state

Canonical fields include a paired-replicate latent abundance representation, including the posterior median log-frequency used downstream as the principal latent conditioning coordinate.

Conceptually:

```text
x_latent = model-based paired-replicate latent log-frequency
```

The latent state is not a directly observed abundance.

### Replicate-resolved observed state

Observed quantities remain separate by technical replicate.

Conceptually:

```text
count_rep1
count_rep2
frequency_rep1
frequency_rep2
x_observed_rep1
x_observed_rep2
```

Observed log-frequency is missing when the corresponding count is zero.

### Observation-model quantities

The state layer also retains quantities such as:

```text
p_value
p_detect_state
p_dropout_state
operational observability at the reference alpha
```

`p_detect_state` and operational observability are not interchangeable.

## 5.4 Full posterior cache

Where requested, Step 2 stores full discrete posterior information under:

```text
per_clone_latent_logP/
```

These caches may be used by Step 5 for endpoint posterior uncertainty propagation.

They are not required for the primary observed AB/BA forward estimator or the primary cross-replicate fluctuation estimator.

---

# 6. Step-4 trajectory contract

Step 4 assembles the longitudinal state table without refitting the latent model.

The canonical output is:

```text
multirepresentation_trajectories_long.parquet
```

The genuine longitudinal state key is:

```text
subject × aaSeqCDR3 × time
```

The trajectory table preserves:

- latent state fields;
- replicate-resolved observed fields;
- observation-model fields;
- operational observation annotations;
- posterior uncertainty fields when available.

Step 4 is an assembly layer.

It does not:

- estimate longitudinal drift;
- estimate fluctuation covariance;
- refit latent states;
- filter by TT/TF/FT/FF;
- remove states using posterior-width or detectability thresholds.

---

# 7. Step-5 transition contract

The canonical transition output is:

```text
longitudinal_transitions.parquet
```

## 7.1 Transition identity

A transition represents:

```text
one biological subject
× one clonotype
× one ordered endpoint pair (t0, t1)
```

with:

```text
t0 < t1
dt = t1 - t0
```

The production genuine-longitudinal workflow uses:

```text
dt = 1, 2, 3, 4, 5
```

A transition is constructed only when the clonotype has Step-4 states at both endpoints.

Step 5 does not synthesize a both-zero state that was absent from the state table.

## 7.2 Latent transition fields

Conceptually important latent fields include:

```text
x0_latent
x1_latent
xmid_latent
xstar_latent
dx_latent
dt
```

`xmid_latent` is the primary conditioning coordinate for Step 11.

`xstar_latent` is retained as a model-based conditioning/sensitivity representation.

`dx_latent` is a representation-level displacement and is not the primary observed forward estimator.

## 7.3 Replicate-resolved transition fields

Observed replicate-resolved fields include endpoint abundances and displacements such as:

```text
x_observed_rep1_t0
x_observed_rep1_t1
x_observed_rep2_t0
x_observed_rep2_t1

dx_observed_rep1
dx_observed_rep2
```

Conceptually:

```text
dx_observed_rep1 = x_observed_rep1_t1 - x_observed_rep1_t0
dx_observed_rep2 = x_observed_rep2_t1 - x_observed_rep2_t0
```

A replicate-specific displacement is defined only when that replicate is positive at both endpoints.

## 7.4 Endpoint operational fields

The transition table preserves the continuous endpoint observation statistics:

```text
p_value_t0
p_value_t1
```

and the reference operational class derived from the reference `alpha`.

The continuous endpoint values are retained so that Step 16 can reconstruct operational classes at alternative alpha values without rebuilding the transition universe.

## 7.5 No primary class filtering during assembly

Step 5 constructs the generic transition universe.

It does not restrict the table to:

```text
TT
T0PLUS
NON_FF
common4
forward AB support
forward BA support
```

Those are downstream, estimand-specific support definitions.

---

# 8. Estimator-support contract

Support flags define whether the observed measurements required by a specific estimator are available.

They are not quality scores.

## 8.1 Forward AB eligibility

The AB fold uses:

```text
conditioning = observed replicate 1 at t0
outcome      = observed displacement in replicate 2
```

Therefore `forward_ab_eligible` requires:

```text
replicate 1 positive at t0
replicate 2 positive at t0
replicate 2 positive at t1
```

## 8.2 Forward BA eligibility

The BA fold uses:

```text
conditioning = observed replicate 2 at t0
outcome      = observed displacement in replicate 1
```

Therefore `forward_ba_eligible` requires:

```text
replicate 2 positive at t0
replicate 1 positive at t0
replicate 1 positive at t1
```

## 8.3 `common4`

`common4` requires both technical replicates to be positive at both transition endpoints:

```text
replicate 1 positive at t0
replicate 1 positive at t1
replicate 2 positive at t0
replicate 2 positive at t1
```

On `common4`, both observed replicate displacements are defined.

This is the primary support for Step 11 cross-replicate fluctuation dynamics.

## 8.4 Support versus operational class

The following concepts must remain separate:

```text
forward_ab_eligible / forward_ba_eligible / common4
    = measurement availability for an estimator

TT / TF / FT / FF
    = operational observation classes defined from p_value and alpha
```

The primary Step-9 and Step-11 analyses are not defined by a TT-only filter.

---

# 9. Step-9 forward-analysis contract

The primary one-week forward estimator uses reciprocal observed folds.

```text
AB:
x_observed_rep1_t0 -> dx_observed_rep2

BA:
x_observed_rep2_t0 -> dx_observed_rep1
```

AB and BA are estimated separately.

They are combined only after abundance binning with exact equal fold weight:

```text
combined = 0.5 × AB + 0.5 × BA
```

AB and BA row counts are not used as estimator weights.

The primary Step-9 forward domain is not restricted by TT/TF/FT/FF class.

---

# 10. Step-11 fluctuation-analysis contract

The primary Step-11 fluctuation quantity is the signed cross-replicate covariance of the two observed replicate displacements, conditioned on:

```text
xmid_latent
dt
common4
```

The corresponding code-level field is:

```text
cross_cov
```

Two complementary quantities are:

```text
same_var_mean
replicate_specific_excess
```

with the conceptual relation:

```text
replicate_specific_excess = same_var_mean - cross_cov
```

`cross_cov` remains signed and is never clipped at zero.

The primary Step-11 analysis does not impose an operational TT-only filter.

---

# 11. Step-13 temporal-analysis contract

Step 13 consumes Step-11 subject-level lag × abundance information.

The primary temporal analysis freezes:

```text
one abundance core
one complete-case biological-subject cohort
```

before fitting temporal summaries across:

```text
dt = 1, 2, 3, 4, 5
```

The primary weighting is:

```text
equal abundance-bin weight within subject
then equal biological-subject weight within lag
```

Step 12 is not subtracted from Step 11 before Step 13.

The pseudo technical reference is a comparison benchmark, not a preprocessing correction.

---

# 12. Positive-control data contract

The positive-control branch generates semi-synthetic count repertoires conditioned on the finalized empirical design.

The generator preserves, according to the frozen calibration:

- subject structure;
- eligible visits;
- replicate-specific sequencing depths;
- empirical-positive clonotype support;
- relevant empirical analysis support.

Generated scenario repertoires use the same basic input schema as ordinary repertoire tables:

```text
aaSeqCDR3
readCount
readFraction
```

Known synthetic trajectories are stored separately under a ground-truth output.

Ground truth is **never** supplied to Step 2 or to the downstream production inference chain.

The production recovery path is:

```text
Step 2 -> Step 4 -> Step 5 -> Step 11 -> Step 13
```

---

# 13. Missingness contract

Not every field is expected to be finite on every transition.

Examples of legitimate missingness include:

- observed log-frequency when the corresponding replicate count is zero;
- replicate displacement when that replicate is not positive at both endpoints;
- forward AB quantities outside `forward_ab_eligible`;
- forward BA quantities outside `forward_ba_eligible`;
- paired displacement covariance outside `common4`.

Completeness must therefore be evaluated **conditional on estimator support**.

A field being missing outside its required support is not, by itself, a data-quality failure.

---

# 14. Interpretation boundaries

The following interpretations are outside the data contract:

| Field / construct | Must not be interpreted as |
|---|---|
| `x_latent` | directly observed abundance or exact biological ground truth |
| zero observed count | biological clonotype absence |
| T/F | biological presence/absence |
| TF/FT | biological extinction/emergence |
| `common4` | quality class |
| pseudo configuration | biological subject |
| pseudo lag | elapsed biological time |
| pseudo randomization percentile | biological confidence interval |
| `cross_cov` | variance with negative values truncated away |
| positive-control ground truth | an input to production state inference |

---

# 15. Downstream stability

Files and fields identified in this document as canonical contracts should be treated as stable interfaces for downstream projects.

In particular, a downstream project such as **ClonoDynamics Reference Space** should consume finalized ClonoDynamics outputs rather than reconstructing internal state definitions independently.

Changes to canonical field semantics, support definitions, or primary output filenames should be treated as interface changes and documented in the repository release notes.
