# ClonoDynamics reproducibility guide

## Purpose

This document defines the reproducibility policy for the current ClonoDynamics repository.

The aim is to make three levels of reproducibility explicit:

1. **workflow reproducibility** — rerunning the same analysis design from the same inputs;
2. **study reproducibility** — reconstructing the analysis reported in the associated manuscript;
3. **software reproducibility** — preserving the exact code, environment, documentation, and provenance corresponding to a public release.

Scientific definitions are documented in [`methods/`](methods/README.md). Input/output interfaces are documented in [`data_contracts.md`](data_contracts.md) and [`outputs_reference.md`](outputs_reference.md).

---

# 1. Reproducibility unit

A ClonoDynamics result is defined by the combination of:

```text
input data
+ scientific code
+ orchestration logic
+ analysis parameters
+ random seeds
+ software environment
+ recorded provenance
+ matching methods documentation
```

A filename alone is not a complete reproducibility record.

---

# 2. Canonical entry point

The recommended entry point is:

```bash
python3 code/orchestrator_clonodynamics.py
```

The master contains the validated orchestration engines for:

```text
core
pseudo
validation
positive_controls
controls
figures
```

The master does not replace the scientific scripts under:

```text
code/01_core/
code/02_pseudo_reference/
code/03_validation/
code/04_controls/
code/05_plotting/
```

The individual scientific scripts remain directly executable for development and methodological inspection.

For routine production runs, use the master so that dependency checks, output policies, provenance, and logs remain consistent.

---

# 3. Living repository versus frozen study release

The default branch should be treated as **living software documentation and implementation**.

Repository paths, documentation, plotting layout, or implementation details may evolve after publication.

The exact software-and-methods state corresponding to a manuscript should therefore be frozen in a tagged release.

Recommended convention:

```text
v1.0.0-paper
```

The tagged release should freeze together:

```text
code/
docs/
README.md
environment specification
LICENSE
CITATION.cff
```

The public manuscript should cite the corresponding release and, when available, an archived DOI.

---

# 4. Input-data immutability

Reproducing a reported analysis requires the same processed clonotype-level input tables or a versioned dataset that is demonstrably equivalent.

For the genuine longitudinal analysis, the canonical processed input is the collection of subject-time-replicate repertoire tables.

For the pseudo reference, the canonical input is the 12-measurement technical dataset.

Input datasets should be archived independently of the code repository when their size or access conditions make direct GitHub storage inappropriate.

For a frozen public study release, retain or publish:

- dataset accession or DOI;
- input file manifest;
- file sizes;
- cryptographic hashes where practical;
- any documented exclusions or unavailable samples.

Do not silently replace a failed or unavailable repertoire with a duplicated or synthetic placeholder.

---

# 5. Environment capture

The public repository should provide one pinned environment specification, for example:

```text
requirements.txt
```

or:

```text
environment.yml
```

At minimum, the environment should constrain the packages used by the production analysis, including:

```text
numpy
pandas
scipy
polars
pyarrow
matplotlib
```

A frozen manuscript release should record:

- Python version;
- package versions;
- operating-system/platform information when relevant;
- exact environment file committed in the release.

The orchestrators additionally record the Python executable and version where implemented.

Exact numerical identity across substantially different scientific-package versions should not be assumed.

---

# 6. Provenance recorded by the workflow

Depending on the block, ClonoDynamics records objects such as:

```text
00_orchestrator_config.json
00_orchestrator_step_manifest.csv
00_run_config.json
00_run_signature.json
00_input_manifest.csv
00_analysis_summary.csv
00_pipeline_manifest.csv
logs/
```

These records should be retained with the analysis outputs.

They may encode:

- input and output paths;
- selected workflow/block;
- operational observation threshold;
- pseudo configuration count;
- bootstrap parameters;
- random seeds;
- script hashes;
- analysis signatures;
- command lines;
- timestamps;
- output-completeness state.

For a reported analysis, these files are part of the scientific provenance record and should not be discarded simply because the principal result tables have already been generated.

---

# 7. Code identity

The current orchestration layer records or checks source hashes where required.

For a frozen run, source identity should be established through both:

```text
Git commit / release tag
```

and:

```text
recorded script hashes
```

Version strings are useful human-readable identifiers but are not a substitute for source identity.

If code hashes differ from a prior run, an orchestrator may reject `resume` even when filenames are unchanged.

This is intentional.

---

# 8. Existing-output policies

ClonoDynamics uses conservative policies for existing outputs.

## 8.1 Resume

Use `resume` when the intention is to continue a compatible analysis without recomputing outputs already verified as complete.

Typical behavior:

```text
complete compatible step -> skip
missing/incomplete step   -> run
compatible checkpoint     -> reuse when supported
```

`resume` is the normal choice after an interruption when the dataset, scientific settings, and code identity are unchanged.

A safe resume may be rejected when provenance does not match.

Do not bypass such a rejection by manually editing config files.

## 8.2 Fresh

Use `fresh` when the intention is to regenerate the selected analysis block from the beginning.

A fresh run removes only outputs owned by that block, according to the corresponding orchestrator policy.

It does not delete the raw input dataset.

Use `fresh` after deliberate changes such as:

- corrected input data;
- scientific-code changes;
- altered primary parameters;
- intentionally different analysis settings.

## 8.3 Repair

The pseudo-reference upstream has an additional repair/resume mode because rebuilding the 66 pair fits and 2,000 compact pseudo configurations can be expensive.

Repair is appropriate when:

- reusable upstream results are scientifically compatible;
- some metadata, manifests, configurations, or downstream aggregate products are incomplete;
- valid expensive intermediate results should be preserved.

Repair must not be used to conceal a scientific contract mismatch.

If the stored pseudo upstream is incompatible with the requested analysis contract, rebuild it deliberately.

---

# 9. Randomness and inferential units

ClonoDynamics uses several forms of resampling/randomization that have different meanings.

They must not be conflated.

## 9.1 Biological subject bootstrap

The genuine longitudinal uncertainty analyses use biological subjects as the resampling unit.

Where the production configuration uses:

```text
2,000 bootstrap resamples
seed = 123
```

the same subject draw is propagated jointly across the relevant lags, domains, or estimators when paired comparison is required.

Always verify the recorded `00_run_config.json` or equivalent output for the specific run.

## 9.2 Pseudo randomization

The pseudo ensemble contains:

```text
2,000 pseudo configurations
```

A pseudo configuration is a randomized technical reassembly of the same biological source.

Across-configuration percentiles are therefore:

```text
technical randomization intervals
```

not biological confidence intervals.

No biological subject bootstrap is defined by treating configuration IDs as subjects.

## 9.3 Positive-control simulation

The production support-conditioned positive-control dose series uses a fixed simulation design.

The reported production configuration includes:

```text
master simulation seed = 20260918
fast component sigma   = 0.20
temporal bootstrap seed = 123
```

Scenario-specific dose parameters are read from the generated oracle scenario table rather than re-entered manually.

The same simulation/random-field contract should be preserved when reproducing the reported validation.

## 9.4 Step-5 posterior propagation

When full-posterior endpoint propagation is enabled, Step 5 uses a deterministic posterior-sampling seed and records it in the run configuration.

The production core orchestrator currently exposes this seed explicitly.

A reproduction should use the value recorded in the original run rather than assuming a default.

---

# 10. Reproducing the main empirical workflow

A typical reported-study reconstruction proceeds in the following order.

## 10.1 Genuine longitudinal core

Launch:

```bash
python3 code/orchestrator_clonodynamics.py
```

Select:

```text
Analysis
core
```

Provide:

```text
genuine longitudinal repertoire dataset
destination for longitudinal results
reference operational alpha
```

The current primary reference alpha is:

```text
0.05
```

This produces Steps 1–6.

## 10.2 Pseudo technical reference

From the same master, select:

```text
Analysis
pseudo
```

Provide the pseudo source dataset and pseudo result root.

This constructs:

- the pseudo design;
- the reusable 66-pair bank;
- the 2,000 compact configurations;
- Step 7;
- the Step-8 transition bank;
- Step 8.

## 10.3 Validation

Select:

```text
Analysis
validation
```

Provide:

```text
completed longitudinal core results
completed pseudo-reference results
destination for validation results
```

This runs the current Steps 9–13 dependency graph.

## 10.4 Positive controls

Select:

```text
Analysis
positive_controls
```

Provide:

```text
corrected longitudinal repertoire dataset
completed core results
completed validation results
positive-control destination
```

This freezes the finalized empirical calibration and runs the support-conditioned validation without rebuilding the finalized empirical analysis.

## 10.5 Controls

Select:

```text
Analysis
controls
```

Provide the finalized longitudinal/validation result root and selected output destination.

This runs Steps 14–16.

## 10.6 Figures

Return to the master and select:

```text
Figures
```

Choose one or more figure blocks.

Figure generation reads finalized analysis outputs and must not alter the scientific analysis tables.

---

# 11. Result-tree preservation

For a manuscript-associated release, retain the result directories used to generate the publication.

Do not retain only final PDFs.

At minimum, preserve:

- canonical analysis tables;
- run configs/signatures;
- manifests;
- logs needed to explain failures or exclusions;
- figure source-data tables;
- final figures;
- positive-control oracle/scenario metadata;
- pseudo randomization metadata.

Large intermediate caches that are reproducible and not needed for audit may be archived separately, but this should be documented rather than silently omitted.

---

# 12. Reproducing figures

Figures are generated from finalized result tables.

The figure workflow should be considered a read-only rendering layer.

Recommended sequence:

```text
finalize analysis
freeze result tree
run figure orchestrator
archive figure source data
assemble manuscript figure if required
```

If a plotter changes but the analysis does not, the figure-only change should be documented separately.

A figure rerender must not be presented as a rerun of the scientific analysis.

---

# 13. Reproducibility boundaries

A successful rerun does not automatically imply byte-for-byte identity with a historical run.

Differences can arise from:

- paths;
- timestamps;
- package versions;
- compression-library versions;
- platform-specific floating-point behavior;
- metadata serialization;
- plotting backends;
- intentionally revised documentation.

Scientific reproducibility should therefore be evaluated at the level appropriate to the output:

```text
same input design
same fitted/analysis contract
same support
same cohort
same abundance geometry
same estimator
same random-seed design
numerically compatible scientific results
```

Exact cryptographic identity is appropriate for frozen source files and archived inputs, but it may not be a meaningful requirement for every regenerated compressed data object.

---

# 14. Manuscript-release checklist

Before creating the paper-associated public release, verify:

- [ ] root `README.md` reflects the final public repository layout;
- [ ] `docs/methods/` matches the scientific analysis being submitted;
- [ ] `data_contracts.md` matches the production field/support definitions;
- [ ] `outputs_reference.md` matches the canonical public outputs;
- [ ] the master orchestrator passes its self-check;
- [ ] all five analysis blocks run under the intended public layout;
- [ ] the figure workflow runs from finalized outputs;
- [ ] a pinned environment file is present;
- [ ] `LICENSE` is present;
- [ ] `CITATION.cff` is present or prepared for the final citation;
- [ ] input data accessions/DOIs are recorded;
- [ ] relevant result/provenance objects are archived;
- [ ] the exact Git commit is tagged;
- [ ] the release tag is cited in the manuscript/repository metadata.

Recommended paper-associated tag:

```text
v1.0.0-paper
```

---

# 15. Relationship to downstream projects

Downstream projects should record the exact ClonoDynamics release used to produce their inputs.

For example:

```text
ClonoDynamics Reference Space
    consumes finalized outputs from ClonoDynamics v1.0.0-paper
```

The downstream project should not silently redefine the upstream ClonoDynamics state, support, or estimator contracts.

This separation preserves a clear provenance chain:

```text
processed repertoire tables
    -> ClonoDynamics
    -> validated longitudinal-dynamics outputs
    -> downstream quantitative reference-space analysis
```
