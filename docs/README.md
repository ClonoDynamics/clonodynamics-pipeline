# ClonoDynamics documentation

This directory contains the scientific, computational, and reproducibility documentation for **ClonoDynamics**.

The root repository [`README.md`](../README.md) provides the project overview, repository architecture, and main execution entry point. The material under `docs/` provides the detailed documentation needed to understand the scientific contracts, canonical outputs, reproducibility rules, and extended methods of the framework.

---

## Documentation map

### 1. Scientific methods

Detailed scientific and computational methods are collected under:

[`methods/`](methods/README.md)

The methods documentation is divided into three complementary parts:

- **[Primary analysis: Steps 1–13](methods/01_primary_analysis.md)**
  Repertoire characterization, paired-replicate state inference, multirepresentation trajectory and transition construction, pseudo-longitudinal technical calibration, replicate-decoupled forward dynamics, cross-replicate fluctuation dynamics, longitudinal-versus-pseudo comparisons, and temporal scaling.

- **[Synthetic positive controls](methods/02_positive_controls.md)**
  Support-conditioned semi-synthetic validation, frozen empirical calibration, oracle dose construction, synthetic repertoire generation, recovery through the production pipeline, verification, and recovery summaries.

- **[Robustness and sensitivity controls: Steps 14–16](methods/03_robustness_controls.md)**
  Calendar-position and interval-composition controls, fixed-threshold operational observation-domain sensitivity, and numerical observation-threshold robustness.

---

### 2. Data contracts

[`data_contracts.md`](data_contracts.md)

Defines the canonical input and intermediate-data contracts used by ClonoDynamics, including:

- repertoire-table requirements;
- longitudinal and pseudo-reference filename conventions;
- biological versus technical units;
- latent versus observed representations;
- Step-5 transition structure;
- AB/BA estimator support;
- `common4`;
- operational T/F and TT/TF/FT/FF classes;
- legitimate support-dependent missingness;
- interpretation boundaries.

This document should be consulted when integrating ClonoDynamics with external or downstream software.

---

### 3. Output reference

[`outputs_reference.md`](outputs_reference.md)

Provides a practical reference to the principal canonical outputs of the pipeline.

For each major analysis block it describes:

- canonical result directories;
- principal output files;
- row/object granularity;
- scientific role;
- main downstream consumers.

It also identifies the outputs that should be preferred as stable interfaces for downstream projects such as **ClonoDynamics Reference Space**.

---

### 4. Reproducibility

[`reproducibility.md`](reproducibility.md)

Describes the reproducibility policy of the project, including:

- the canonical master entry point;
- living repository versus frozen manuscript release;
- input-data immutability;
- environment capture;
- provenance records;
- source-code identity;
- `resume`, `fresh`, and pseudo `repair`;
- bootstrap, pseudo-randomization, and positive-control seeds;
- recommended reconstruction order;
- preservation of result trees;
- manuscript-release checklist.

---

### 5. Troubleshooting

[`troubleshooting.md`](troubleshooting.md)

Reserved for verified issues encountered during real ClonoDynamics use.

This document is intentionally minimal at present. Problems should be added only after their cause and resolution have been established.

---

## Recommended reading paths

### New user

```text
root README.md
    -> docs/data_contracts.md
    -> docs/methods/README.md
    -> docs/outputs_reference.md
```

### Reproducing the published study

```text
root README.md
    -> docs/reproducibility.md
    -> docs/methods/
    -> docs/outputs_reference.md
```

### Developing a downstream analysis

```text
docs/data_contracts.md
    -> docs/outputs_reference.md
    -> relevant scientific method document
```

---

## Documentation versioning

The files in `docs/` are the **living documentation** of the current ClonoDynamics software.

They may evolve when:

- repository paths are reorganized;
- orchestration is improved;
- implementation details are clarified;
- new public outputs are introduced;
- explanatory text is improved.

The exact documentation corresponding to a submitted or published manuscript should be frozen together with the matching source code in a tagged release, for example:

```text
v1.0.0-paper
```

The release should preserve, at minimum:

```text
README.md
code/
docs/
environment specification
LICENSE
CITATION.cff
```

This keeps the active repository maintainable while preserving an immutable computational and methodological record for the reported study.

---

## Project relationship

ClonoDynamics provides the validated upstream framework:

```text
processed repertoire data
    -> replicate-resolved state inference
    -> longitudinal transitions
    -> validated longitudinal dynamics
```

A separate downstream project, **ClonoDynamics Reference Space**, is intended to consume finalized ClonoDynamics outputs to construct:

```text
quantitative dynamic encoding
    -> healthy physiological reference space
    -> projection and analysis of perturbed repertoire dynamics
```

The downstream project should record the exact ClonoDynamics release from which its inputs were generated.
