# ClonoDynamics plotting package

This directory is the canonical plotting layer for the reorganized ClonoDynamics repository.
Plotters are read-only with respect to scientific analysis outputs: they consume finalized
analysis tables and render publication figures/source-data exports without rerunning the
underlying inferential pipeline.

The directory mirrors the four analysis blocks:

```
05_plotting/
├── 01_core/
│   ├── plot_step01_repertoire_characterization.py
│   ├── plot_step03_latent_state_diagnostics.py
│   └── plot_step06_transition_support.py
├── 02_pseudo_reference/
│   └── plot_steps07_08_pseudo_reference.py
├── 03_validation/
│   ├── plot_steps09_10_forward_dynamics.py
│   ├── plot_steps11_12_fluctuation_dynamics.py
│   ├── plot_step13_temporal_scaling.py
│   └── positive_controls/
│       └── plot_positive_controls.py
└── 04_controls/
    ├── plot_steps14_16_robustness_summary.py
    └── details/
        ├── plot_step14_interval_position_details.py
        ├── plot_step15_observation_domain_details.py
        └── plot_step16_threshold_robustness_details.py
```

## Canonical default plotting workflow

The default figure orchestrator runs the following plotting stages:

1. Step 1 repertoire characterization.
2. Step 3 latent-state / observation-model diagnostics.
3. Step 6 transition-support characterization.
4. Pseudo-reference Steps 7–8.
5. Validation Steps 9–10 (forward dynamics).
6. Validation Steps 11–12 (fluctuation dynamics).
7. Validation Step 13 (temporal scaling).
8. Synthetic positive-control summary, when available.
9. Integrated Steps 14–16 robustness summary.

The integrated controls plotter is the canonical default for the current final controls
figure architecture. It produces the integrated main Figure 8 and assembled/separate
Supplementary Figures 9–11. The three single-step control plotters are retained in
`04_controls/details/` because they expose additional diagnostic panels, but they overlap
substantially with the integrated output and are not run by default.

## Preservation policy

The supplied plotting logic, figure dimensions, axis limits, fonts, line widths, markers,
colours, annotations, legends, support thresholds, export DPI and default output formats
have been preserved. Repository cleanup changes only filenames, code location and example
paths; no scientific or graphical setting is intentionally altered.

## Figure outputs

The orchestrator writes into a user-selected figure root using the same block structure:

```
figures/
├── 01_core/
├── 02_pseudo_reference/
├── 03_validation/
└── 04_controls/
```

Each plotter retains its own internal output organization (for example `main_figure/`,
`supplementary_figure/`, `figure_source_data/`, or Step-3 `main/` and `supplementary/`).

## Not generated here

A study-design / methodological schematic, if used as manuscript Figure 1, is not generated
by these scripts. No standalone plotting stage is required for Steps 2, 4 or 5: their relevant
state/transition information is represented downstream by the Step-3 and Step-6/dynamics
plotters.
