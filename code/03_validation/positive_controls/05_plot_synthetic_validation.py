#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot synthetic temporal-validation summaries and export publication source data.

Purpose and workflow position
-----------------------------
Read the numerical summaries from 04_summarize_temporal_recovery.py. This program
does not generate synthetic data, fit temporal slopes or rerun bootstrap draws.
Use it as the single figure-export stage; suppress the summarizer's legacy
figure files with --no-png --no-pdf --no-svg.

Inputs in --summary-dir
-----------------------
00_validation_summary.csv
01_temporal_profiles.csv
02_incremental_recovery_summary.csv
03_recovery_fit_summary.csv
The script checks expected columns and scenario membership. It assumes the
upstream paired-bootstrap audit has already established valid matching.

Current implemented display (preserved unchanged)
------------------------------------------------
A: recovered cross-replicate covariance profiles versus lag for all four doses.
B: BASELINE-ADJUSTED incremental recovery versus the oracle slope increment;
   intervals are paired scenario-minus-plateau subject-bootstrap intervals.
Optional supplementary panel: absolute oracle versus recovered slopes,
   including the plateau and a descriptive linear fit.
IMPORTANT: the proposed alternative panel B (four absolute recovered slopes
versus prescribed DOSE) is NOT implemented in this source. Renaming and
rewriting its docstring must not be mistaken for that later plotting revision.
A separate panel update and source-data mapping are still required to adopt it.

Outputs under --outdir
----------------------
Table_validation_dose_response_publication.csv
SourceData_Figure_validation_A_temporal_profiles.csv
SourceData_Figure_validation_B_incremental_recovery.csv
Figure_validation_A_temporal_profiles.[pdf/png/svg]
Figure_validation_B_incremental_recovery.[pdf/png/svg]
With --include-supplementary:
  SuppFigure_validation_oracle_vs_recovered.[pdf/png/svg]
PDF and PNG are on by default; SVG is off. Panels are saved separately; the
program does not assemble a composite figure. The optional supplementary plot
has no dedicated source CSV in the unchanged original implementation.

Interpretation
--------------
Do not interpret the plateau's raw recovered/oracle ratio or label the recovery
factor as general accuracy. Negative confidence limits are retained. Oracle
and recovered quantities refer to different support-specific populations.
Subject-bootstrap intervals refer to the supplied simulated realization.

Execution
---------
python 05_plot_synthetic_validation.py --summary-dir SUMMARY_DIR \
    --outdir NEW_FIGURE_DIR --include-supplementary
Dependencies: Python >=3.9, NumPy, pandas and Matplotlib. Named files may be
overwritten; use a new output directory. Aesthetics remain editable below.

Provenance
----------
Documentation-only edition of 11-plot_support_conditioned_validation_v2.py.
All plotting functions, statistical inputs, defaults and output names remain
unchanged. The filename is now independent of version/figure numbering;
SCRIPT_MAP.csv retains the original identity and both SHA-256 hashes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


# =============================================================================
# USER-EDITABLE PUBLICATION AESTHETICS
# =============================================================================

FONT_FAMILY = "Arial"
FONT_SIZE = 8
AXIS_LABEL_SIZE = 8
TICK_LABEL_SIZE = 7
LEGEND_SIZE = 7
ANNOTATION_SIZE = 7

FIGURE_SIZE = (3.45, 2.85)  # inches
PNG_DPI = 600

LINE_WIDTH = 1.2
ERRORBAR_LINE_WIDTH = 1.0
MARKER_SIZE = 4.8
CAPSIZE = 2.0

# Distinct symbols; colors are left to the active Matplotlib color cycle.
SCENARIO_MARKERS = {
    "R0p00": "o",
    "R0p25": "s",
    "R0p50": "^",
    "R1p00": "D",
}

# Panel A
PANEL_A_LEGEND_NCOL = 2
PANEL_A_SHOW_ZERO_LINE = False
PANEL_A_YMIN = None
PANEL_A_YMAX = None

# Panel B
PANEL_B_SHOW_ZERO_LINE = True
PANEL_B_SHOW_IDENTITY = True
PANEL_B_SHOW_RECOVERY_FIT = True
PANEL_B_LABEL_POINTS = True

# Supplementary absolute recovery panel
SUPP_SHOW_R2_IN_LEGEND = False


# =============================================================================
# CONSTANTS
# =============================================================================

PRIMARY_ESTIMATOR = "equal_bin_equal_subject"
PRIMARY_METRIC = "cross_cov"

REQUIRED_FILES = {
    "summary": "00_validation_summary.csv",
    "profiles": "01_temporal_profiles.csv",
    "incremental": "02_incremental_recovery_summary.csv",
    "fits": "03_recovery_fit_summary.csv",
}


# =============================================================================
# HELPERS
# =============================================================================

def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_required(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return pd.read_csv(path)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": FONT_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def standardize_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    stem: str,
    *,
    pdf: bool,
    png: bool,
    svg: bool,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if pdf:
        fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    if png:
        fig.savefig(outdir / f"{stem}.png", dpi=PNG_DPI, bbox_inches="tight")
    if svg:
        fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def scenario_order(summary: pd.DataFrame) -> list[str]:
    d = summary.copy()
    d["target_R"] = pd.to_numeric(d["target_R"], errors="raise")
    return (
        d.sort_values(["target_R", "scenario_label"])["scenario_label"]
        .astype(str)
        .tolist()
    )


def marker_for(scenario: str, index: int) -> str:
    if scenario in SCENARIO_MARKERS:
        return SCENARIO_MARKERS[scenario]
    fallback = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]
    return fallback[index % len(fallback)]


def finite_numeric(x: Sequence[float]) -> np.ndarray:
    a = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    return a[np.isfinite(a)]


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def load_inputs(summary_dir: Path):
    data = {
        key: read_required(summary_dir / filename, filename)
        for key, filename in REQUIRED_FILES.items()
    }

    summary = data["summary"]
    profiles = data["profiles"]
    incremental = data["incremental"]
    fits = data["fits"]

    summary_required = {
        "scenario_label",
        "target_R",
        "q_bio",
        "realized_R",
        "oracle_slope_per_week",
        "recovered_slope_per_week",
        "recovered_bootstrap_q025",
        "recovered_bootstrap_q975",
        "delta_oracle_slope_vs_baseline",
        "delta_recovered_slope_vs_baseline",
        "incremental_recovery_fraction",
        "paired_bootstrap_q025_delta",
        "paired_bootstrap_q975_delta",
        "paired_bootstrap_positive_fraction",
        "incremental_recovery_fraction_boot_q025",
        "incremental_recovery_fraction_boot_q975",
        "n_complete_case_subjects",
        "n_core_bins",
        "core_x_min",
        "core_x_max",
    }
    require(
        summary_required.issubset(summary.columns),
        "00_validation_summary.csv missing columns: "
        + str(sorted(summary_required - set(summary.columns))),
    )

    profile_required = {
        "scenario_label",
        "target_R",
        "dt",
        "value",
        "bootstrap_q025",
        "bootstrap_q975",
        "n_subjects",
    }
    require(
        profile_required.issubset(profiles.columns),
        "01_temporal_profiles.csv missing columns: "
        + str(sorted(profile_required - set(profiles.columns))),
    )

    inc_required = {
        "scenario_label",
        "baseline_scenario",
        "target_R",
        "delta_oracle_slope_vs_baseline",
        "delta_recovered_slope_vs_baseline",
        "paired_bootstrap_q025_delta",
        "paired_bootstrap_q975_delta",
        "paired_bootstrap_positive_fraction",
        "incremental_recovery_fraction_point",
        "incremental_recovery_fraction_boot_q025",
        "incremental_recovery_fraction_boot_q975",
        "n_paired_bootstrap",
        "paired_bootstrap_status",
    }
    require(
        inc_required.issubset(incremental.columns),
        "02_incremental_recovery_summary.csv missing columns: "
        + str(sorted(inc_required - set(incremental.columns))),
    )

    fit_required = {"fit", "intercept", "slope", "r_squared"}
    require(
        fit_required.issubset(fits.columns),
        "03_recovery_fit_summary.csv missing columns: "
        + str(sorted(fit_required - set(fits.columns))),
    )

    return summary, profiles, incremental, fits


# =============================================================================
# PUBLICATION TABLE
# =============================================================================

def build_publication_table(
    summary: pd.DataFrame,
    incremental: pd.DataFrame,
) -> pd.DataFrame:
    d = summary.merge(
        incremental[
            [
                "scenario_label",
                "baseline_scenario",
                "paired_bootstrap_q025_delta",
                "paired_bootstrap_q975_delta",
                "paired_bootstrap_positive_fraction",
                "incremental_recovery_fraction_boot_q025",
                "incremental_recovery_fraction_boot_q975",
                "n_paired_bootstrap",
            ]
        ],
        on="scenario_label",
        how="left",
        suffixes=("", "_inc"),
        validate="one_to_one",
    )

    # Prefer incremental-table copies if merge created duplicates.
    for c in [
        "paired_bootstrap_q025_delta",
        "paired_bootstrap_q975_delta",
        "paired_bootstrap_positive_fraction",
        "incremental_recovery_fraction_boot_q025",
        "incremental_recovery_fraction_boot_q975",
        "n_paired_bootstrap",
    ]:
        inc_col = f"{c}_inc"
        if inc_col in d.columns:
            d[c] = d[inc_col]
            d = d.drop(columns=[inc_col])

    d["incremental_recovery_percent"] = 100.0 * pd.to_numeric(
        d["incremental_recovery_fraction"], errors="coerce"
    )
    d["incremental_recovery_ci025_percent"] = 100.0 * pd.to_numeric(
        d["incremental_recovery_fraction_boot_q025"], errors="coerce"
    )
    d["incremental_recovery_ci975_percent"] = 100.0 * pd.to_numeric(
        d["incremental_recovery_fraction_boot_q975"], errors="coerce"
    )

    columns = [
        "scenario_label",
        "target_R",
        "q_bio",
        "realized_R",
        "oracle_slope_per_week",
        "recovered_slope_per_week",
        "recovered_bootstrap_q025",
        "recovered_bootstrap_q975",
        "delta_oracle_slope_vs_baseline",
        "delta_recovered_slope_vs_baseline",
        "paired_bootstrap_q025_delta",
        "paired_bootstrap_q975_delta",
        "paired_bootstrap_positive_fraction",
        "incremental_recovery_fraction",
        "incremental_recovery_percent",
        "incremental_recovery_fraction_boot_q025",
        "incremental_recovery_fraction_boot_q975",
        "incremental_recovery_ci025_percent",
        "incremental_recovery_ci975_percent",
        "n_paired_bootstrap",
        "n_complete_case_subjects",
        "n_core_bins",
        "core_x_min",
        "core_x_max",
    ]

    return (
        d[columns]
        .sort_values(["target_R", "scenario_label"])
        .reset_index(drop=True)
    )


# =============================================================================
# FIGURE A — TEMPORAL DOSE RESPONSE
# =============================================================================

def plot_temporal_profiles(
    summary: pd.DataFrame,
    profiles: pd.DataFrame,
    outdir: Path,
    *,
    pdf: bool,
    png: bool,
    svg: bool,
) -> None:
    order = scenario_order(summary)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    for i, scenario in enumerate(order):
        d = profiles[profiles["scenario_label"].astype(str).eq(scenario)].copy()
        d["dt"] = pd.to_numeric(d["dt"], errors="raise").astype(int)
        d = d.sort_values("dt")

        require(
            d["dt"].tolist() == [1, 2, 3, 4, 5],
            f"{scenario}: expected dt=1..5 in temporal profile",
        )

        x = d["dt"].to_numpy(float)
        y = pd.to_numeric(d["value"], errors="raise").to_numpy(float)
        lo = pd.to_numeric(d["bootstrap_q025"], errors="raise").to_numpy(float)
        hi = pd.to_numeric(d["bootstrap_q975"], errors="raise").to_numpy(float)

        yerr = np.vstack([y - lo, hi - y])
        target_R = float(d["target_R"].iloc[0])

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            marker=marker_for(scenario, i),
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            elinewidth=ERRORBAR_LINE_WIDTH,
            capsize=CAPSIZE,
            label=f"{scenario} (R={target_R:.2f})",
        )

    if PANEL_A_SHOW_ZERO_LINE:
        ax.axhline(0.0, linestyle=":", linewidth=0.8)

    ax.set_xlabel(r"Temporal lag, $\Delta t$ (weeks)")
    ax.set_ylabel("Cross-replicate displacement covariance")
    ax.set_xticks([1, 2, 3, 4, 5])

    if PANEL_A_YMIN is not None or PANEL_A_YMAX is not None:
        ax.set_ylim(PANEL_A_YMIN, PANEL_A_YMAX)
    else:
        # Preserve a biologically interpretable zero baseline without drawing it.
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(0.0, lo), hi)

    ax.legend(
        frameon=False,
        ncol=PANEL_A_LEGEND_NCOL,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.01, 1.0, 0.15),
        mode="expand",
        borderaxespad=0.0,
        handlelength=2.0,
        columnspacing=1.0,
    )

    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "Figure_validation_A_temporal_profiles",
        pdf=pdf,
        png=png,
        svg=svg,
    )


# =============================================================================
# FIGURE B — BASELINE-ADJUSTED INCREMENTAL RECOVERY
# =============================================================================

def plot_incremental_recovery(
    summary: pd.DataFrame,
    incremental: pd.DataFrame,
    fits: pd.DataFrame,
    outdir: Path,
    *,
    pdf: bool,
    png: bool,
    svg: bool,
) -> None:
    baseline_rows = incremental[
        incremental["paired_bootstrap_status"].astype(str).eq("baseline")
    ]
    require(len(baseline_rows) == 1, "Expected exactly one baseline scenario")
    baseline = str(baseline_rows.iloc[0]["scenario_label"])

    d = incremental[
        ~incremental["scenario_label"].astype(str).eq(baseline)
    ].copy()
    d["target_R"] = pd.to_numeric(d["target_R"], errors="raise")
    d = d.sort_values(["target_R", "scenario_label"]).reset_index(drop=True)

    fit = fits[
        fits["fit"].astype(str).eq("baseline_adjusted_through_origin")
    ]
    require(
        len(fit) == 1,
        "03_recovery_fit_summary.csv must contain one baseline_adjusted_through_origin row",
    )
    attenuation = float(fit.iloc[0]["slope"])

    x = pd.to_numeric(
        d["delta_oracle_slope_vs_baseline"], errors="raise"
    ).to_numpy(float)
    y = pd.to_numeric(
        d["delta_recovered_slope_vs_baseline"], errors="raise"
    ).to_numpy(float)
    lo = pd.to_numeric(
        d["paired_bootstrap_q025_delta"], errors="coerce"
    ).to_numpy(float)
    hi = pd.to_numeric(
        d["paired_bootstrap_q975_delta"], errors="coerce"
    ).to_numpy(float)

    require(np.all(np.isfinite(x)), "Non-finite oracle increments")
    require(np.all(np.isfinite(y)), "Non-finite recovered increments")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    max_x = max(float(np.max(x)), 1e-8)
    max_y = max(float(np.max(y)), float(np.nanmax(hi)), 1e-8)
    grid_max = 1.08 * max(max_x, max_y)
    grid = np.linspace(0.0, grid_max, 200)

    if PANEL_B_SHOW_IDENTITY:
        ax.plot(
            grid,
            grid,
            linestyle="--",
            linewidth=0.9,
            label="Identity",
        )

    if PANEL_B_SHOW_RECOVERY_FIT:
        ax.plot(
            grid,
            attenuation * grid,
            linewidth=LINE_WIDTH,
            label=f"Baseline-adjusted recovery = {attenuation:.3f} × oracle",
        )

    for i, row in d.iterrows():
        yy = float(row["delta_recovered_slope_vs_baseline"])
        qlo = float(row["paired_bootstrap_q025_delta"])
        qhi = float(row["paired_bootstrap_q975_delta"])

        if np.isfinite(qlo) and np.isfinite(qhi):
            yerr = np.array([[yy - qlo], [qhi - yy]])
        else:
            yerr = None

        ax.errorbar(
            [float(row["delta_oracle_slope_vs_baseline"])],
            [yy],
            yerr=yerr,
            marker=marker_for(str(row["scenario_label"]), i + 1),
            markersize=MARKER_SIZE + 0.5,
            linewidth=0,
            elinewidth=ERRORBAR_LINE_WIDTH,
            capsize=CAPSIZE,
        )

        if PANEL_B_LABEL_POINTS:
            ax.annotate(
                str(row["scenario_label"]),
                (
                    float(row["delta_oracle_slope_vs_baseline"]),
                    yy,
                ),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=ANNOTATION_SIZE,
            )

    if PANEL_B_SHOW_ZERO_LINE:
        ax.axhline(0.0, linestyle=":", linewidth=0.8)

    ax.set_xlabel(f"Increment in oracle temporal slope vs {baseline}")
    ax.set_ylabel(f"Increment in recovered temporal slope vs {baseline}")

    # Preserve any negative lower CI (important for R0p25).
    finite_lows = finite_numeric(lo)
    ymin_data = float(np.min(finite_lows)) if finite_lows.size else 0.0
    ymin = min(0.0, ymin_data)
    ymax = max(grid_max * attenuation, max_y)
    ypad = 0.08 * max(ymax - ymin, 1e-6)
    xpad = 0.04 * max(grid_max, 1e-6)

    ax.set_xlim(-xpad, grid_max)
    ax.set_ylim(ymin - ypad, ymax + ypad)

    ax.legend(
        frameon=False,
        loc="upper left",
        handlelength=2.3,
    )

    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "Figure_validation_B_incremental_recovery",
        pdf=pdf,
        png=png,
        svg=svg,
    )


# =============================================================================
# OPTIONAL SUPPLEMENTARY ABSOLUTE ORACLE RECOVERY
# =============================================================================

def plot_supp_absolute_recovery(
    summary: pd.DataFrame,
    fits: pd.DataFrame,
    outdir: Path,
    *,
    pdf: bool,
    png: bool,
    svg: bool,
) -> None:
    d = summary.sort_values(["target_R", "scenario_label"]).reset_index(drop=True)

    x = pd.to_numeric(d["oracle_slope_per_week"], errors="raise").to_numpy(float)
    y = pd.to_numeric(d["recovered_slope_per_week"], errors="raise").to_numpy(float)
    lo = pd.to_numeric(d["recovered_bootstrap_q025"], errors="raise").to_numpy(float)
    hi = pd.to_numeric(d["recovered_bootstrap_q975"], errors="raise").to_numpy(float)

    fit = fits[fits["fit"].astype(str).eq("all_scenarios_ols")]
    require(len(fit) == 1, "Expected one all_scenarios_ols fit row")
    fr = fit.iloc[0]
    intercept = float(fr["intercept"])
    slope = float(fr["slope"])
    r2 = float(fr["r_squared"])

    lower = min(float(np.min(x)), float(np.min(y)), float(np.min(lo)), 0.0)
    upper = max(float(np.max(x)), float(np.max(y)), float(np.max(hi)), 0.0)
    pad = 0.08 * max(upper - lower, 1e-6)
    grid = np.linspace(lower - pad, upper + pad, 200)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    ax.plot(grid, grid, linestyle="--", linewidth=0.9, label="Identity")

    fit_label = f"Descriptive linear fit: slope = {slope:.3f}"
    if SUPP_SHOW_R2_IN_LEGEND:
        fit_label += f"; $R^2$ = {r2:.3f}"

    ax.plot(
        grid,
        intercept + slope * grid,
        linewidth=LINE_WIDTH,
        label=fit_label,
    )

    for i, row in d.iterrows():
        yy = float(row["recovered_slope_per_week"])
        qlo = float(row["recovered_bootstrap_q025"])
        qhi = float(row["recovered_bootstrap_q975"])
        yerr = np.array([[yy - qlo], [qhi - yy]])

        ax.errorbar(
            [float(row["oracle_slope_per_week"])],
            [yy],
            yerr=yerr,
            marker=marker_for(str(row["scenario_label"]), i),
            markersize=MARKER_SIZE + 0.5,
            linewidth=0,
            elinewidth=ERRORBAR_LINE_WIDTH,
            capsize=CAPSIZE,
        )
        ax.annotate(
            str(row["scenario_label"]),
            (float(row["oracle_slope_per_week"]), yy),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=ANNOTATION_SIZE,
        )

    ax.axhline(0.0, linestyle=":", linewidth=0.8)
    ax.axvline(0.0, linestyle=":", linewidth=0.8)

    ax.set_xlabel("Oracle temporal slope per week")
    ax.set_ylabel("Recovered Step-13 slope per week")
    ax.legend(frameon=False, loc="upper left")

    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "SuppFigure_validation_oracle_vs_recovered",
        pdf=pdf,
        png=png,
        svg=svg,
    )


# =============================================================================
# CLI / MAIN
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Publication plotting v2 for ClonoDynamics support-conditioned "
            "synthetic dose-response validation."
        ),
    )

    p.add_argument(
        "--summary-dir",
        type=Path,
        default=Path(
            "./results_validation/support_conditioned_validation_summary_v5_1"
        ),
        help="Directory produced by the v1 validation summarizer.",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path(
            "./results_validation/support_conditioned_validation_publication_v2"
        ),
    )

    p.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Also write the absolute oracle-versus-recovered supplementary panel.",
    )

    p.add_argument("--pdf", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--png", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--svg", action=argparse.BooleanOptionalAction, default=False)

    return p


def main() -> int:
    args = build_parser().parse_args()

    summary_dir = args.summary_dir.expanduser().resolve(strict=True)
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    configure_matplotlib()

    summary, profiles, incremental, fits = load_inputs(summary_dir)

    order = scenario_order(summary)

    # Audit: every scenario must appear in the profile table.
    require(
        set(order) == set(profiles["scenario_label"].astype(str).unique()),
        "Scenario mismatch between validation summary and temporal profiles.",
    )

    # Publication table.
    publication = build_publication_table(summary, incremental)
    publication.to_csv(
        outdir / "Table_validation_dose_response_publication.csv",
        index=False,
    )

    # Source data for reproducible figure deposition.
    a_source = (
        profiles[
            [
                "scenario_label",
                "target_R",
                "q_bio",
                "oracle_slope_per_week",
                "dt",
                "value",
                "bootstrap_q025",
                "bootstrap_q975",
                "n_subjects",
            ]
        ]
        .sort_values(["target_R", "scenario_label", "dt"])
        .reset_index(drop=True)
    )
    a_source.to_csv(
        outdir / "SourceData_Figure_validation_A_temporal_profiles.csv",
        index=False,
    )

    baseline = str(
        incremental.loc[
            incremental["paired_bootstrap_status"].astype(str).eq("baseline"),
            "scenario_label",
        ].iloc[0]
    )

    b_source = incremental[
        ~incremental["scenario_label"].astype(str).eq(baseline)
    ][
        [
            "scenario_label",
            "baseline_scenario",
            "target_R",
            "delta_oracle_slope_vs_baseline",
            "delta_recovered_slope_vs_baseline",
            "paired_bootstrap_q025_delta",
            "paired_bootstrap_q975_delta",
            "paired_bootstrap_positive_fraction",
            "incremental_recovery_fraction_point",
            "incremental_recovery_fraction_boot_q025",
            "incremental_recovery_fraction_boot_q975",
            "n_paired_bootstrap",
        ]
    ].copy()
    b_source.to_csv(
        outdir / "SourceData_Figure_validation_B_incremental_recovery.csv",
        index=False,
    )

    # Main panels.
    plot_temporal_profiles(
        summary,
        profiles,
        outdir,
        pdf=args.pdf,
        png=args.png,
        svg=args.svg,
    )
    plot_incremental_recovery(
        summary,
        incremental,
        fits,
        outdir,
        pdf=args.pdf,
        png=args.png,
        svg=args.svg,
    )

    # Optional supplementary panel.
    if args.include_supplementary:
        plot_supp_absolute_recovery(
            summary,
            fits,
            outdir,
            pdf=args.pdf,
            png=args.png,
            svg=args.svg,
        )

    # Console summary.
    fit = fits[
        fits["fit"].astype(str).eq("baseline_adjusted_through_origin")
    ].iloc[0]
    attenuation = float(fit["slope"])

    print("\n=== Validation plotting v2 complete ===")
    print("Primary figures:")
    print(" - Figure_validation_A_temporal_profiles")
    print(" - Figure_validation_B_incremental_recovery")
    print(
        f"Baseline-adjusted recovery factor: {attenuation:.6f} × oracle increment"
    )
    print("Publication table:")
    print(" - Table_validation_dose_response_publication.csv")
    print("Figure source data:")
    print(" - SourceData_Figure_validation_A_temporal_profiles.csv")
    print(" - SourceData_Figure_validation_B_incremental_recovery.csv")
    if args.include_supplementary:
        print("Supplementary figure:")
        print(" - SuppFigure_validation_oracle_vs_recovered")
    print(f"Output directory: {outdir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
