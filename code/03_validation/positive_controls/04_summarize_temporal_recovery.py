#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize oracle-calibrated doses and recovered temporal-covariance slopes.

Purpose and workflow position
-----------------------------
Run AFTER synthetic repertoires have passed through the existing production
state/transition pipeline, Step 11 and Step 13. It does not generate data, refit
states or recompute cross-replicate covariance from count files. It is the
numerical source for 05_plot_synthetic_validation.py.

Inputs
------
--oracle-dir/01_scenario_table.csv.
For every scenario selected by --scenarios:
  --dose-root/<scenario>/13-temporal_fluctuation_scaling/
    02_complete_case_subjects.csv
    04_cohort_core_metric_by_dt.csv
    05_temporal_slope_summary.csv
    10_joint_subject_bootstrap.csv (optional in code; required for paired CIs).
The selected primary estimator is equal_bin_equal_subject and metric cross_cov.
Five ordered lags (1..5), a unique primary slope row and consistent within-scenario
subject counts are required. The default scenarios are R0p00/R0p25/R0p50/R1p00.

Computed quantities
-------------------
Extract marginal recovered slopes, subject-bootstrap intervals and temporal
profiles. Use the smallest-target scenario as baseline. Compute absolute
oracle/recovered slopes, scenario-minus-baseline increments, incremental
recovery fractions and descriptive least-squares fits. If primary bootstrap
rows exist, match IDs one-to-one, calculate scenario-minus-baseline slopes and
report their percentile intervals. Oracle increments are held fixed in these
ratios; their simulation uncertainty is not bootstrapped here.

Outputs under --outdir
----------------------
00_validation_summary.csv
01_temporal_profiles.csv
02_incremental_recovery_summary.csv
03_recovery_fit_summary.csv
04_validation_metadata.json
The historical plotting functions also produce, by default, individual PDF/PNG
panels: Figure_validation_A_temporal_profiles,
Figure_validation_B_oracle_vs_recovered, and
Figure_validation_C_incremental_recovery. SVG is optional.

Avoiding duplicate figure exports
--------------------------------
For a single publication plotter, run this summarizer with
--no-png --no-pdf --no-svg, then run 05_plot_synthetic_validation.py in a separate
output directory. This suppresses figure FILES; the original code still imports
Matplotlib and constructs/closes the figures internally. No numerical or
plotting branch has been removed in this documentation-only edition.

Interpretation and audit requirements
------------------------------------
The script joins bootstrap IDs but does NOT verify equal seeds, identical
ordered subjects or identical cores ACROSS scenarios. Those conditions require
an external run-config/draw audit before interpreting intervals as paired.
Missing bootstrap files yield missing paired intervals, not invented values.
The full oracle and support-conditioned recovered covariance use different
support geometries. Raw recovered/oracle ratios are unstable near a zero oracle
slope (especially the plateau); retain them only as historical output fields.
Recovery fits are descriptive, not universal calibration or multi-seed power.

Execution and provenance
------------------------
Dependencies: Python >=3.9, NumPy, pandas and Matplotlib. Existing named output
files can be overwritten; use a fresh --outdir. Original versioned default
DIRECTORY names and output schemas are retained for compatibility.
Documentation-only edition of
10-summarize_support_conditioned_dose_response_v5_1.py; see SCRIPT_MAP.csv.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


# =============================================================================
# USER-EDITABLE AESTHETICS
# =============================================================================

FONT_FAMILY = "Arial"
FONT_SIZE = 8
AXIS_LABEL_SIZE = 8
TICK_LABEL_SIZE = 7
LEGEND_SIZE = 7
ANNOTATION_SIZE = 7

FIGURE_SIZE = (3.45, 2.75)  # inches
PNG_DPI = 600

LINE_WIDTH = 1.2
MARKER_SIZE = 4.5
CAPSIZE = 2.0

SCENARIO_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


# =============================================================================
# CONSTANTS
# =============================================================================

PRIMARY_ESTIMATOR = "equal_bin_equal_subject"
PRIMARY_METRIC = "cross_cov"

STEP13_FILES = {
    "complete": "02_complete_case_subjects.csv",
    "cohort": "04_cohort_core_metric_by_dt.csv",
    "slope": "05_temporal_slope_summary.csv",
    "bootstrap": "10_joint_subject_bootstrap.csv",
}


# =============================================================================
# HELPERS
# =============================================================================

def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_csv_required(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return pd.read_csv(path)


def bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.copy()
    x = s.astype(str).str.strip().str.lower()
    mapping = {
        "true": True, "1": True, "1.0": True, "yes": True, "y": True,
        "false": False, "0": False, "0.0": False, "no": False, "n": False,
    }
    out = x.map(mapping)
    if out.isna().any():
        bad = sorted(x[out.isna()].unique().tolist())
        raise ValueError(f"Unrecognized boolean values: {bad}")
    return out.astype(bool)


def finite(x) -> np.ndarray:
    a = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    return a[np.isfinite(a)]


def percentile_ci(values, low=0.025, high=0.975) -> Tuple[float, float, float, int]:
    x = finite(values)
    if x.size == 0:
        return np.nan, np.nan, np.nan, 0
    qlo, med, qhi = np.quantile(x, [low, 0.5, high])
    return float(qlo), float(med), float(qhi), int(x.size)


def primary_row(df: pd.DataFrame, label: str) -> pd.Series:
    required = {"estimator", "metric"}
    require(required.issubset(df.columns),
            f"{label}: missing columns {sorted(required - set(df.columns))}")
    d = df[
        df["estimator"].astype(str).eq(PRIMARY_ESTIMATOR)
        & df["metric"].astype(str).eq(PRIMARY_METRIC)
    ]
    require(len(d) == 1,
            f"{label}: expected exactly one {PRIMARY_ESTIMATOR}/{PRIMARY_METRIC} row, found {len(d)}")
    return d.iloc[0]


def configure_matplotlib() -> None:
    mpl.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": FONT_SIZE,
        "axes.labelsize": AXIS_LABEL_SIZE,
        "xtick.labelsize": TICK_LABEL_SIZE,
        "ytick.labelsize": TICK_LABEL_SIZE,
        "legend.fontsize": LEGEND_SIZE,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })


def save_figure(fig: plt.Figure, outdir: Path, stem: str, png: bool, pdf: bool, svg: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if pdf:
        fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    if png:
        fig.savefig(outdir / f"{stem}.png", dpi=PNG_DPI, bbox_inches="tight")
    if svg:
        fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def standardize_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)


# =============================================================================
# LOAD AND VALIDATE
# =============================================================================

def load_oracle(oracle_dir: Path) -> pd.DataFrame:
    path = oracle_dir / "01_scenario_table.csv"
    df = read_csv_required(path, "Oracle scenario table")

    required = {
        "scenario_label",
        "target_R",
        "sigma_fast",
        "q_bio",
        "realized_R",
        "oracle_slope_per_week",
        "oracle_V1",
        "oracle_V5",
    }
    require(required.issubset(df.columns),
            f"Oracle table missing columns: {sorted(required - set(df.columns))}")

    for c in required - {"scenario_label"}:
        df[c] = pd.to_numeric(df[c], errors="raise")

    df["scenario_label"] = df["scenario_label"].astype(str)
    df = df.sort_values(["target_R", "scenario_label"]).reset_index(drop=True)

    require(df["scenario_label"].is_unique, "Oracle scenario labels are not unique.")
    return df


def load_step13_scenario(step13_dir: Path, scenario: str) -> Dict[str, pd.DataFrame]:
    complete = read_csv_required(step13_dir / STEP13_FILES["complete"],
                                 f"{scenario} complete-case table")
    cohort = read_csv_required(step13_dir / STEP13_FILES["cohort"],
                               f"{scenario} cohort temporal profile")
    slope = read_csv_required(step13_dir / STEP13_FILES["slope"],
                              f"{scenario} temporal slope summary")

    boot_path = step13_dir / STEP13_FILES["bootstrap"]
    bootstrap = pd.read_csv(boot_path) if boot_path.is_file() else pd.DataFrame()

    return {
        "complete": complete,
        "cohort": cohort,
        "slope": slope,
        "bootstrap": bootstrap,
    }


def extract_scenario_summary(
    scenario: str,
    oracle_row: pd.Series,
    step13: Dict[str, pd.DataFrame],
) -> Tuple[dict, pd.DataFrame, pd.DataFrame]:

    # ---- Complete-case subjects
    complete = step13["complete"].copy()
    require("selected_for_primary_complete_case" in complete.columns,
            f"{scenario}: complete-case table lacks selected_for_primary_complete_case")
    selected = bool_series(complete["selected_for_primary_complete_case"])
    selected_subjects = sorted(
        pd.to_numeric(complete.loc[selected, "subject"], errors="raise").astype(int).tolist()
    )

    # ---- Primary slope
    slope = step13["slope"].copy()
    sr = primary_row(slope, f"{scenario} slope summary")

    slope_required = {
        "observed_slope_per_week",
        "bootstrap_median_slope",
        "bootstrap_q025",
        "bootstrap_q975",
        "bootstrap_positive_slope_fraction",
        "bootstrap_negative_slope_fraction",
        "n_complete_case_subjects",
        "n_core_bins",
        "core_x_min",
        "core_x_max",
    }
    require(slope_required.issubset(slope.columns),
            f"{scenario}: slope table missing {sorted(slope_required - set(slope.columns))}")

    # ---- Primary temporal profile
    cohort = step13["cohort"].copy()
    required_profile = {
        "estimator", "metric", "dt", "value",
        "bootstrap_median", "bootstrap_q025", "bootstrap_q975",
        "n_subjects",
    }
    require(required_profile.issubset(cohort.columns),
            f"{scenario}: cohort table missing {sorted(required_profile - set(cohort.columns))}")

    profile = cohort[
        cohort["estimator"].astype(str).eq(PRIMARY_ESTIMATOR)
        & cohort["metric"].astype(str).eq(PRIMARY_METRIC)
    ].copy()
    profile["dt"] = pd.to_numeric(profile["dt"], errors="raise").astype(int)
    profile = profile.sort_values("dt").reset_index(drop=True)
    require(profile["dt"].tolist() == [1, 2, 3, 4, 5],
            f"{scenario}: expected dt 1..5, found {profile['dt'].tolist()}")

    profile.insert(0, "scenario_label", scenario)
    profile.insert(1, "target_R", float(oracle_row["target_R"]))
    profile.insert(2, "q_bio", float(oracle_row["q_bio"]))
    profile.insert(3, "oracle_slope_per_week", float(oracle_row["oracle_slope_per_week"]))

    # ---- Primary bootstrap slope distribution, when available
    boot = step13["bootstrap"].copy()
    primary_boot = pd.DataFrame()
    if not boot.empty:
        boot_required = {"bootstrap", "estimator", "metric", "slope_per_week"}
        require(boot_required.issubset(boot.columns),
                f"{scenario}: bootstrap table missing {sorted(boot_required - set(boot.columns))}")
        primary_boot = boot[
            boot["estimator"].astype(str).eq(PRIMARY_ESTIMATOR)
            & boot["metric"].astype(str).eq(PRIMARY_METRIC)
        ][["bootstrap", "slope_per_week"]].copy()
        primary_boot["bootstrap"] = pd.to_numeric(primary_boot["bootstrap"], errors="raise").astype(int)
        primary_boot["slope_per_week"] = pd.to_numeric(
            primary_boot["slope_per_week"], errors="coerce"
        )
        require(primary_boot["bootstrap"].is_unique,
                f"{scenario}: duplicate primary bootstrap IDs")
        primary_boot["scenario_label"] = scenario

    summary = {
        "scenario_label": scenario,
        "target_R": float(oracle_row["target_R"]),
        "sigma_fast": float(oracle_row["sigma_fast"]),
        "q_bio": float(oracle_row["q_bio"]),
        "realized_R": float(oracle_row["realized_R"]),
        "oracle_V1": float(oracle_row["oracle_V1"]),
        "oracle_V5": float(oracle_row["oracle_V5"]),
        "oracle_slope_per_week": float(oracle_row["oracle_slope_per_week"]),
        "recovered_intercept_at_dt1": float(sr["observed_intercept_at_dt1"]),
        "recovered_slope_per_week": float(sr["observed_slope_per_week"]),
        "recovered_bootstrap_median_slope": float(sr["bootstrap_median_slope"]),
        "recovered_bootstrap_q025": float(sr["bootstrap_q025"]),
        "recovered_bootstrap_q975": float(sr["bootstrap_q975"]),
        "bootstrap_positive_slope_fraction": float(sr["bootstrap_positive_slope_fraction"]),
        "bootstrap_negative_slope_fraction": float(sr["bootstrap_negative_slope_fraction"]),
        "positive_slope_resolved": bool(float(sr["bootstrap_q025"]) > 0),
        "negative_slope_resolved": bool(float(sr["bootstrap_q975"]) < 0),
        "n_complete_case_subjects": int(sr["n_complete_case_subjects"]),
        "complete_case_subjects": ";".join(map(str, selected_subjects)),
        "n_core_bins": int(sr["n_core_bins"]),
        "core_x_min": float(sr["core_x_min"]),
        "core_x_max": float(sr["core_x_max"]),
    }

    # Cross-check selected subject count
    require(summary["n_complete_case_subjects"] == len(selected_subjects),
            f"{scenario}: slope summary and complete-case table disagree on subject count")

    return summary, profile, primary_boot


# =============================================================================
# RECOVERY STATISTICS
# =============================================================================

def add_recovery_columns(summary: pd.DataFrame) -> pd.DataFrame:
    d = summary.sort_values(["target_R", "scenario_label"]).reset_index(drop=True).copy()

    baseline = d.iloc[0]
    baseline_oracle = float(baseline["oracle_slope_per_week"])
    baseline_recovered = float(baseline["recovered_slope_per_week"])

    d["raw_recovery_ratio"] = np.where(
        np.abs(d["oracle_slope_per_week"]) > 1e-15,
        d["recovered_slope_per_week"] / d["oracle_slope_per_week"],
        np.nan,
    )

    d["delta_oracle_slope_vs_baseline"] = (
        d["oracle_slope_per_week"] - baseline_oracle
    )
    d["delta_recovered_slope_vs_baseline"] = (
        d["recovered_slope_per_week"] - baseline_recovered
    )

    denom = d["delta_oracle_slope_vs_baseline"].to_numpy(float)
    numer = d["delta_recovered_slope_vs_baseline"].to_numpy(float)
    d["incremental_recovery_fraction"] = np.where(
        np.abs(denom) > 1e-15,
        numer / denom,
        np.nan,
    )

    return d


def recovery_fit_summary(summary: pd.DataFrame) -> pd.DataFrame:
    x = summary["oracle_slope_per_week"].to_numpy(float)
    y = summary["recovered_slope_per_week"].to_numpy(float)

    slope, intercept = np.polyfit(x, y, 1)
    yhat = intercept + slope * x
    sse = float(np.sum((y - yhat) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - sse / sst) if sst > 0 else np.nan

    p = summary[np.abs(summary["delta_oracle_slope_vs_baseline"]) > 1e-15].copy()
    dx = p["delta_oracle_slope_vs_baseline"].to_numpy(float)
    dy = p["delta_recovered_slope_vs_baseline"].to_numpy(float)

    attenuation = float(np.sum(dx * dy) / np.sum(dx * dx))
    fractions = p["incremental_recovery_fraction"].to_numpy(float)

    rows = [
        {
            "fit": "all_scenarios_ols",
            "n_points": len(summary),
            "intercept": float(intercept),
            "slope": float(slope),
            "r_squared": r2,
            "description": "Recovered Step-13 slope versus oracle slope",
        },
        {
            "fit": "baseline_adjusted_through_origin",
            "n_points": len(p),
            "intercept": 0.0,
            "slope": attenuation,
            "r_squared": np.nan,
            "description": "Delta recovered versus delta oracle relative to zero-accumulation baseline",
        },
        {
            "fit": "incremental_recovery_fraction_summary",
            "n_points": len(p),
            "intercept": np.nan,
            "slope": float(np.mean(fractions)),
            "r_squared": np.nan,
            "description": (
                f"Mean incremental recovery fraction; min={np.min(fractions):.6f}; "
                f"max={np.max(fractions):.6f}"
            ),
        },
    ]
    return pd.DataFrame(rows)


def paired_incremental_bootstrap(
    summary: pd.DataFrame,
    boot_by_scenario: Dict[str, pd.DataFrame],
) -> pd.DataFrame:

    baseline = summary.sort_values(["target_R", "scenario_label"]).iloc[0]
    baseline_label = str(baseline["scenario_label"])
    baseline_boot = boot_by_scenario.get(baseline_label, pd.DataFrame())

    rows = []

    for r in summary.sort_values(["target_R", "scenario_label"]).itertuples(index=False):
        scenario = str(r.scenario_label)
        delta_oracle = float(r.delta_oracle_slope_vs_baseline)
        delta_point = float(r.delta_recovered_slope_vs_baseline)

        if scenario == baseline_label:
            rows.append({
                "scenario_label": scenario,
                "baseline_scenario": baseline_label,
                "target_R": float(r.target_R),
                "delta_oracle_slope_vs_baseline": 0.0,
                "delta_recovered_slope_vs_baseline": 0.0,
                "paired_bootstrap_median_delta": 0.0,
                "paired_bootstrap_q025_delta": 0.0,
                "paired_bootstrap_q975_delta": 0.0,
                "paired_bootstrap_positive_fraction": np.nan,
                "incremental_recovery_fraction_point": np.nan,
                "incremental_recovery_fraction_boot_median": np.nan,
                "incremental_recovery_fraction_boot_q025": np.nan,
                "incremental_recovery_fraction_boot_q975": np.nan,
                "n_paired_bootstrap": 0,
                "paired_bootstrap_status": "baseline",
            })
            continue

        this_boot = boot_by_scenario.get(scenario, pd.DataFrame())

        if baseline_boot.empty or this_boot.empty:
            rows.append({
                "scenario_label": scenario,
                "baseline_scenario": baseline_label,
                "target_R": float(r.target_R),
                "delta_oracle_slope_vs_baseline": delta_oracle,
                "delta_recovered_slope_vs_baseline": delta_point,
                "paired_bootstrap_median_delta": np.nan,
                "paired_bootstrap_q025_delta": np.nan,
                "paired_bootstrap_q975_delta": np.nan,
                "paired_bootstrap_positive_fraction": np.nan,
                "incremental_recovery_fraction_point": (
                    delta_point / delta_oracle if abs(delta_oracle) > 1e-15 else np.nan
                ),
                "incremental_recovery_fraction_boot_median": np.nan,
                "incremental_recovery_fraction_boot_q025": np.nan,
                "incremental_recovery_fraction_boot_q975": np.nan,
                "n_paired_bootstrap": 0,
                "paired_bootstrap_status": "bootstrap_file_missing",
            })
            continue

        m = baseline_boot[["bootstrap", "slope_per_week"]].merge(
            this_boot[["bootstrap", "slope_per_week"]],
            on="bootstrap",
            how="inner",
            suffixes=("_baseline", "_scenario"),
            validate="one_to_one",
        )
        require(len(m) > 0, f"{scenario}: no matched bootstrap IDs with baseline")

        m["delta"] = m["slope_per_week_scenario"] - m["slope_per_week_baseline"]
        x = finite(m["delta"])
        qlo, med, qhi, n = percentile_ci(x)

        if abs(delta_oracle) > 1e-15:
            rr = x / delta_oracle
            rqlo, rmed, rqhi, _ = percentile_ci(rr)
            point_ratio = delta_point / delta_oracle
        else:
            rqlo = rmed = rqhi = point_ratio = np.nan

        rows.append({
            "scenario_label": scenario,
            "baseline_scenario": baseline_label,
            "target_R": float(r.target_R),
            "delta_oracle_slope_vs_baseline": delta_oracle,
            "delta_recovered_slope_vs_baseline": delta_point,
            "paired_bootstrap_median_delta": med,
            "paired_bootstrap_q025_delta": qlo,
            "paired_bootstrap_q975_delta": qhi,
            "paired_bootstrap_positive_fraction": float(np.mean(x > 0)) if n else np.nan,
            "incremental_recovery_fraction_point": point_ratio,
            "incremental_recovery_fraction_boot_median": rmed,
            "incremental_recovery_fraction_boot_q025": rqlo,
            "incremental_recovery_fraction_boot_q975": rqhi,
            "n_paired_bootstrap": n,
            "paired_bootstrap_status": "paired_by_bootstrap_id",
        })

    return pd.DataFrame(rows)


# =============================================================================
# FIGURES
# =============================================================================

def figure_temporal_profiles(
    profiles: pd.DataFrame,
    scenario_order: List[str],
    outdir: Path,
    png: bool,
    pdf: bool,
    svg: bool,
) -> None:

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    for i, scenario in enumerate(scenario_order):
        d = profiles[profiles["scenario_label"].eq(scenario)].sort_values("dt")
        require(len(d) == 5, f"{scenario}: expected 5 temporal-profile rows")

        x = d["dt"].to_numpy(float)
        y = d["value"].to_numpy(float)
        lo = d["bootstrap_q025"].to_numpy(float)
        hi = d["bootstrap_q975"].to_numpy(float)
        yerr = np.vstack([y - lo, hi - y])

        target_R = float(d["target_R"].iloc[0])

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            marker=SCENARIO_MARKERS[i % len(SCENARIO_MARKERS)],
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            capsize=CAPSIZE,
            label=f"{scenario}  (target R={target_R:.2f})",
        )

    ax.axhline(0.0, linewidth=0.8, linestyle="--")
    ax.set_xlabel(r"Temporal lag, $\Delta t$ (weeks)")
    ax.set_ylabel("Cross-replicate covariance")
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.legend(frameon=False)
    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig, outdir, "Figure_validation_A_temporal_profiles",
        png=png, pdf=pdf, svg=svg
    )


def figure_oracle_vs_recovered(
    summary: pd.DataFrame,
    fit_summary: pd.DataFrame,
    outdir: Path,
    png: bool,
    pdf: bool,
    svg: bool,
) -> None:

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    d = summary.sort_values(["target_R", "scenario_label"]).reset_index(drop=True)

    x = d["oracle_slope_per_week"].to_numpy(float)
    y = d["recovered_slope_per_week"].to_numpy(float)
    lo = d["recovered_bootstrap_q025"].to_numpy(float)
    hi = d["recovered_bootstrap_q975"].to_numpy(float)

    xmin = min(float(np.min(x)), float(np.min(y)), 0.0)
    xmax = max(float(np.max(x)), float(np.max(y)), 0.0)
    pad = 0.08 * max(xmax - xmin, 1e-6)
    grid = np.linspace(xmin - pad, xmax + pad, 200)

    # Identity
    ax.plot(grid, grid, linestyle="--", linewidth=0.9, label="Identity")

    # OLS fit
    fr = fit_summary[fit_summary["fit"].eq("all_scenarios_ols")].iloc[0]
    fit_y = float(fr["intercept"]) + float(fr["slope"]) * grid
    ax.plot(
        grid,
        fit_y,
        linewidth=LINE_WIDTH,
        label=(
            f"OLS: y={float(fr['intercept']):+.4f}"
            f"+{float(fr['slope']):.3f}x; "
            f"$R^2$={float(fr['r_squared']):.3f}"
        ),
    )

    for i, r in enumerate(d.itertuples(index=False)):
        yerr = np.array([
            [r.recovered_slope_per_week - r.recovered_bootstrap_q025],
            [r.recovered_bootstrap_q975 - r.recovered_slope_per_week],
        ])
        ax.errorbar(
            [r.oracle_slope_per_week],
            [r.recovered_slope_per_week],
            yerr=yerr,
            marker=SCENARIO_MARKERS[i % len(SCENARIO_MARKERS)],
            markersize=MARKER_SIZE + 0.5,
            linewidth=0,
            elinewidth=LINE_WIDTH,
            capsize=CAPSIZE,
        )
        ax.annotate(
            str(r.scenario_label),
            (r.oracle_slope_per_week, r.recovered_slope_per_week),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=ANNOTATION_SIZE,
        )

    ax.axhline(0.0, linewidth=0.7, linestyle=":")
    ax.axvline(0.0, linewidth=0.7, linestyle=":")
    ax.set_xlabel("Oracle temporal slope per week")
    ax.set_ylabel("Recovered Step-13 slope per week")
    ax.legend(frameon=False)
    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig, outdir, "Figure_validation_B_oracle_vs_recovered",
        png=png, pdf=pdf, svg=svg
    )


def figure_incremental_recovery(
    summary: pd.DataFrame,
    incremental: pd.DataFrame,
    fit_summary: pd.DataFrame,
    outdir: Path,
    png: bool,
    pdf: bool,
    svg: bool,
) -> None:

    d = summary[
        np.abs(summary["delta_oracle_slope_vs_baseline"]) > 1e-15
    ].sort_values("target_R").copy()

    inc = incremental.set_index("scenario_label")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    x = d["delta_oracle_slope_vs_baseline"].to_numpy(float)
    y = d["delta_recovered_slope_vs_baseline"].to_numpy(float)

    xmax = max(float(np.max(x)), float(np.max(y)), 1e-6)
    grid = np.linspace(0.0, xmax * 1.08, 200)

    ax.plot(grid, grid, linestyle="--", linewidth=0.9, label="Identity")

    fr = fit_summary[
        fit_summary["fit"].eq("baseline_adjusted_through_origin")
    ].iloc[0]
    attenuation = float(fr["slope"])
    ax.plot(
        grid,
        attenuation * grid,
        linewidth=LINE_WIDTH,
        label=f"Through-origin fit: {attenuation:.3f}× oracle increment",
    )

    for i, r in enumerate(d.itertuples(index=False)):
        ir = inc.loc[str(r.scenario_label)]

        if (
            str(ir["paired_bootstrap_status"]) == "paired_by_bootstrap_id"
            and np.isfinite(ir["paired_bootstrap_q025_delta"])
            and np.isfinite(ir["paired_bootstrap_q975_delta"])
        ):
            lo = float(ir["paired_bootstrap_q025_delta"])
            hi = float(ir["paired_bootstrap_q975_delta"])
            yy = float(r.delta_recovered_slope_vs_baseline)
            yerr = np.array([[yy - lo], [hi - yy]])
        else:
            yerr = None

        ax.errorbar(
            [r.delta_oracle_slope_vs_baseline],
            [r.delta_recovered_slope_vs_baseline],
            yerr=yerr,
            marker=SCENARIO_MARKERS[(i + 1) % len(SCENARIO_MARKERS)],
            markersize=MARKER_SIZE + 0.5,
            linewidth=0,
            elinewidth=LINE_WIDTH,
            capsize=CAPSIZE,
        )

        ax.annotate(
            f"{r.scenario_label}\n{100*r.incremental_recovery_fraction:.1f}%",
            (r.delta_oracle_slope_vs_baseline, r.delta_recovered_slope_vs_baseline),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=ANNOTATION_SIZE,
        )

    ax.axhline(0.0, linewidth=0.7, linestyle=":")
    ax.set_xlabel(r"Oracle slope increment vs R0p00")
    ax.set_ylabel(r"Recovered slope increment vs R0p00")
    ax.legend(frameon=False)
    standardize_axes(ax)
    fig.tight_layout()

    save_figure(
        fig, outdir, "Figure_validation_C_incremental_recovery",
        png=png, pdf=pdf, svg=svg
    )


# =============================================================================
# MAIN
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Summarize and plot ClonoDynamics empirical-support-conditioned "
            "synthetic dose-response validation."
        ),
    )

    p.add_argument(
        "--oracle-dir",
        type=Path,
        default=Path("./results_validation/qbio_oracle_v5"),
    )
    p.add_argument(
        "--dose-root",
        type=Path,
        default=Path("./results_validation/support_conditioned_dose_v5_1"),
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("./results_validation/support_conditioned_validation_summary_v5_1"),
    )
    p.add_argument(
        "--scenarios",
        nargs="*",
        default=["R0p00", "R0p25", "R0p50", "R1p00"],
    )

    p.add_argument("--png", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--pdf", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--svg", action=argparse.BooleanOptionalAction, default=False)

    return p


def main() -> int:
    args = build_parser().parse_args()

    oracle_dir = args.oracle_dir.expanduser().resolve()
    dose_root = args.dose_root.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    configure_matplotlib()

    oracle = load_oracle(oracle_dir)

    requested = [str(x) for x in args.scenarios]
    missing_oracle = sorted(set(requested) - set(oracle["scenario_label"]))
    require(not missing_oracle,
            f"Requested scenarios absent from oracle table: {missing_oracle}")

    oracle = oracle[oracle["scenario_label"].isin(requested)].copy()
    oracle = oracle.sort_values(["target_R", "scenario_label"]).reset_index(drop=True)

    scenario_rows = []
    profile_frames = []
    boot_by_scenario: Dict[str, pd.DataFrame] = {}

    for orow in oracle.itertuples(index=False):
        scenario = str(orow.scenario_label)
        step13_dir = dose_root / scenario / "13-temporal_fluctuation_scaling"

        step13 = load_step13_scenario(step13_dir, scenario)
        summary_row, profile, pboot = extract_scenario_summary(
            scenario,
            pd.Series(orow._asdict()),
            step13,
        )

        scenario_rows.append(summary_row)
        profile_frames.append(profile)
        boot_by_scenario[scenario] = pboot

    summary = pd.DataFrame(scenario_rows)
    summary = add_recovery_columns(summary)

    profiles = pd.concat(profile_frames, ignore_index=True)
    fits = recovery_fit_summary(summary)
    incremental = paired_incremental_bootstrap(summary, boot_by_scenario)

    # Add paired-bootstrap incremental columns into the main summary.
    merge_cols = [
        "scenario_label",
        "paired_bootstrap_q025_delta",
        "paired_bootstrap_q975_delta",
        "paired_bootstrap_positive_fraction",
        "incremental_recovery_fraction_boot_median",
        "incremental_recovery_fraction_boot_q025",
        "incremental_recovery_fraction_boot_q975",
        "n_paired_bootstrap",
        "paired_bootstrap_status",
    ]
    summary = summary.merge(
        incremental[merge_cols],
        on="scenario_label",
        how="left",
        validate="one_to_one",
    )

    # Outputs
    summary.to_csv(outdir / "00_validation_summary.csv", index=False)
    profiles.to_csv(outdir / "01_temporal_profiles.csv", index=False)
    incremental.to_csv(outdir / "02_incremental_recovery_summary.csv", index=False)
    fits.to_csv(outdir / "03_recovery_fit_summary.csv", index=False)

    metadata = {
        "primary_estimator": PRIMARY_ESTIMATOR,
        "primary_metric": PRIMARY_METRIC,
        "oracle_dir": str(oracle_dir),
        "dose_root": str(dose_root),
        "scenario_order": summary.sort_values("target_R")["scenario_label"].tolist(),
        "baseline_scenario": str(summary.sort_values("target_R").iloc[0]["scenario_label"]),
        "incremental_recovery_definition": (
            "(recovered_slope_s - recovered_slope_baseline) / "
            "(oracle_slope_s - oracle_slope_baseline)"
        ),
        "paired_bootstrap_note": (
            "Paired by Step-13 bootstrap ID when 10_joint_subject_bootstrap.csv "
            "is available for baseline and target scenarios."
        ),
    }
    (outdir / "04_validation_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    scenario_order = summary.sort_values("target_R")["scenario_label"].tolist()

    figure_temporal_profiles(
        profiles, scenario_order, outdir,
        png=args.png, pdf=args.pdf, svg=args.svg,
    )
    figure_oracle_vs_recovered(
        summary, fits, outdir,
        png=args.png, pdf=args.pdf, svg=args.svg,
    )
    figure_incremental_recovery(
        summary, incremental, fits, outdir,
        png=args.png, pdf=args.pdf, svg=args.svg,
    )

    # Console summary
    print("\n=== ClonoDynamics validation dose-response summary ===")
    display_cols = [
        "scenario_label",
        "target_R",
        "q_bio",
        "oracle_slope_per_week",
        "recovered_slope_per_week",
        "recovered_bootstrap_q025",
        "recovered_bootstrap_q975",
        "delta_oracle_slope_vs_baseline",
        "delta_recovered_slope_vs_baseline",
        "incremental_recovery_fraction",
    ]
    with pd.option_context(
        "display.max_columns", None,
        "display.width", 180,
        "display.float_format", lambda x: f"{x:.6f}",
    ):
        print(summary[display_cols].to_string(index=False))

    ols = fits[fits["fit"].eq("all_scenarios_ols")].iloc[0]
    incfit = fits[fits["fit"].eq("baseline_adjusted_through_origin")].iloc[0]

    print(
        f"\nOracle -> recovered OLS: "
        f"y = {ols['intercept']:+.6f} + {ols['slope']:.6f} x; "
        f"R^2 = {ols['r_squared']:.6f}"
    )
    print(
        f"Baseline-adjusted through-origin attenuation factor: "
        f"{incfit['slope']:.6f}"
    )
    print(f"\nWrote outputs to: {outdir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
