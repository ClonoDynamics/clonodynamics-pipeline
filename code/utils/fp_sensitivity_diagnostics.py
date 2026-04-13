#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fp_sensitivity_diagnostics_v6.py

Publication-oriented sensitivity diagnostics for Fokker–Planck analyses across:
    1) detectability threshold alpha
    2) transition / observability class (e.g. TT vs ALL, weighted vs unweighted)
    3) fitted model (e.g. hinge_plateau_smooth, OU_linear)

This version is designed for paper-ready plotting and adds a crucial distinction:

It also adds BD OVERLAY plots, allowing drift and diffusion from multiple loaded
conditions (obs_class / alpha combinations from the manifest or auto-discovery)
to be drawn on the same axes for direct comparison.

    DATA FILTERING RANGE
        --x_min / --x_max
        These remove rows from the data before metrics and plotting.

    PLOTTING RANGE
        --x_plot_min / --x_plot_max
        These only change the visible x-axis limits.
        They DO NOT change the data used for summaries, integrals, ratios, or diagnostics.

This distinction is especially useful when:
    - you want full-range numerical analysis
    - but cleaner “bulk-only” visualization in the figure

===============================================================================
INPUT MODES
===============================================================================

A) Manifest mode
----------------
Pass:

    --manifest_csv /path/to/manifest.csv

Required columns:
    obs_class, alpha, grid_csv, summary_json

Optional column:
    model

If the model column is absent, the script assigns:
    model = "unknown_model"

Use manifest mode when you want:
    - explicit reproducibility
    - manual curation of conditions
    - stable ordering independent of directory layout

B) Auto-discovery mode
----------------------
Pass:

    --input_root /path/to/results_root

The script recursively scans the root directory for:
    *.grid.csv

For each grid file, it looks for the paired summary JSON:
    same_name.summary.json
derived from:
    same_name.grid.csv

It then infers:
    - alpha
    - obs_class
    - model

Inference sources:
    1) summary_json, when available
    2) directory names and filenames

Examples of recognized patterns:
    alpha:
        p_005
        p_0_02
        alpha_0.02
        p0.05

    obs_class:
        TT
        ALL
        TF
        FT
        FF
        ALL_unweighted
        ALL_w_clone
        TT_unweighted
        TT_w_clone

    model:
        OU_linear
        hinge_plateau_smooth
        tanh_saturating
        rational_saturating

IMPORTANT:
    Use either --manifest_csv OR --input_root, but not both.

===============================================================================
EXPECTED GRID CONTENT
===============================================================================

Required columns (with common-name fallbacks):
    x
    p_emp
    p_closed
    J
    S

Optional columns (used if present):
    ratio_pemp_over_pclosed
        If absent, ratio is computed as p_emp / p_closed

    b / drift
        Used for inferred-function plots

    D_final / diffusion
        Used for inferred-function plots

Column-name resolution is flexible and supports common alternatives.

===============================================================================
WHAT THE SCRIPT CAN DRAW
===============================================================================

1) FP diagnostics
-----------------
Main Fokker–Planck diagnostics:
    - empirical vs stationary density
    - density ratio p_emp / p_closed
    - probability current J(x)
    - source term S(x)

2) BD diagnostics
-----------------
If drift and diffusion columns are available, and if requested with --also_bd,
the script can also plot:
    - drift b(x)
    - diffusion D(x)

These are useful for:
    - TT vs ALL controls
    - alpha sensitivity
    - OU vs hinge comparisons
    - supplementary figures showing the inferred ingredients of the FP model

===============================================================================
OUTPUTS
===============================================================================

For each model, the script writes:

Tables:
    fp_sensitivity_summary__<model>.csv
    fp_sensitivity_summary_pretty__<model>.csv
    fp_sensitivity_report__<model>.txt

Figures:
    fig_fp_sensitivity_metrics_by_obsclass__<model>.<fmt>
    fig_fp_sensitivity_multipanel_4xN__<model>.<fmt>

Optionally, when --also_bd is active and BD columns exist:
    fig_bd_sensitivity_multipanel_2xN__<model>.<fmt>

Single-condition figures are also written unless disabled:
    fig_fp_pdf__<model>__<condition>.<fmt>
    fig_fp_ratio__<model>__<condition>.<fmt>
    fig_fp_flux__<model>__<condition>.<fmt>
    fig_fp_source__<model>__<condition>.<fmt>

Optionally:
    fig_bd_drift__<model>__<condition>.<fmt>
    fig_bd_diffusion__<model>__<condition>.<fmt>

Overlay figures (new in v6, when requested):
    fig_bd_overlay_drift__<model>.<fmt>
    fig_bd_overlay_diffusion__<model>.<fmt>

===============================================================================
PLOT SELECTION MODES
===============================================================================

By default
----------
The script writes:
    - summary tables
    - metrics plot
    - single-condition FP plots
    - FP multipanel
    - optionally BD single plots and BD multipanel

--only_multipanel
-----------------
If this flag is used, the script suppresses all single-condition plots and
writes only multipanel figures:
    - FP multipanel
    - optionally BD multipanel (if --also_bd is active)

--only_fp_multipanel
--------------------
If this flag is used, the script writes only:
    - the FP multipanel
and skips:
    - all single-condition plots
    - BD multipanel
    - BD single-condition plots

--only_bd_multipanel
--------------------
If this flag is used, the script writes only:
    - the BD multipanel
and skips:
    - all single-condition plots
    - the FP multipanel
    - FP single-condition plots

MUTUAL EXCLUSIVITY
------------------
The following flags are mutually exclusive:
    --only_multipanel
    --only_fp_multipanel
    --only_bd_multipanel

===============================================================================
DETAILED DESCRIPTION OF OPTIONS
===============================================================================

INPUT OPTIONS
-------------
--manifest_csv
    Path to manifest CSV with columns:
        obs_class, alpha, grid_csv, summary_json
    optional:
        model

--input_root
    Root directory for recursive auto-discovery of grid/summary pairs.

--out_dir
    Output directory where all tables and figures will be written.

--fmt
    Figure format:
        pdf
        png

ORDERING OPTIONS
----------------
--obs_order
    Optional comma-separated order for observability classes / conditions.
    Example:
        --obs_order ALL_unweighted,ALL_w_clone,TT_unweighted,TT_w_clone

--alpha_order
    Optional comma-separated order for alpha values.
    Example:
        --alpha_order 0.005,0.01,0.02,0.05

--model_order
    Optional comma-separated order for models.
    Example:
        --model_order OU_linear,hinge_plateau_smooth

--facet_mode
    Controls column order inside each model-specific multipanel.
    Allowed values:
        obs_then_alpha
        alpha_then_obs

X-AXIS OPTIONS
--------------
--x_min
    Global lower x bound applied during data loading.
    This FILTERS the data:
        rows with x < x_min are removed before metrics and plotting.

--x_max
    Global upper x bound applied during data loading.
    This FILTERS the data:
        rows with x > x_max are removed before metrics and plotting.

--x_mode
    Controls x-limits in the final plots.
    Allowed values:
        common
        individual

    common:
        all panels for a given model share the same x-limits
        recommended for direct comparison

    individual:
        each panel keeps its own natural x range

--x_plot_min
    Lower x limit for PLOTTING ONLY.
    Does not filter the data.
    Useful for visually focusing on the bulk while keeping full-range analysis.

--x_plot_max
    Upper x limit for PLOTTING ONLY.
    Does not filter the data.

Practical distinction:
    --x_min / --x_max       change the analysis domain
    --x_plot_min / x_plot_max change only the visible plotting window

SMOOTHING OPTIONS
-----------------
--ratio_smooth_sigma_pts
    Gaussian smoothing sigma (in grid points) for the density ratio.

--flux_smooth_sigma_pts
    Gaussian smoothing sigma (in grid points) for both J(x) and S(x).

--drift_smooth_sigma_pts
    Optional additional smoothing sigma for drift b(x) in BD plots.

--diffusion_smooth_sigma_pts
    Optional additional smoothing sigma for diffusion D(x) in BD plots.

Y-LIMIT OPTIONS
---------------
You may manually set y-limits for any row:

    --pdf_ymin / --pdf_ymax
    --ratio_ymin / --ratio_ymax
    --flux_ymin / --flux_ymax
    --source_ymin / --source_ymax
    --drift_ymin / --drift_ymax
    --diffusion_ymin / --diffusion_ymax

If omitted, robust automatic limits are computed.

AUTOMATIC Y-LIMIT OPTIONS
-------------------------
--auto_qlo
    Lower quantile for robust automatic y-limits.

--auto_qhi
    Upper quantile for robust automatic y-limits.

--auto_pad_frac
    Padding fraction added around the robust automatic limits.

--flux_symmetric
    Forces J(x) limits to be symmetric around zero.

--source_symmetric
    Forces S(x) limits to be symmetric around zero.

FIGURE SIZE AND LINE STYLE OPTIONS
----------------------------------
--single_width / --single_height
    Size of single-condition plots in inches.

--multi_width / --multi_height
    Size of FP multipanel in inches.

--bd_multi_height
    Height of BD multipanel in inches.

--dpi
    Output resolution.

--line_w_pdf
    Line width for density curves.

--line_w_ratio
    Line width for ratio curves.

--line_w_flux
    Line width for J(x).

--line_w_source
    Line width for S(x).

--line_w_bd
    Line width for drift and diffusion plots.

--refline_w
    Line width for horizontal reference lines.

--refline_alpha
    Transparency for horizontal reference lines.

DISPLAY / STYLING FLAGS
-----------------------
--also_bd
    If BD columns are available, also produce drift/diffusion plots.

--show_legends
    Show legends in single-condition plots where appropriate.

--panel_titles
    Show per-panel condition titles.

--figure_title
    Show figure-level title in multipanels.

VISUAL GUIDE BANDS
------------------
--show_ratio_band
    Shade a band around ratio=1.

--ratio_band_halfwidth
    Half-width of the shaded ratio band.

--show_zero_band
    Shade a band around zero in J(x) and S(x).

--flux_zero_band_halfwidth
    Half-width of the shaded band around J(x)=0.

--source_zero_band_halfwidth
    Half-width of the shaded band around S(x)=0.

PLOT SELECTION FLAGS
--------------------
--only_multipanel
    Produce only multipanel figures (FP multipanel, and BD multipanel if
    --also_bd is active).

--only_fp_multipanel
    Produce only the FP multipanel.

--only_bd_multipanel
    Produce only the BD multipanel.

BD OVERLAY OPTIONS
------------------
--bd_overlay
    Produce overlay plots for drift and diffusion, with all loaded conditions from the
    manifest (or auto-discovered set) drawn on the same axes for each model.

--bd_overlay_by
    Group overlay curves by:
        condition   -> one curve per loaded condition (default)
        obs_class   -> overlays all conditions but legend emphasizes obs_class labels

--bd_overlay_legend_outside
    Place overlay legend outside the plotting area on the right.

===============================================================================
EXAMPLES
===============================================================================

1) Main paper style: only FP multipanel, bulk-focused view
----------------------------------------------------------
python fp_sensitivity_diagnostics_v5.py \
  --input_root ./results \
  --out_dir ./fp_sensitivity \
  --obs_order ALL_unweighted,ALL_w_clone,TT_unweighted,TT_w_clone \
  --alpha_order 0.02 \
  --model_order hinge_plateau_smooth \
  --facet_mode obs_then_alpha \
  --panel_titles \
  --figure_title \
  --show_ratio_band \
  --show_zero_band \
  --flux_symmetric \
  --source_symmetric \
  --x_plot_min -10.5 \
  --x_plot_max -6.0 \
  --only_fp_multipanel

2) Supplementary inferred functions only
----------------------------------------
python fp_sensitivity_diagnostics_v5.py \
  --input_root ./results \
  --out_dir ./fp_sensitivity_bd \
  --obs_order ALL_unweighted,ALL_w_clone,TT_unweighted,TT_w_clone \
  --alpha_order 0.02 \
  --model_order hinge_plateau_smooth \
  --also_bd \
  --panel_titles \
  --figure_title \
  --x_plot_min -10.5 \
  --x_plot_max -6.0 \
  --only_bd_multipanel

3) Full-range numerical analysis but zoomed plotting
----------------------------------------------------
python fp_sensitivity_diagnostics_v5.py \
  --manifest_csv ./fp_sensitivity_manifest.csv \
  --out_dir ./fp_sensitivity_full \
  --obs_order TT,ALL \
  --alpha_order 0.005,0.01,0.02,0.05 \
  --model_order OU_linear,hinge_plateau_smooth \
  --also_bd \
  --x_plot_min -10.5 \
  --x_plot_max -6.0
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


KNOWN_MODELS = [
    "OU_linear",
    "hinge_plateau_smooth",
    "tanh_saturating",
    "rational_saturating",
]
KNOWN_OBS = ["TT", "ALL", "TF", "FT", "FF"]


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Publication-oriented FP sensitivity diagnostics across alpha, obs_class, and model."
    )

    ap.add_argument("--manifest_csv", default=None,
                    help="CSV manifest with columns: obs_class, alpha, grid_csv, summary_json, optional model")
    ap.add_argument("--input_root", default=None,
                    help="Root directory to auto-discover *.grid.csv / *.summary.json pairs")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--fmt", default="pdf", choices=["pdf", "png"], help="Figure format")

    ap.add_argument("--obs_order", default="", help="Optional comma-separated obs_class order")
    ap.add_argument("--alpha_order", default="", help="Optional comma-separated alpha order")
    ap.add_argument("--model_order", default="", help="Optional comma-separated model order")
    ap.add_argument("--facet_mode", default="obs_then_alpha",
                    choices=["obs_then_alpha", "alpha_then_obs"],
                    help="Column ordering inside each model-specific multipanel.")

    ap.add_argument("--x_min", type=float, default=None, help="Optional lower x bound for data filtering")
    ap.add_argument("--x_max", type=float, default=None, help="Optional upper x bound for data filtering")
    ap.add_argument("--x_mode", default="common", choices=["common", "individual"],
                    help="Use common x-limits across all columns of a model, or per-panel x-limits.")
    ap.add_argument("--x_plot_min", type=float, default=None,
                    help="Lower x limit for plotting only; does not filter the data.")
    ap.add_argument("--x_plot_max", type=float, default=None,
                    help="Upper x limit for plotting only; does not filter the data.")

    ap.add_argument("--ratio_smooth_sigma_pts", type=float, default=1.0)
    ap.add_argument("--flux_smooth_sigma_pts", type=float, default=1.5)
    ap.add_argument("--drift_smooth_sigma_pts", type=float, default=0.0)
    ap.add_argument("--diffusion_smooth_sigma_pts", type=float, default=0.0)

    for key in ["pdf", "ratio", "flux", "source", "drift", "diffusion"]:
        ap.add_argument(f"--{key}_ymin", type=float, default=None)
        ap.add_argument(f"--{key}_ymax", type=float, default=None)

    ap.add_argument("--auto_qlo", type=float, default=0.02)
    ap.add_argument("--auto_qhi", type=float, default=0.98)
    ap.add_argument("--auto_pad_frac", type=float, default=0.05)
    ap.add_argument("--flux_symmetric", action="store_true",
                    help="Force symmetric y-limits around zero for J(x).")
    ap.add_argument("--source_symmetric", action="store_true",
                    help="Force symmetric y-limits around zero for S(x).")

    ap.add_argument("--single_width", type=float, default=5.1)
    ap.add_argument("--single_height", type=float, default=3.1)
    ap.add_argument("--multi_width", type=float, default=8)
    ap.add_argument("--multi_height", type=float, default=5)
    ap.add_argument("--bd_multi_height", type=float, default=3)
    ap.add_argument("--dpi", type=int, default=300)

    ap.add_argument("--line_w_pdf", type=float, default=1.7)
    ap.add_argument("--line_w_ratio", type=float, default=1.5)
    ap.add_argument("--line_w_flux", type=float, default=1.0)
    ap.add_argument("--line_w_source", type=float, default=1.0)
    ap.add_argument("--line_w_bd", type=float, default=1.3)
    ap.add_argument("--refline_w", type=float, default=1.0)
    ap.add_argument("--refline_alpha", type=float, default=0.9)

    ap.add_argument("--also_bd", action="store_true",
                    help="If BD columns are available, produce inferred-function plots too.")
    ap.add_argument("--show_legends", action="store_true",
                    help="Show legends in single-condition plots where appropriate.")
    ap.add_argument("--panel_titles", action="store_true",
                    help="Show per-panel condition titles.")
    ap.add_argument("--figure_title", action="store_true",
                    help="Show figure-level title in multipanels.")

    ap.add_argument("--show_ratio_band", action="store_true",
                    help="Shade a band around ratio=1.")
    ap.add_argument("--ratio_band_halfwidth", type=float, default=0.10,
                    help="Half-width of shaded ratio band around 1.")
    ap.add_argument("--show_zero_band", action="store_true",
                    help="Shade a band around zero in J(x) and S(x).")
    ap.add_argument("--flux_zero_band_halfwidth", type=float, default=0.03,
                    help="Half-width of shaded zero band for J(x).")
    ap.add_argument("--source_zero_band_halfwidth", type=float, default=1.0,
                    help="Half-width of shaded zero band for S(x).")

    ap.add_argument("--only_multipanel", action="store_true",
                    help="Produce only multipanel figures (FP multipanel, and BD multipanel if --also_bd is active).")
    ap.add_argument("--only_fp_multipanel", action="store_true",
                    help="Produce only the FP multipanel. Skip all single plots and all BD plots.")
    ap.add_argument("--only_bd_multipanel", action="store_true",
                    help="Produce only the BD multipanel. Skip all single plots and all FP plots.")
    ap.add_argument("--bd_overlay", action="store_true",
                    help="Produce overlay plots for drift and diffusion with all loaded conditions on the same axes.")
    ap.add_argument("--bd_overlay_by", default="condition", choices=["condition", "obs_class"],
                    help="Legend labeling mode for BD overlay plots.")
    ap.add_argument("--bd_overlay_legend_outside", action="store_true",
                    help="Place BD overlay legend outside the axes on the right.")

    return ap


def setup_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.titlesize": 11,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": True,
        "axes.spines.right": True,
    })


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def parse_ordered_strings(s: str) -> List[str]:
    if not s:
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def parse_ordered_floats(s: str) -> List[float]:
    if not s:
        return []
    out = []
    for x in str(s).split(","):
        x = x.strip()
        if x:
            out.append(float(x))
    return out


def alpha_to_label(alpha: float) -> str:
    if alpha >= 0.01:
        return f"p={alpha:g}"
    return f"p={alpha:.3f}".rstrip("0").rstrip(".")


def alpha_to_filename(alpha: float) -> str:
    s = f"{alpha:.6g}"
    return f"p_{s.replace('.', '_')}"


def model_to_filename(model: str) -> str:
    return str(model).replace("/", "_").replace(" ", "_")


def condition_label(obs_class: str, alpha: float) -> str:
    return f"{obs_class} | {alpha_to_label(alpha)}"


def condition_filename(obs_class: str, alpha: float) -> str:
    obs = str(obs_class).strip().replace(" ", "_")
    return f"{obs}__{alpha_to_filename(alpha)}"


def prettify_obs_class(obs: str) -> str:
    s = str(obs).strip()
    mapping = {
        "ALL_unweighted": "ALL, unweighted",
        "ALL_w_clone": "ALL, clone-weighted",
        "TT_unweighted": "TT, unweighted",
        "TT_w_clone": "TT, clone-weighted",
        "TT": "TT",
        "ALL": "ALL",
        "TF": "TF",
        "FT": "FT",
        "FF": "FF",
    }
    if s in mapping:
        return mapping[s]
    s = s.replace("_w_clone", ", clone-weighted")
    s = s.replace("_unweighted", ", unweighted")
    s = s.replace("_", " ")
    return s


def pretty_condition_title(obs_class: str, alpha: float) -> str:
    return f"{prettify_obs_class(obs_class)} | {alpha_to_label(alpha)}"


def normalize_ylim(ymin: Optional[float], ymax: Optional[float]) -> Optional[Tuple[Optional[float], Optional[float]]]:
    if ymin is None or ymax is None:
        return None
    if ymin is not None and not np.isfinite(ymin):
        ymin = None
    if ymax is not None and not np.isfinite(ymax):
        ymax = None
    if ymin is not None and ymax is not None and ymax <= ymin:
        return None

    return (ymin, ymax)


def load_json(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def gaussian_kernel1d(sigma_pts: float, truncate: float = 4.0) -> np.ndarray:
    sigma = float(sigma_pts)
    if (not np.isfinite(sigma)) or sigma <= 0:
        return np.array([1.0], dtype=float)
    radius = max(1, int(truncate * sigma + 0.5))
    x = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-(x * x) / (2.0 * sigma * sigma))
    k /= np.sum(k)
    return k


def smooth_series(y: np.ndarray, sigma_pts: float) -> np.ndarray:
    y = np.asarray(y, float)
    if y.size == 0:
        return y.copy()
    if (not np.isfinite(sigma_pts)) or sigma_pts <= 0:
        return y.copy()
    k = gaussian_kernel1d(sigma_pts)
    pad = len(k) // 2
    yp = np.pad(y, pad, mode="edge")
    ys = np.convolve(yp, k, mode="same")[pad:-pad]
    return ys


def pick_first_present(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def finite_mask(*arrs: np.ndarray) -> np.ndarray:
    m = None
    for a in arrs:
        aa = np.asarray(a, float)
        mm = np.isfinite(aa)
        m = mm if m is None else (m & mm)
    return m if m is not None else np.array([], dtype=bool)


def apply_x_filter(
    x: np.ndarray,
    arrays: List[np.ndarray],
    x_min: Optional[float],
    x_max: Optional[float],
) -> Tuple[np.ndarray, List[np.ndarray]]:
    x = np.asarray(x, float)
    m = np.isfinite(x)
    if x_min is not None:
        m &= (x >= float(x_min))
    if x_max is not None:
        m &= (x <= float(x_max))
    x2 = x[m]
    out = [np.asarray(a)[m] for a in arrays]
    return x2, out


def robust_limits_from_arrays(
    arrays: List[np.ndarray],
    qlo: float = 0.02,
    qhi: float = 0.98,
    symmetric: bool = False,
    pad_frac: float = 0.05,
) -> Optional[Tuple[float, float]]:
    vals = []
    for a in arrays:
        aa = np.asarray(a, float)
        aa = aa[np.isfinite(aa)]
        if aa.size:
            vals.append(aa)

    if not vals:
        return None

    x = np.concatenate(vals)
    if x.size == 0:
        return None

    lo = float(np.quantile(x, qlo))
    hi = float(np.quantile(x, qhi))

    if symmetric:
        m = max(abs(lo), abs(hi))
        lo, hi = -m, m

    if hi <= lo:
        m = float(np.nanmean(x)) if np.isfinite(np.nanmean(x)) else 0.0
        lo, hi = m - 1.0, m + 1.0

    span = hi - lo
    pad = pad_frac * span if span > 0 else 1.0
    return (lo - pad, hi + pad)


def maybe_set_ylim(ax, lim: Optional[Tuple[float, float]]) -> None:
    if lim is not None:
        ax.set_ylim(lim[0], lim[1])


def maybe_set_xlim(ax, lim: Optional[Tuple[float, float]]) -> None:
    if lim is not None:
        ax.set_xlim(lim[0], lim[1])


def infer_alpha_from_text(text: str) -> Optional[float]:
    patterns = [
        r'(?:^|[/_ \-])p[_\-]?(\d+(?:[._]\d+)?)',
        r'(?:^|[/_ \-])alpha[_\-]?(\d+(?:[._]\d+)?)',
        r'(?:^|[/_ \-])a[_\-]?(\d+(?:[._]\d+)?)',
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            raw = m.group(1).replace("_", ".")
            try:
                val = float(raw)
                if val > 1:
                    digits = m.group(1).replace("_", "").replace(".", "")
                    if digits.isdigit():
                        val = float("0." + digits)
                return val
            except Exception:
                pass

    m = re.search(r'(?:^|[/_ \-])p[_\-]?0*([1-9]\d*)$', text, flags=re.IGNORECASE)
    if m:
        digits = m.group(1)
        try:
            return float("0." + digits)
        except Exception:
            pass
    return None


def infer_obs_class_from_text(text: str) -> Optional[str]:
    for obs in ["ALL_unweighted", "ALL_w_clone", "TT_unweighted", "TT_w_clone"]:
        if obs in text:
            return obs
    for obs in KNOWN_OBS:
        if re.search(rf'(^|[/_ \-]){obs}($|[/_ \-])', text, flags=re.IGNORECASE):
            return obs
    return None


def infer_model_from_text(text: str) -> Optional[str]:
    for model in KNOWN_MODELS:
        if model in text:
            return model
    return None


def paired_summary_from_grid(grid_csv: Path) -> Optional[Path]:
    p = grid_csv.parent / grid_csv.name.replace(".grid.csv", ".summary.json")
    if p.exists():
        return p
    p2 = grid_csv.with_suffix("").with_suffix(".summary.json")
    if p2.exists():
        return p2
    return None


def discover_pairs_from_root(root: Path) -> pd.DataFrame:
    rows = []
    grid_files = list(root.rglob("*.grid.csv"))

    for grid_csv in grid_files:
        summary_json = paired_summary_from_grid(grid_csv)
        if summary_json is None:
            continue

        try:
            meta = load_json(summary_json)
        except Exception:
            meta = {}

        text = str(grid_csv.parent) + " " + grid_csv.name

        alpha = None
        if meta:
            a = safe_float(meta.get("alpha", np.nan))
            if np.isfinite(a):
                alpha = a
        if alpha is None:
            alpha = infer_alpha_from_text(text)

        obs_class = None
        if meta:
            obs_class = meta.get("obs_class_requested", None)
            if obs_class is None:
                obs_class = meta.get("obs_class", None)
            if obs_class is not None:
                obs_class = str(obs_class)
        if obs_class is None:
            obs_class = infer_obs_class_from_text(text)

        model = None
        if meta:
            model = meta.get("model", None)
            if model is not None:
                model = str(model)
        if model is None:
            model = infer_model_from_text(text)
        if model is None:
            model = "unknown_model"

        if alpha is None:
            continue
        if obs_class is None:
            obs_class = "UNKNOWN"

        rows.append({
            "obs_class": str(obs_class),
            "alpha": float(alpha),
            "model": str(model),
            "grid_csv": str(grid_csv.resolve()),
            "summary_json": str(summary_json.resolve()),
        })

    if not rows:
        raise ValueError(f"No compatible *.grid.csv / *.summary.json pairs found under: {root}")

    df = pd.DataFrame(rows).drop_duplicates(subset=["obs_class", "alpha", "model", "grid_csv"])
    return df.sort_values(["model", "obs_class", "alpha"]).reset_index(drop=True)


def resolve_columns(g: pd.DataFrame) -> Dict[str, Optional[str]]:
    return {
        "x": pick_first_present(g, ["x", "x_mid", "log_freq", "logf"]),
        "p_emp": pick_first_present(g, ["p_emp", "p_empirical", "pdf_emp", "p_observed"]),
        "p_closed": pick_first_present(g, ["p_closed", "p_model", "pdf_model", "p_stationary"]),
        "ratio": pick_first_present(g, ["ratio_pemp_over_pclosed", "ratio_emp_over_model", "ratio", "density_ratio"]),
        "J": pick_first_present(g, ["J_smooth", "J_emp", "J_model", "J", "flux", "probability_flux"]),
        "S": pick_first_present(g, ["S_smooth", "S_emp", "S_model", "S", "source", "source_term"]),
        "b": pick_first_present(g, ["b", "b_model", "drift", "b_hinge_plateau_smooth", "b_OU_linear"]),
        "D": pick_first_present(g, ["D_final", "D", "diffusion", "D_smooth", "D_raw_interp"]),
    }


def load_dataset(model: str, obs_class: str, alpha: float, grid_csv: Path, summary_json: Path, args) -> Dict:
    g = pd.read_csv(grid_csv)
    meta = load_json(summary_json)
    col = resolve_columns(g)

    if col["x"] is None:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] Missing x column in grid: {grid_csv}")
    if col["p_emp"] is None:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] Missing empirical density column in grid: {grid_csv}")
    if col["p_closed"] is None:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] Missing stationary/model density column in grid: {grid_csv}")
    if col["J"] is None:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] Missing flux column in grid: {grid_csv}")
    if col["S"] is None:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] Missing source column in grid: {grid_csv}")

    x = g[col["x"]].to_numpy(float)
    p_emp = g[col["p_emp"]].to_numpy(float)
    p_closed = g[col["p_closed"]].to_numpy(float)

    if col["ratio"] is not None:
        ratio = g[col["ratio"]].to_numpy(float)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = p_emp / p_closed

    J = g[col["J"]].to_numpy(float)
    S = g[col["S"]].to_numpy(float)

    b = g[col["b"]].to_numpy(float) if col["b"] is not None else None
    D = g[col["D"]].to_numpy(float) if col["D"] is not None else None

    if b is not None and D is not None:
        m = finite_mask(x, p_emp, p_closed, ratio, J, S, b, D)
        x = x[m]; p_emp = p_emp[m]; p_closed = p_closed[m]
        ratio = ratio[m]; J = J[m]; S = S[m]; b = b[m]; D = D[m]
        x, [p_emp, p_closed, ratio, J, S, b, D] = apply_x_filter(
            x, [p_emp, p_closed, ratio, J, S, b, D], args.x_min, args.x_max
        )
    else:
        m = finite_mask(x, p_emp, p_closed, ratio, J, S)
        x = x[m]; p_emp = p_emp[m]; p_closed = p_closed[m]
        ratio = ratio[m]; J = J[m]; S = S[m]
        x, [p_emp, p_closed, ratio, J, S] = apply_x_filter(
            x, [p_emp, p_closed, ratio, J, S], args.x_min, args.x_max
        )

    if x.size == 0:
        raise ValueError(f"[{model}, {obs_class}, alpha={alpha}] No rows remain after x filtering.")

    ratio_s = smooth_series(ratio, args.ratio_smooth_sigma_pts)
    J_s = smooth_series(J, args.flux_smooth_sigma_pts)
    S_s = smooth_series(S, args.flux_smooth_sigma_pts)

    b_s = smooth_series(b, args.drift_smooth_sigma_pts) if b is not None else None
    D_s = smooth_series(D, args.diffusion_smooth_sigma_pts) if D is not None else None

    if model == "unknown_model":
        model = str(meta.get("model", "unknown_model"))

    return {
        "model": str(model),
        "model_file": model_to_filename(str(model)),
        "obs_class": str(obs_class),
        "pretty_obs_class": prettify_obs_class(str(obs_class)),
        "alpha": float(alpha),
        "alpha_label": alpha_to_label(float(alpha)),
        "alpha_file": alpha_to_filename(float(alpha)),
        "condition_label": condition_label(str(obs_class), float(alpha)),
        "pretty_condition_title": pretty_condition_title(str(obs_class), float(alpha)),
        "condition_file": condition_filename(str(obs_class), float(alpha)),
        "grid_csv": str(grid_csv),
        "summary_json": str(summary_json),
        "x": x,
        "xlim_data": (float(np.min(x)), float(np.max(x))),
        "p_emp": p_emp,
        "p_closed": p_closed,
        "ratio": ratio_s,
        "J": J_s,
        "S": S_s,
        "b": b_s,
        "D": D_s,
        "JS": safe_float(meta.get("JS", np.nan)),
        "KS": safe_float(meta.get("KS", np.nan)),
        "tag": meta.get("tag", ""),
    }


def ordered_conditions(datasets: List[Dict], obs_order: List[str], alpha_order: List[float], facet_mode: str) -> List[Dict]:
    obs_seen = sorted({d["obs_class"] for d in datasets})
    alpha_seen = sorted({float(d["alpha"]) for d in datasets})

    obs_final = obs_order[:] if obs_order else obs_seen
    for o in obs_seen:
        if o not in obs_final:
            obs_final.append(o)

    alpha_final = alpha_order[:] if alpha_order else alpha_seen
    for a in alpha_seen:
        if a not in alpha_final:
            alpha_final.append(a)

    lookup = {(d["obs_class"], float(d["alpha"])): d for d in datasets}
    out = []

    if facet_mode == "obs_then_alpha":
        for obs in obs_final:
            for a in alpha_final:
                if (obs, a) in lookup:
                    out.append(lookup[(obs, a)])
    else:
        for a in alpha_final:
            for obs in obs_final:
                if (obs, a) in lookup:
                    out.append(lookup[(obs, a)])

    return out


def resolve_global_ylims(
    datasets: List[Dict],
    pdf_ylim_user: Optional[Tuple[float, float]],
    ratio_ylim_user: Optional[Tuple[float, float]],
    flux_ylim_user: Optional[Tuple[float, float]],
    source_ylim_user: Optional[Tuple[float, float]],
    drift_ylim_user: Optional[Tuple[float, float]],
    diffusion_ylim_user: Optional[Tuple[float, float]],
    auto_qlo: float,
    auto_qhi: float,
    auto_pad_frac: float,
    flux_symmetric: bool,
    source_symmetric: bool,
) -> Dict[str, Optional[Tuple[float, float]]]:
    pdf_ylim = pdf_ylim_user
    if pdf_ylim is None:
        pdf_ylim = robust_limits_from_arrays(
            [d["p_emp"] for d in datasets] + [d["p_closed"] for d in datasets],
            qlo=0.0, qhi=1.0, symmetric=False, pad_frac=auto_pad_frac
        )

    ratio_ylim = ratio_ylim_user
    if ratio_ylim is None:
        ratio_ylim = robust_limits_from_arrays(
            [d["ratio"] for d in datasets],
            qlo=auto_qlo, qhi=auto_qhi, symmetric=False, pad_frac=auto_pad_frac
        )

    flux_ylim = flux_ylim_user
    if flux_ylim is None:
        flux_ylim = robust_limits_from_arrays(
            [d["J"] for d in datasets],
            qlo=auto_qlo, qhi=auto_qhi, symmetric=flux_symmetric, pad_frac=auto_pad_frac
        )

    source_ylim = source_ylim_user
    if source_ylim is None:
        source_ylim = robust_limits_from_arrays(
            [d["S"] for d in datasets],
            qlo=auto_qlo, qhi=auto_qhi, symmetric=source_symmetric, pad_frac=auto_pad_frac
        )

    drift_arrays = [d["b"] for d in datasets if d["b"] is not None]
    drift_ylim = drift_ylim_user
    if drift_ylim is None and drift_arrays:
        drift_ylim = robust_limits_from_arrays(
            drift_arrays, qlo=auto_qlo, qhi=auto_qhi, symmetric=False, pad_frac=auto_pad_frac
        )

    diffusion_arrays = [d["D"] for d in datasets if d["D"] is not None]
    diffusion_ylim = diffusion_ylim_user
    if diffusion_ylim is None and diffusion_arrays:
        diffusion_ylim = robust_limits_from_arrays(
            diffusion_arrays, qlo=auto_qlo, qhi=auto_qhi, symmetric=False, pad_frac=auto_pad_frac
        )

    return {
        "pdf": pdf_ylim,
        "ratio": ratio_ylim,
        "flux": flux_ylim,
        "source": source_ylim,
        "drift": drift_ylim,
        "diffusion": diffusion_ylim,
    }


def resolve_base_xlim(datasets: List[Dict], x_mode: str) -> Optional[Tuple[float, float]]:
    if x_mode == "individual":
        return None
    vals = []
    for d in datasets:
        x = np.asarray(d["x"], float)
        x = x[np.isfinite(x)]
        if x.size:
            vals.append((float(np.min(x)), float(np.max(x))))
    if not vals:
        return None
    return (min(v[0] for v in vals), max(v[1] for v in vals))


def resolve_plot_xlim(base_xlim: Optional[Tuple[float, float]], args) -> Optional[Tuple[float, float]]:
    if base_xlim is None and args.x_plot_min is None and args.x_plot_max is None:
        return None

    if base_xlim is None:
        xmin = args.x_plot_min
        xmax = args.x_plot_max
        if xmin is None or xmax is None:
            return None
        return (float(xmin), float(xmax))

    xmin = base_xlim[0] if args.x_plot_min is None else float(args.x_plot_min)
    xmax = base_xlim[1] if args.x_plot_max is None else float(args.x_plot_max)
    if xmax <= xmin:
        raise ValueError(f"Invalid plotting x-range: ({xmin}, {xmax})")
    return (xmin, xmax)


def save_summary_table(datasets: List[Dict], out_dir: Path, model: str) -> None:
    rows = []
    for d in datasets:
        x = d["x"]
        J = d["J"]
        S = d["S"]
        ratio = d["ratio"]

        rows.append({
            "model": d["model"],
            "obs_class": d["obs_class"],
            "alpha": d["alpha"],
            "alpha_label": d["alpha_label"],
            "condition_label": d["condition_label"],
            "JS": d["JS"],
            "KS": d["KS"],
            "int_abs_J": float(np.trapezoid(np.abs(J), x)) if len(x) > 1 else np.nan,
            "int_abs_S": float(np.trapezoid(np.abs(S), x)) if len(x) > 1 else np.nan,
            "mean_abs_J": float(np.mean(np.abs(J))) if len(x) else np.nan,
            "mean_abs_S": float(np.mean(np.abs(S))) if len(x) else np.nan,
            "mean_ratio": float(np.mean(ratio)) if len(x) else np.nan,
            "median_ratio": float(np.median(ratio)) if len(x) else np.nan,
            "n_grid": int(len(x)),
            "grid_csv": d["grid_csv"],
            "summary_json": d["summary_json"],
            "tag": d["tag"],
        })

    df = pd.DataFrame(rows).sort_values(["obs_class", "alpha"]).reset_index(drop=True)
    suffix = model_to_filename(model)
    df.to_csv(out_dir / f"fp_sensitivity_summary__{suffix}.csv", index=False)

    pretty = df.copy()
    for c in ["alpha", "JS", "KS", "int_abs_J", "int_abs_S", "mean_abs_J", "mean_abs_S", "mean_ratio", "median_ratio"]:
        pretty[c] = pretty[c].map(lambda v: f"{v:.6g}" if pd.notna(v) else "")
    pretty.to_csv(out_dir / f"fp_sensitivity_summary_pretty__{suffix}.csv", index=False)


def write_report(datasets: List[Dict], out_dir: Path, args, ylims: Dict[str, Optional[Tuple[float, float]]], model: str, base_xlim: Optional[Tuple[float, float]], plot_xlim: Optional[Tuple[float, float]]) -> None:
    suffix = model_to_filename(model)
    lines = [
        f"FP sensitivity diagnostics report | model={model}",
        "===============================================",
        "",
        f"Input manifest_csv: {args.manifest_csv}",
        f"Input input_root: {args.input_root}",
        f"Output dir: {args.out_dir}",
        f"Format: {args.fmt}",
        "",
        "Global settings",
        f"  x_min = {args.x_min}",
        f"  x_max = {args.x_max}",
        f"  x_mode = {args.x_mode}",
        f"  base xlim = {base_xlim}",
        f"  x_plot_min = {args.x_plot_min}",
        f"  x_plot_max = {args.x_plot_max}",
        f"  final plot xlim = {plot_xlim}",
        f"  obs_order = {args.obs_order}",
        f"  alpha_order = {args.alpha_order}",
        f"  model_order = {args.model_order}",
        f"  facet_mode = {args.facet_mode}",
        f"  ratio_smooth_sigma_pts = {args.ratio_smooth_sigma_pts}",
        f"  flux_smooth_sigma_pts = {args.flux_smooth_sigma_pts}",
        f"  drift_smooth_sigma_pts = {args.drift_smooth_sigma_pts}",
        f"  diffusion_smooth_sigma_pts = {args.diffusion_smooth_sigma_pts}",
        f"  only_multipanel = {args.only_multipanel}",
        f"  only_fp_multipanel = {args.only_fp_multipanel}",
        f"  only_bd_multipanel = {args.only_bd_multipanel}",
        f"  bd_overlay = {args.bd_overlay}",
        f"  bd_overlay_by = {args.bd_overlay_by}",
        f"  bd_overlay_legend_outside = {args.bd_overlay_legend_outside}",
        "",
        "Resolved common y-limits",
        f"  pdf       = {ylims['pdf']}",
        f"  ratio     = {ylims['ratio']}",
        f"  flux      = {ylims['flux']}",
        f"  source    = {ylims['source']}",
        f"  drift     = {ylims['drift']}",
        f"  diffusion = {ylims['diffusion']}",
        "",
        "Datasets",
    ]

    for d in datasets:
        lines.append(
            f"  obs_class={d['obs_class']}  alpha={d['alpha']:.6g}  tag={d['tag']}  "
            f"n_grid={len(d['x'])}  JS={d['JS']:.6g}  KS={d['KS']:.6g}"
        )
        lines.append(f"    grid:    {d['grid_csv']}")
        lines.append(f"    summary: {d['summary_json']}")

    (out_dir / f"fp_sensitivity_report__{suffix}.txt").write_text("\n".join(lines), encoding="utf-8")


def add_ratio_band(ax, args):
    if args.show_ratio_band:
        ax.axhspan(1.0 - args.ratio_band_halfwidth, 1.0 + args.ratio_band_halfwidth,
                   color="0.9", zorder=0)


def add_zero_band(ax, halfwidth: float, args):
    if args.show_zero_band:
        ax.axhspan(-halfwidth, halfwidth, color="0.92", zorder=0)


def plot_metrics_by_obsclass(datasets: List[Dict], out_dir: Path, fmt: str, dpi: int, model: str) -> None:
    obs_classes = sorted({d["obs_class"] for d in datasets})
    suffix = model_to_filename(model)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    metrics = [("JS", axes[0]), ("KS", axes[1])]

    for metric, ax in metrics:
        for obs in obs_classes:
            ds = sorted([d for d in datasets if d["obs_class"] == obs], key=lambda z: z["alpha"])
            alphas = [d["alpha"] for d in ds]
            vals = [d[metric] for d in ds]
            ax.plot(alphas, vals, "-o", label=prettify_obs_class(obs), lw=1.6, ms=4)
        ax.set_xscale("log")
        ax.set_xlabel("alpha")
        ax.set_ylabel(metric)
        ax.set_title(metric)
        ax.legend(frameon=False)

    fig.savefig(out_dir / f"fig_fp_sensitivity_metrics_by_obsclass__{suffix}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_pdf(d: Dict, out_dir: Path, fmt: str, dpi: int, w: float, h: float,
                    pdf_ylim, plot_xlim, panel_titles: bool, show_legends: bool, args) -> None:
    fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
    ax.plot(d["x"], d["p_emp"], label="empirical", lw=args.line_w_pdf)
    ax.plot(d["x"], d["p_closed"], label="stationary", lw=args.line_w_pdf)
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("density")
    if panel_titles:
        ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
    maybe_set_ylim(ax, pdf_ylim)
    maybe_set_xlim(ax, plot_xlim)
    if show_legends:
        ax.legend(frameon=False)
    fig.savefig(out_dir / f"fig_fp_pdf__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_ratio(d: Dict, out_dir: Path, fmt: str, dpi: int, w: float, h: float,
                      ratio_ylim, plot_xlim, panel_titles: bool, args) -> None:
    fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
    add_ratio_band(ax, args)
    ax.plot(d["x"], d["ratio"], lw=args.line_w_ratio)
    ax.axhline(1.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("p_emp / p_closed")
    if panel_titles:
        ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
    maybe_set_ylim(ax, ratio_ylim)
    maybe_set_xlim(ax, plot_xlim)
    fig.savefig(out_dir / f"fig_fp_ratio__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_flux(d: Dict, out_dir: Path, fmt: str, dpi: int, w: float, h: float,
                     flux_ylim, plot_xlim, panel_titles: bool, args) -> None:
    fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
    add_zero_band(ax, args.flux_zero_band_halfwidth, args)
    ax.plot(d["x"], d["J"], lw=args.line_w_flux)
    ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("J(x)")
    if panel_titles:
        ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
    maybe_set_ylim(ax, flux_ylim)
    maybe_set_xlim(ax, plot_xlim)
    fig.savefig(out_dir / f"fig_fp_flux__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_source(d: Dict, out_dir: Path, fmt: str, dpi: int, w: float, h: float,
                       source_ylim, plot_xlim, panel_titles: bool, args) -> None:
    fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
    add_zero_band(ax, args.source_zero_band_halfwidth, args)
    ax.plot(d["x"], d["S"], lw=args.line_w_source)
    ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("S(x)")
    if panel_titles:
        ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
    maybe_set_ylim(ax, source_ylim)
    maybe_set_xlim(ax, plot_xlim)
    fig.savefig(out_dir / f"fig_fp_source__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_bd(d: Dict, out_dir: Path, fmt: str, dpi: int, w: float, h: float,
                   drift_ylim, diffusion_ylim, plot_xlim, panel_titles: bool, args) -> None:
    if d["b"] is not None:
        fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
        ax.plot(d["x"], d["b"], lw=args.line_w_bd)
        ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
        ax.set_xlabel("x = ln f")
        ax.set_ylabel("b(x)")
        if panel_titles:
            ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
        maybe_set_ylim(ax, drift_ylim)
        maybe_set_xlim(ax, plot_xlim)
        fig.savefig(out_dir / f"fig_bd_drift__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
        plt.close(fig)

    if d["D"] is not None:
        fig, ax = plt.subplots(figsize=(w, h), constrained_layout=True)
        ax.plot(d["x"], d["D"], lw=args.line_w_bd)
        ax.set_xlabel("x = ln f")
        ax.set_ylabel("D(x)")
        if panel_titles:
            ax.set_title(f"{d['model']} | {d['pretty_condition_title']}")
        maybe_set_ylim(ax, diffusion_ylim)
        maybe_set_xlim(ax, plot_xlim)
        fig.savefig(out_dir / f"fig_bd_diffusion__{d['model_file']}__{d['condition_file']}.{fmt}", dpi=dpi)
        plt.close(fig)


def plot_fp_multipanel(
    conditions: List[Dict],
    out_dir: Path,
    fmt: str,
    dpi: int,
    fig_w: float,
    fig_h: float,
    pdf_ylim,
    ratio_ylim,
    flux_ylim,
    source_ylim,
    plot_xlim,
    args,
    model: str,
) -> None:
    ncols = len(conditions)
    suffix = model_to_filename(model)
    fig, axes = plt.subplots(4, ncols, figsize=(fig_w, fig_h), sharex=False, constrained_layout=True)
    if ncols == 1:
        axes = np.array(axes).reshape(4, 1)

    row_labels = ["Density", "Density ratio", "Probability current", "Source term"]

    for col, d in enumerate(conditions):
        ax = axes[0, col]
        ax.plot(d["x"], d["p_emp"], label="empirical", lw=args.line_w_pdf)
        ax.plot(d["x"], d["p_closed"], label="stationary", lw=args.line_w_pdf)
        if args.panel_titles:
            ax.set_title(d["pretty_condition_title"])
        if col == 0:
            ax.set_ylabel(row_labels[0])
        maybe_set_ylim(ax, pdf_ylim)
        maybe_set_xlim(ax, plot_xlim)
        if col == ncols - 1:
            ax.legend(frameon=False, loc="upper right")

        ax = axes[1, col]
        add_ratio_band(ax, args)
        ax.plot(d["x"], d["ratio"], lw=args.line_w_ratio)
        ax.axhline(1.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
        if col == 0:
            ax.set_ylabel(row_labels[1])
        maybe_set_ylim(ax, ratio_ylim)
        maybe_set_xlim(ax, plot_xlim)

        ax = axes[2, col]
        add_zero_band(ax, args.flux_zero_band_halfwidth, args)
        ax.plot(d["x"], d["J"], lw=args.line_w_flux)
        ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
        if col == 0:
            ax.set_ylabel(row_labels[2])
        maybe_set_ylim(ax, flux_ylim)
        maybe_set_xlim(ax, plot_xlim)

        ax = axes[3, col]
        add_zero_band(ax, args.source_zero_band_halfwidth, args)
        ax.plot(d["x"], d["S"], lw=args.line_w_source)
        ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
        if col == 0:
            ax.set_ylabel(row_labels[3])
        ax.set_xlabel("x = ln f")
        maybe_set_ylim(ax, source_ylim)
        maybe_set_xlim(ax, plot_xlim)

    if args.figure_title:
        fig.suptitle(f"{model}: FP sensitivity diagnostics", y=1.02)

    fig.savefig(out_dir / f"fig_fp_sensitivity_multipanel_4x{ncols}__{suffix}.{fmt}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_bd_multipanel(
    conditions: List[Dict],
    out_dir: Path,
    fmt: str,
    dpi: int,
    fig_w: float,
    fig_h: float,
    drift_ylim,
    diffusion_ylim,
    plot_xlim,
    args,
    model: str,
) -> None:
    valid = [d for d in conditions if d["b"] is not None and d["D"] is not None]
    if not valid:
        return

    suffix = model_to_filename(model)
    ncols = len(valid)
    fig, axes = plt.subplots(2, ncols, figsize=(fig_w, fig_h), sharex=False, constrained_layout=True)
    if ncols == 1:
        axes = np.array(axes).reshape(2, 1)

    row_labels = ["Drift", "Diffusion"]

    for col, d in enumerate(valid):
        ax = axes[0, col]
        ax.plot(d["x"], d["b"], lw=args.line_w_bd)
        ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
        if args.panel_titles:
            ax.set_title(d["pretty_condition_title"])
        if col == 0:
            ax.set_ylabel(row_labels[0])
        maybe_set_ylim(ax, drift_ylim)
        maybe_set_xlim(ax, plot_xlim)

        ax = axes[1, col]
        ax.plot(d["x"], d["D"], lw=args.line_w_bd)
        if col == 0:
            ax.set_ylabel(row_labels[1])
        ax.set_xlabel("x = ln f")
        maybe_set_ylim(ax, diffusion_ylim)
        maybe_set_xlim(ax, plot_xlim)

    if args.figure_title:
        fig.suptitle(f"{model}: inferred functions", y=1.02)

    fig.savefig(out_dir / f"fig_bd_sensitivity_multipanel_2x{ncols}__{suffix}.{fmt}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)



def plot_bd_overlay(
    conditions: List[Dict],
    out_dir: Path,
    fmt: str,
    dpi: int,
    fig_w: float,
    fig_h: float,
    drift_ylim,
    diffusion_ylim,
    plot_xlim,
    args,
    model: str,
) -> None:
    """
    Overlay drift and diffusion curves for all loaded conditions of one model.
    """
    valid = [d for d in conditions if d["b"] is not None and d["D"] is not None]
    if not valid:
        return

    suffix = model_to_filename(model)

    def make_label(d: Dict) -> str:
        if args.bd_overlay_by == "obs_class":
            return d["pretty_obs_class"] + f" | {d['alpha_label']}"
        return d["pretty_condition_title"]

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    for d in valid:
        ax.plot(d["x"], d["b"], lw=args.line_w_bd, label=make_label(d))
    ax.axhline(0.0, ls="--", lw=args.refline_w, alpha=args.refline_alpha)
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("Drift")
    if args.figure_title:
        ax.set_title(f"{model}: drift overlay")
    maybe_set_ylim(ax, drift_ylim)
    maybe_set_xlim(ax, plot_xlim)
    if args.bd_overlay_legend_outside:
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)
        fig.savefig(out_dir / f"fig_bd_overlay_drift__{suffix}.{fmt}", dpi=dpi, bbox_inches="tight")
    else:
        ax.legend(frameon=False)
        fig.savefig(out_dir / f"fig_bd_overlay_drift__{suffix}.{fmt}", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    for d in valid:
        ax.plot(d["x"], d["D"], lw=args.line_w_bd, label=make_label(d))
    ax.set_xlabel("x = ln f")
    ax.set_ylabel("Diffusion")
    if args.figure_title:
        ax.set_title(f"{model}: diffusion overlay")
    maybe_set_ylim(ax, diffusion_ylim)
    maybe_set_xlim(ax, plot_xlim)
    if args.bd_overlay_legend_outside:
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)
        fig.savefig(out_dir / f"fig_bd_overlay_diffusion__{suffix}.{fmt}", dpi=dpi, bbox_inches="tight")
    else:
        ax.legend(frameon=False)
        fig.savefig(out_dir / f"fig_bd_overlay_diffusion__{suffix}.{fmt}", dpi=dpi)
    plt.close(fig)

def main() -> None:
    args = build_argparser().parse_args()
    setup_matplotlib()

    selected = [args.only_multipanel, args.only_fp_multipanel, args.only_bd_multipanel]
    if sum(bool(x) for x in selected) > 1:
        raise ValueError("Use at most one of: --only_multipanel, --only_fp_multipanel, --only_bd_multipanel")

    if not args.manifest_csv and not args.input_root:
        raise ValueError("Provide either --manifest_csv or --input_root.")
    if args.manifest_csv and args.input_root:
        raise ValueError("Provide only one of --manifest_csv or --input_root.")

    out_dir = ensure_dir(Path(args.out_dir).resolve())

    pdf_ylim_user = normalize_ylim(args.pdf_ymin, args.pdf_ymax)
    ratio_ylim_user = normalize_ylim(args.ratio_ymin, args.ratio_ymax)
    flux_ylim_user = normalize_ylim(args.flux_ymin, args.flux_ymax)
    source_ylim_user = normalize_ylim(args.source_ymin, args.source_ymax)
    drift_ylim_user = normalize_ylim(args.drift_ymin, args.drift_ymax)
    diffusion_ylim_user = normalize_ylim(args.diffusion_ymin, args.diffusion_ymax)

    if args.input_root:
        manifest = discover_pairs_from_root(Path(args.input_root).expanduser().resolve())
    else:
        manifest = pd.read_csv(Path(args.manifest_csv).expanduser().resolve())
        required = {"obs_class", "alpha", "grid_csv", "summary_json"}
        missing = required - set(manifest.columns)
        if missing:
            raise ValueError(f"manifest_csv missing required columns: {sorted(missing)}")
        if "model" not in manifest.columns:
            manifest["model"] = "unknown_model"

    manifest = manifest.copy()
    manifest["alpha"] = pd.to_numeric(manifest["alpha"], errors="coerce")
    if manifest["alpha"].isna().any():
        raise ValueError("Column 'alpha' contains non-numeric values.")
    manifest["obs_class"] = manifest["obs_class"].astype(str)
    manifest["model"] = manifest["model"].astype(str)

    datasets: List[Dict] = []
    for _, row in manifest.iterrows():
        model = str(row["model"])
        obs_class = str(row["obs_class"])
        alpha = float(row["alpha"])
        grid_csv = Path(str(row["grid_csv"])).expanduser().resolve()
        summary_json = Path(str(row["summary_json"])).expanduser().resolve()

        if not grid_csv.exists():
            raise FileNotFoundError(f"Grid CSV not found: {grid_csv}")
        if not summary_json.exists():
            raise FileNotFoundError(f"Summary JSON not found: {summary_json}")

        datasets.append(load_dataset(model, obs_class, alpha, grid_csv, summary_json, args))

    model_order = parse_ordered_strings(args.model_order)
    models_seen = [d["model"] for d in datasets]
    unique_models = []
    for m in model_order:
        if m in models_seen and m not in unique_models:
            unique_models.append(m)
    for m in models_seen:
        if m not in unique_models:
            unique_models.append(m)

    obs_order = parse_ordered_strings(args.obs_order)
    alpha_order = parse_ordered_floats(args.alpha_order)

    for model in unique_models:
        ds_model = [d for d in datasets if d["model"] == model]
        conditions = ordered_conditions(ds_model, obs_order, alpha_order, args.facet_mode)
        base_xlim = resolve_base_xlim(conditions, args.x_mode)
        plot_xlim = resolve_plot_xlim(base_xlim, args)

        ylims = resolve_global_ylims(
            datasets=conditions,
            pdf_ylim_user=pdf_ylim_user,
            ratio_ylim_user=ratio_ylim_user,
            flux_ylim_user=flux_ylim_user,
            source_ylim_user=source_ylim_user,
            drift_ylim_user=drift_ylim_user,
            diffusion_ylim_user=diffusion_ylim_user,
            auto_qlo=args.auto_qlo,
            auto_qhi=args.auto_qhi,
            auto_pad_frac=args.auto_pad_frac,
            flux_symmetric=args.flux_symmetric,
            source_symmetric=args.source_symmetric,
        )

        save_summary_table(conditions, out_dir, model)
        write_report(conditions, out_dir, args, ylims, model, base_xlim, plot_xlim)
        plot_metrics_by_obsclass(conditions, out_dir, args.fmt, args.dpi, model)

        make_single = not (args.only_multipanel or args.only_fp_multipanel or args.only_bd_multipanel)
        make_fp_multi = not args.only_bd_multipanel
        make_bd_multi = args.also_bd and not args.only_fp_multipanel

        if make_single:
            for d in conditions:
                plot_single_pdf(
                    d, out_dir, args.fmt, args.dpi,
                    args.single_width, args.single_height,
                    ylims["pdf"], plot_xlim, args.panel_titles, args.show_legends, args
                )
                plot_single_ratio(
                    d, out_dir, args.fmt, args.dpi,
                    args.single_width, args.single_height,
                    ylims["ratio"], plot_xlim, args.panel_titles, args
                )
                plot_single_flux(
                    d, out_dir, args.fmt, args.dpi,
                    args.single_width, args.single_height,
                    ylims["flux"], plot_xlim, args.panel_titles, args
                )
                plot_single_source(
                    d, out_dir, args.fmt, args.dpi,
                    args.single_width, args.single_height,
                    ylims["source"], plot_xlim, args.panel_titles, args
                )
                if args.also_bd:
                    plot_single_bd(
                        d, out_dir, args.fmt, args.dpi,
                        args.single_width, args.single_height,
                        ylims["drift"], ylims["diffusion"], plot_xlim, args.panel_titles, args
                    )

        if make_fp_multi:
            plot_fp_multipanel(
                conditions=conditions,
                out_dir=out_dir,
                fmt=args.fmt,
                dpi=args.dpi,
                fig_w=args.multi_width,
                fig_h=args.multi_height,
                pdf_ylim=ylims["pdf"],
                ratio_ylim=ylims["ratio"],
                flux_ylim=ylims["flux"],
                source_ylim=ylims["source"],
                plot_xlim=plot_xlim,
                args=args,
                model=model,
            )

        if make_bd_multi:
            plot_bd_multipanel(
                conditions=conditions,
                out_dir=out_dir,
                fmt=args.fmt,
                dpi=args.dpi,
                fig_w=args.multi_width,
                fig_h=args.bd_multi_height,
                drift_ylim=ylims["drift"],
                diffusion_ylim=ylims["diffusion"],
                plot_xlim=plot_xlim,
                args=args,
                model=model,
            )

        if args.bd_overlay and args.also_bd:
            plot_bd_overlay(
                conditions=conditions,
                out_dir=out_dir,
                fmt=args.fmt,
                dpi=args.dpi,
                fig_w=args.multi_width * 0.55,
                fig_h=args.bd_multi_height * 0.9,
                drift_ylim=ylims["drift"],
                diffusion_ylim=ylims["diffusion"],
                plot_xlim=plot_xlim,
                args=args,
                model=model,
            )

    print(f"[DONE] Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()


# Overlay example:
# python fp_sensitivity_diagnostics_v6.py \
#   --manifest_csv ./fp_sensitivity_manifest.csv \
#   --out_dir ./fp_sensitivity_overlay \
#   --model_order hinge_plateau_smooth \
#   --also_bd \
#   --bd_overlay \
#   --bd_overlay_legend_outside \
#   --x_plot_min -11.0 \
#   --x_plot_max -6.5
