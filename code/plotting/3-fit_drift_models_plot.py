#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-fit_drift_models_plot.py

Plot-only script producing SINGLE figures (no multipanel).

Reads outputs from the analysis script:
  - drift_function_grid.csv   (must contain column 'x' and one or more 'b_*' columns)
  - binned_medians.csv        (must contain series='dx' and optionally series='resid_*')
Optional:
  - bootstrap_reps.csv        (if present, will overlay bootstrap curves when possible)

Outputs (in --fig_dir)
----------------------
- fig_drift_compare.png (+ optional PDF)
- fig_resid_median_by_bin.png (+ optional PDF)
- fig_resid_compare.png (+ optional PDF; backward-compatible alias)
- fig_bootstrap_overlay.png (+ optional PDF; only if bootstrap_reps.csv exists)

Main fix in this version
------------------------
Colors are now fixed consistently across all figures:
- OU_linear
- tanh_saturating
- rational_saturating
- hinge_plateau_smooth

Each model always keeps the same color in:
1) drift compare
2) residual median by bin
3) bootstrap overlay

Example
-------
python 5-fit_drift_models_plot.py \
  --data_dir ./results/step6_drift_dt1_TT_weighted \
  --fig_dir  ./figures/step6_drift_dt1_TT_weighted \
  --dpi 300 --fig_w 4.8 --fig_h 3.2 \
  --font_family Arial --font_size 10 \
  --legend_outside --legend_loc "upper left"
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Callable, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------
# IO / utils
# ---------------------------------------------------------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def _as_float_col(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        raise SystemExit(f"Missing required column '{col}' in input.")
    return pd.to_numeric(df[col], errors="coerce").to_numpy(float)


def _to_float_array(s: pd.Series | np.ndarray) -> np.ndarray:
    return pd.to_numeric(pd.Series(s), errors="coerce").to_numpy(float)


# ---------------------------------------------------------
# Model predictions (for bootstrap overlay)
# ---------------------------------------------------------
def ou_predict(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return a + b * x


def tanh_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    return A * np.tanh((x_star - x) / Delta)


def rational_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    z = (x_star - x) / Delta
    return A * (z / np.sqrt(1.0 + z * z))


def _sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def hinge_plateau_predict(x: np.ndarray, a: float, b: float, A: float, x_c: float, Delta: float) -> np.ndarray:
    s = 1.0 - _sigmoid((x - x_c) / Delta)  # left=1, right=0
    return (a + b * x) * s + (-A) * (1.0 - s)


# ---------------------------------------------------------
# Plot aesthetics
# ---------------------------------------------------------
def apply_rcparams(font_family: str, font_size: float, axes_w: float, tick_w: float, tick_len: float) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": font_size,
        "axes.linewidth": axes_w,
        "xtick.major.width": tick_w,
        "ytick.major.width": tick_w,
        "xtick.major.size": tick_len,
        "ytick.major.size": tick_len,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _maybe_set_limits(ax, xlim_lo, xlim_hi, ylim_lo, ylim_hi):
    if np.isfinite(xlim_lo) or np.isfinite(xlim_hi):
        lo = xlim_lo if np.isfinite(xlim_lo) else None
        hi = xlim_hi if np.isfinite(xlim_hi) else None
        ax.set_xlim(lo, hi)
    if np.isfinite(ylim_lo) or np.isfinite(ylim_hi):
        lo = ylim_lo if np.isfinite(ylim_lo) else None
        hi = ylim_hi if np.isfinite(ylim_hi) else None
        ax.set_ylim(lo, hi)


# ---------------------------------------------------------
# Legend + layout + save helpers
# ---------------------------------------------------------
def apply_legend(ax: plt.Axes, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        ax.legend(
            loc=args.legend_loc,
            bbox_to_anchor=(args.legend_bbox_x, args.legend_bbox_y),
            frameon=bool(args.legend_frame),
            fontsize=args.legend_fontsize,
            borderaxespad=0.0,
        )
    else:
        ax.legend(
            loc=args.legend_loc,
            frameon=bool(args.legend_frame),
            fontsize=args.legend_fontsize,
        )


def apply_layout(fig: plt.Figure, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.tight_layout(rect=[0.0, 0.0, float(args.legend_pad_right), 1.0])
    else:
        fig.tight_layout()


def save_png(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.savefig(path, dpi=int(args.dpi), bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path, dpi=int(args.dpi))


def save_pdf(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.savefig(path, bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path)


# ---------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------
def discover_b_cols(drift_grid: pd.DataFrame) -> List[str]:
    return [c for c in drift_grid.columns if str(c).startswith("b_")]


def discover_resid_series(bmed: pd.DataFrame) -> List[str]:
    if "series" not in bmed.columns:
        return []
    series = bmed["series"].astype(str).unique().tolist()
    return [s for s in series if s.startswith("resid_")]


def pretty_label(name: str) -> str:
    x = str(name)
    if x.startswith("b_"):
        x = x[2:]
    if x.startswith("resid_"):
        x = x[6:]
    return x.replace("_", " ")


def parse_show_list(arg: str) -> Optional[List[str]]:
    s = (arg or "").strip()
    if not s:
        return None
    return [t.strip() for t in s.split(",") if t.strip()]


# ---------------------------------------------------------
# Fixed model order + color mapping
# ---------------------------------------------------------
MODEL_ORDER = [
    "OU_linear",
    "tanh_saturating",
    "rational_saturating",
    "hinge_plateau_smooth",
]

MODEL_COLOR_MAP = {
    "OU_linear": "C0",
    "tanh_saturating": "C1",
    "rational_saturating": "C2",
    "hinge_plateau_smooth": "C3",
}


def model_name_from_bcol(bcol: str) -> str:
    return bcol[2:] if bcol.startswith("b_") else bcol


def model_name_from_resid_series(series_name: str) -> str:
    return series_name[6:] if series_name.startswith("resid_") else series_name


def sort_models(model_names: List[str]) -> List[str]:
    rank = {m: i for i, m in enumerate(MODEL_ORDER)}
    return sorted(model_names, key=lambda x: (rank.get(x, 999), x))


def sort_b_cols(b_cols: List[str]) -> List[str]:
    return [f"b_{m}" for m in sort_models([model_name_from_bcol(c) for c in b_cols]) if f"b_{m}" in b_cols]


def sort_resid_series(resid_series: List[str]) -> List[str]:
    return [f"resid_{m}" for m in sort_models([model_name_from_resid_series(s) for s in resid_series]) if f"resid_{m}" in resid_series]


def get_model_color(model_name: str) -> str:
    return MODEL_COLOR_MAP.get(model_name, "black")


# ---------------------------------------------------------
# Bootstrap mapping
# ---------------------------------------------------------
def build_bootstrap_predictors(reps_df: pd.DataFrame) -> Dict[str, Callable[[pd.Series, np.ndarray], np.ndarray]]:
    preds: Dict[str, Callable[[pd.Series, np.ndarray], np.ndarray]] = {}
    cols = set(reps_df.columns.astype(str))

    if {"ou_a", "ou_b"}.issubset(cols):
        preds["OU_linear"] = lambda r, xs: ou_predict(xs, float(r["ou_a"]), float(r["ou_b"]))
    elif {"OU_a", "OU_b"}.issubset(cols):
        preds["OU_linear"] = lambda r, xs: ou_predict(xs, float(r["OU_a"]), float(r["OU_b"]))

    if {"log_A", "log_x_star", "log_Delta"}.issubset(cols):
        preds["tanh_saturating"] = lambda r, xs: tanh_predict(xs, float(r["log_A"]), float(r["log_x_star"]), float(r["log_Delta"]))
    if {"tanh_A", "tanh_x_star", "tanh_Delta"}.issubset(cols):
        preds["tanh_saturating"] = lambda r, xs: tanh_predict(xs, float(r["tanh_A"]), float(r["tanh_x_star"]), float(r["tanh_Delta"]))

    if {"rat_A", "rat_x_star", "rat_Delta"}.issubset(cols):
        preds["rational_saturating"] = lambda r, xs: rational_predict(xs, float(r["rat_A"]), float(r["rat_x_star"]), float(r["rat_Delta"]))

    if {"hinge_a", "hinge_b", "hinge_A", "hinge_x_c", "hinge_Delta"}.issubset(cols):
        preds["hinge_plateau_smooth"] = lambda r, xs: hinge_plateau_predict(
            xs, float(r["hinge_a"]), float(r["hinge_b"]), float(r["hinge_A"]), float(r["hinge_x_c"]), float(r["hinge_Delta"])
        )

    return preds


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--fig_dir", required=True)

    # Typography
    ap.add_argument("--font_family", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)

    # Figure geometry
    ap.add_argument("--fig_w", type=float, default=7.2)
    ap.add_argument("--fig_h", type=float, default=4.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--pdf", action="store_true")

    # Styling
    ap.add_argument("--line_w", type=float, default=1.6)
    ap.add_argument("--axline_w", type=float, default=1.0)
    ap.add_argument("--tick_w", type=float, default=1.0)
    ap.add_argument("--tick_len", type=float, default=4.0)
    ap.add_argument("--axes_w", type=float, default=1.0)
    ap.add_argument("--marker_s", type=float, default=26.0)
    ap.add_argument("--marker_alpha", type=float, default=0.85)

    # Legends
    ap.add_argument("--legend_loc", default="best")
    ap.add_argument("--legend_frame", type=int, default=0)
    ap.add_argument("--legend_fontsize", type=float, default=10.0)

    # Optional legend outside
    ap.add_argument("--legend_outside", action="store_true")
    ap.add_argument("--legend_bbox_x", type=float, default=1.02)
    ap.add_argument("--legend_bbox_y", type=float, default=1.0)
    ap.add_argument("--legend_pad_right", type=float, default=0.78)
    ap.add_argument("--legend_pad_inches", type=float, default=0.05)

    # Labels
    ap.add_argument("--xlabel", default="x0")
    ap.add_argument("--ylabel_dx", default="dx")
    ap.add_argument("--ylabel_resid", default="median residual (dx - fit)")

    # Optional axis limits
    ap.add_argument("--xlim_lo", type=float, default=np.nan)
    ap.add_argument("--xlim_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_dx_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_dx_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_resid_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_resid_hi", type=float, default=np.nan)

    # Titles
    ap.add_argument("--title_prefix", default="")
    ap.add_argument("--no_titles", action="store_true")

    # Which models to show
    ap.add_argument(
        "--show_models",
        default="",
        help="Comma-separated model names (e.g. OU_linear,tanh_saturating). Empty=all."
    )

    # Bootstrap overlay
    ap.add_argument("--bootstrap_draw", type=int, default=80)
    ap.add_argument("--bootstrap_alpha", type=float, default=0.12)
    ap.add_argument("--bootstrap_line_w", type=float, default=0.8)

    return ap.parse_args()


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def main() -> None:
    args = parse_args()

    data_dir = Path(args.data_dir)
    fig_dir = ensure_dir(Path(args.fig_dir))

    apply_rcparams(args.font_family, args.font_size, args.axes_w, args.tick_w, args.tick_len)

    grid_path = data_dir / "drift_function_grid.csv"
    bmed_path = data_dir / "binned_medians.csv"
    boot_path = data_dir / "bootstrap_reps.csv"

    if not grid_path.exists():
        raise SystemExit(f"Missing: {grid_path}")
    if not bmed_path.exists():
        raise SystemExit(f"Missing: {bmed_path}")

    drift_grid = pd.read_csv(grid_path)
    bmed = pd.read_csv(bmed_path)

    if "x" not in drift_grid.columns:
        raise SystemExit("drift_function_grid.csv must contain column 'x'.")
    if "series" not in bmed.columns:
        raise SystemExit("binned_medians.csv must contain column 'series'.")
    for c in ["x0_mid", "y_med"]:
        if c not in bmed.columns:
            raise SystemExit(f"binned_medians.csv must contain column '{c}'.")

    has_boot = boot_path.exists()
    reps_df = pd.read_csv(boot_path) if has_boot else pd.DataFrame()

    # Meta for titles
    dt = int(drift_grid["dt"].iloc[0]) if "dt" in drift_grid.columns else None
    obs_class = str(drift_grid["obs_class"].iloc[0]) if "obs_class" in drift_grid.columns else ""
    weight_col = str(drift_grid["weight_col"].iloc[0]) if "weight_col" in drift_grid.columns else ""
    n = int(drift_grid["n"].iloc[0]) if "n" in drift_grid.columns else None

    tcls = obs_class if obs_class else "ALL"
    tw = weight_col if weight_col else "unweighted"
    base_title = f"dt={dt}, class={tcls}, {tw}, n={n}"
    if args.title_prefix.strip():
        base_title = f"{args.title_prefix.strip()}  ({base_title})"

    # Discover + filter + sort drift models
    b_cols = discover_b_cols(drift_grid)
    if not b_cols:
        raise SystemExit("No drift curves found (need columns starting with 'b_').")

    show = parse_show_list(args.show_models)
    if show is not None:
        keep = []
        for model_name in show:
            c = f"b_{model_name}" if not model_name.startswith("b_") else model_name
            if c in b_cols:
                keep.append(c)
        if not keep:
            raise SystemExit(f"--show_models requested {show} but none found among: {b_cols}")
        b_cols = keep

    b_cols = sort_b_cols(b_cols)

    # Discover + filter + sort residual series
    resid_series = discover_resid_series(bmed)
    if show is not None:
        resid_series = [f"resid_{m}" for m in show if f"resid_{m}" in resid_series]
    resid_series = sort_resid_series(resid_series)

    xs = _to_float_array(drift_grid["x"])

    # ---------------------------------------------------------
    # Figure 1: Drift compare
    # ---------------------------------------------------------
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(1, 1, 1)

    bdx = bmed[bmed["series"].astype(str) == "dx"].copy()
    if not bdx.empty:
        ax.scatter(
            _as_float_col(bdx, "x0_mid"),
            _as_float_col(bdx, "y_med"),
            s=args.marker_s,
            alpha=args.marker_alpha,
            color="black",
            label="binned median(dx)",
            zorder=3,
        )

    for c in b_cols:
        model_name = model_name_from_bcol(c)
        y = _to_float_array(drift_grid[c])
        mask = np.isfinite(xs) & np.isfinite(y)
        ax.plot(
            xs[mask],
            y[mask],
            linewidth=args.line_w,
            color=get_model_color(model_name),
            label=pretty_label(c),
            zorder=2,
        )

    ax.axhline(0, linewidth=args.axline_w, color="black")
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel_dx)
    _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_dx_lo, args.ylim_dx_hi)

    apply_legend(ax, args)
    if not args.no_titles:
        ax.set_title(f"Finite-time drift fit: {base_title}")

    apply_layout(fig, args)
    save_png(fig, fig_dir / "fig_drift_compare.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "fig_drift_compare.pdf", args)
    plt.close(fig)

    # ---------------------------------------------------------
    # Figure 2: Residual median by x0-bin
    # ---------------------------------------------------------
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(1, 1, 1)

    if not resid_series:
        ax.text(0.5, 0.5, "No resid_* series found in binned_medians.csv", ha="center", va="center")
    else:
        for sname in resid_series:
            part = bmed[bmed["series"].astype(str) == sname].copy()
            if part.empty:
                continue

            model_name = model_name_from_resid_series(sname)
            x_part = _as_float_col(part, "x0_mid")
            y_part = _as_float_col(part, "y_med")
            mask = np.isfinite(x_part) & np.isfinite(y_part)

            ax.plot(
                x_part[mask],
                y_part[mask],
                linewidth=args.line_w,
                color=get_model_color(model_name),
                label=pretty_label(sname),
            )

    ax.axhline(0, linewidth=args.axline_w, color="black")
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel_resid)
    _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_resid_lo, args.ylim_resid_hi)

    apply_legend(ax, args)
    if not args.no_titles:
        ax.set_title(f"Median residual by x0-bin: {base_title}")

    apply_layout(fig, args)
    save_png(fig, fig_dir / "fig_resid_median_by_bin.png", args)
    save_png(fig, fig_dir / "fig_resid_compare.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "fig_resid_median_by_bin.pdf", args)
        save_pdf(fig, fig_dir / "fig_resid_compare.pdf", args)
    plt.close(fig)

    # ---------------------------------------------------------
    # Figure 3: Bootstrap overlay
    # ---------------------------------------------------------
    if has_boot and (not reps_df.empty):
        fig = plt.figure(figsize=(args.fig_w, args.fig_h))
        ax = fig.add_subplot(1, 1, 1)

        if not bdx.empty:
            ax.scatter(
                _as_float_col(bdx, "x0_mid"),
                _as_float_col(bdx, "y_med"),
                s=args.marker_s,
                alpha=args.marker_alpha,
                color="black",
                label="binned median(dx)",
                zorder=3,
            )

        preds = build_bootstrap_predictors(reps_df)

        wanted_models = [model_name_from_bcol(c) for c in b_cols]
        preds = {k: v for k, v in preds.items() if k in wanted_models}

        max_draw = int(min(max(args.bootstrap_draw, 0), len(reps_df)))
        if max_draw > 0 and preds:
            draw_idx = np.linspace(0, len(reps_df) - 1, max_draw).astype(int)
            for model_name in sort_models(list(preds.keys())):
                fn = preds[model_name]
                col = get_model_color(model_name)
                for j in draw_idx:
                    r = reps_df.iloc[j]
                    try:
                        yb = fn(r, xs)
                        yb = np.asarray(yb, dtype=float)
                        mask = np.isfinite(xs) & np.isfinite(yb)
                        if np.any(mask):
                            ax.plot(
                                xs[mask],
                                yb[mask],
                                linewidth=args.bootstrap_line_w,
                                alpha=args.bootstrap_alpha,
                                color=col,
                                zorder=1,
                            )
                    except Exception:
                        continue

        for c in b_cols:
            model_name = model_name_from_bcol(c)
            y = _to_float_array(drift_grid[c])
            mask = np.isfinite(xs) & np.isfinite(y)
            ax.plot(
                xs[mask],
                y[mask],
                linewidth=args.line_w,
                color=get_model_color(model_name),
                label=pretty_label(c),
                zorder=2,
            )

        ax.axhline(0, linewidth=args.axline_w, color="black")
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel_dx)
        _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_dx_lo, args.ylim_dx_hi)

        apply_legend(ax, args)
        if not args.no_titles:
            ax.set_title(f"Bootstrap overlay: {base_title}")

        apply_layout(fig, args)
        save_png(fig, fig_dir / "fig_bootstrap_overlay.png", args)
        if args.pdf:
            save_pdf(fig, fig_dir / "fig_bootstrap_overlay.pdf", args)
        plt.close(fig)

    print("✅ Completed SINGLE-FIGURE plotting (multi-model, consistent colors).")
    print(f"[PATH] fig_dir: {fig_dir.resolve()}")
    print("[OUT]  fig_drift_compare.png")
    print("[OUT]  fig_resid_median_by_bin.png  (also saved as fig_resid_compare.png)")
    if has_boot and (not reps_df.empty):
        print("[OUT]  fig_bootstrap_overlay.png")
    else:
        print("[NOTE] bootstrap_reps.csv missing or empty -> no bootstrap overlay.")


if __name__ == "__main__":
    main()