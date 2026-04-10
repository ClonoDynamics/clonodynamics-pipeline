#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fp_alpha_diagnostics_v5.py

Diagnostics for Fokker–Planck fits across detectability thresholds alpha.

Input
-----
CSV manifest with columns:
    alpha,grid_csv,summary_json

Expected grid columns
---------------------
Required:
    x
Preferred:
    p_emp
    p_closed
    ratio_pemp_over_pclosed
    J_emp
    S_emp

Fallbacks are implemented for common alternative names.

Outputs
-------
Summary tables:
    fp_alpha_summary.csv
    fp_alpha_summary_pretty.csv
    fp_alpha_report.txt

Metrics:
    fig_fp_alpha_metrics.{fmt}

Single-panel plots for each alpha:
    fig_fp_alpha_pdf_p_<alpha>.{fmt}
    fig_fp_alpha_ratio_p_<alpha>.{fmt}
    fig_fp_alpha_flux_p_<alpha>.{fmt}
    fig_fp_alpha_source_p_<alpha>.{fmt}

Multipanel (if exactly 4 alphas are present):
    fig_fp_alpha_multipanel_4x4.{fmt}

Style
-----
Arial, 10 pt.

Example
-------
python fp_alpha_diagnostics_v5.py \
  --pairs_csv ./alpha_pairs.csv \
  --out_dir ./fp_alpha_diagnostics \
  --fmt pdf \
  --x_min -10.8 \
  --x_max -6.2 \
  --ratio_ymin 0.6 --ratio_ymax 1.8 \
  --flux_ymin -0.15 --flux_ymax 0.15 \
  --source_ymin -2 --source_ymax 2 \
  --ratio_smooth_sigma_pts 1.0 \
  --flux_smooth_sigma_pts 1.5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================================================
# CLI
# =========================================================
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="FP alpha diagnostics: singles + 4x4 multipanel.")

    ap.add_argument("--pairs_csv", required=True,
                    help="CSV manifest with columns: alpha,grid_csv,summary_json")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--fmt", default="pdf", choices=["pdf", "png"], help="Figure format")

    # x filtering
    ap.add_argument("--x_min", type=float, default=None, help="Optional lower x bound")
    ap.add_argument("--x_max", type=float, default=None, help="Optional upper x bound")

    # smoothing
    ap.add_argument("--ratio_smooth_sigma_pts", type=float, default=1.0,
                    help="Gaussian smoothing sigma in points for ratio")
    ap.add_argument("--flux_smooth_sigma_pts", type=float, default=1.5,
                    help="Gaussian smoothing sigma in points for flux and source")

    # manual y-limits
    ap.add_argument("--ratio_ymin", type=float, default=None)
    ap.add_argument("--ratio_ymax", type=float, default=None)

    ap.add_argument("--flux_ymin", type=float, default=None)
    ap.add_argument("--flux_ymax", type=float, default=None)

    ap.add_argument("--source_ymin", type=float, default=None)
    ap.add_argument("--source_ymax", type=float, default=None)

    # robust auto-limits
    ap.add_argument("--auto_qlo", type=float, default=0.02,
                    help="Lower quantile for robust automatic y-limits")
    ap.add_argument("--auto_qhi", type=float, default=0.98,
                    help="Upper quantile for robust automatic y-limits")
    ap.add_argument("--auto_pad_frac", type=float, default=0.05,
                    help="Padding fraction added to automatic y-limits")

    # figure sizes
    ap.add_argument("--single_width", type=float, default=3.1)
    ap.add_argument("--single_height", type=float, default=3.1)
    ap.add_argument("--multi_width", type=float, default=12.8)
    ap.add_argument("--multi_height", type=float, default=11.0)
    ap.add_argument("--dpi", type=int, default=300)

    return ap


# =========================================================
# Style
# =========================================================
def setup_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.titlesize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


# =========================================================
# Helpers
# =========================================================
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def alpha_to_label(alpha: float) -> str:
    if alpha >= 0.01:
        return f"p_{alpha:g}"
    return f"p_{alpha:.3f}".rstrip("0").rstrip(".")


def alpha_to_filename(alpha: float) -> str:
    s = f"{alpha:.6g}"
    return f"p_{s.replace('.', '_')}"


def normalize_ylim(ymin: Optional[float], ymax: Optional[float]) -> Optional[Tuple[float, float]]:
    if ymin is None or ymax is None:
        return None
    if not (np.isfinite(ymin) and np.isfinite(ymax)):
        return None
    if ymax <= ymin:
        return None
    return (float(ymin), float(ymax))


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


# =========================================================
# Column mapping
# =========================================================
def resolve_columns(g: pd.DataFrame) -> Dict[str, Optional[str]]:
    return {
        "x": pick_first_present(g, ["x", "x_mid", "log_freq", "logf"]),
        "p_emp": pick_first_present(g, ["p_emp", "p_empirical", "pdf_emp", "p_observed"]),
        "p_closed": pick_first_present(g, ["p_closed", "p_model", "pdf_model", "p_stationary"]),
        "ratio": pick_first_present(g, ["ratio_pemp_over_pclosed", "ratio_emp_over_model", "ratio", "density_ratio"]),
        "J": pick_first_present(g, ["J_emp", "J_model", "J", "flux", "probability_flux"]),
        "S": pick_first_present(g, ["S_emp", "S_model", "S", "source", "source_term"]),
    }


# =========================================================
# Data loading
# =========================================================
def load_dataset(alpha: float, grid_csv: Path, summary_json: Path, args) -> Dict:
    g = pd.read_csv(grid_csv)
    meta = load_json(summary_json)
    col = resolve_columns(g)

    if col["x"] is None:
        raise ValueError(f"[alpha={alpha}] Missing x column in grid: {grid_csv}")
    if col["p_emp"] is None:
        raise ValueError(f"[alpha={alpha}] Missing empirical density column in grid: {grid_csv}")
    if col["p_closed"] is None:
        raise ValueError(f"[alpha={alpha}] Missing stationary/model density column in grid: {grid_csv}")
    if col["J"] is None:
        raise ValueError(f"[alpha={alpha}] Missing flux column in grid: {grid_csv}")
    if col["S"] is None:
        raise ValueError(f"[alpha={alpha}] Missing source column in grid: {grid_csv}")

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

    m = finite_mask(x, p_emp, p_closed, ratio, J, S)
    x = x[m]
    p_emp = p_emp[m]
    p_closed = p_closed[m]
    ratio = ratio[m]
    J = J[m]
    S = S[m]

    x, [p_emp, p_closed, ratio, J, S] = apply_x_filter(
        x, [p_emp, p_closed, ratio, J, S], args.x_min, args.x_max
    )

    if x.size == 0:
        raise ValueError(f"[alpha={alpha}] No rows remain after x filtering.")

    ratio_s = smooth_series(ratio, args.ratio_smooth_sigma_pts)
    J_s = smooth_series(J, args.flux_smooth_sigma_pts)
    S_s = smooth_series(S, args.flux_smooth_sigma_pts)

    return {
        "alpha": float(alpha),
        "alpha_label": alpha_to_label(float(alpha)),
        "alpha_file": alpha_to_filename(float(alpha)),
        "grid_csv": str(grid_csv),
        "summary_json": str(summary_json),
        "x": x,
        "p_emp": p_emp,
        "p_closed": p_closed,
        "ratio": ratio_s,
        "J": J_s,
        "S": S_s,
        "JS": safe_float(meta.get("JS", np.nan)),
        "KS": safe_float(meta.get("KS", np.nan)),
        "tag": meta.get("tag", ""),
    }


def resolve_global_ylims(
    datasets: List[Dict],
    ratio_ylim_user: Optional[Tuple[float, float]],
    flux_ylim_user: Optional[Tuple[float, float]],
    source_ylim_user: Optional[Tuple[float, float]],
    auto_qlo: float,
    auto_qhi: float,
    auto_pad_frac: float,
) -> Dict[str, Optional[Tuple[float, float]]]:
    pdf_ylim = robust_limits_from_arrays(
        [d["p_emp"] for d in datasets] + [d["p_closed"] for d in datasets],
        qlo=0.0,
        qhi=1.0,
        symmetric=False,
        pad_frac=auto_pad_frac,
    )

    ratio_ylim = ratio_ylim_user
    if ratio_ylim is None:
        ratio_ylim = robust_limits_from_arrays(
            [d["ratio"] for d in datasets],
            qlo=auto_qlo,
            qhi=auto_qhi,
            symmetric=False,
            pad_frac=auto_pad_frac,
        )

    flux_ylim = flux_ylim_user
    if flux_ylim is None:
        flux_ylim = robust_limits_from_arrays(
            [d["J"] for d in datasets],
            qlo=auto_qlo,
            qhi=auto_qhi,
            symmetric=True,
            pad_frac=auto_pad_frac,
        )

    source_ylim = source_ylim_user
    if source_ylim is None:
        source_ylim = robust_limits_from_arrays(
            [d["S"] for d in datasets],
            qlo=auto_qlo,
            qhi=auto_qhi,
            symmetric=True,
            pad_frac=auto_pad_frac,
        )

    return {
        "pdf": pdf_ylim,
        "ratio": ratio_ylim,
        "flux": flux_ylim,
        "source": source_ylim,
    }


# =========================================================
# Tables / report
# =========================================================
def save_summary_table(datasets: List[Dict], out_dir: Path) -> None:
    rows = []
    for d in datasets:
        x = d["x"]
        J = d["J"]
        S = d["S"]

        rows.append({
            "alpha": d["alpha"],
            "alpha_label": d["alpha_label"],
            "JS": d["JS"],
            "KS": d["KS"],
            "int_abs_J": float(np.trapz(np.abs(J), x)) if len(x) > 1 else np.nan,
            "int_abs_S": float(np.trapz(np.abs(S), x)) if len(x) > 1 else np.nan,
            "mean_abs_J": float(np.mean(np.abs(J))) if len(x) else np.nan,
            "mean_abs_S": float(np.mean(np.abs(S))) if len(x) else np.nan,
            "n_grid": int(len(x)),
            "grid_csv": d["grid_csv"],
            "summary_json": d["summary_json"],
            "tag": d["tag"],
        })

    df = pd.DataFrame(rows).sort_values("alpha").reset_index(drop=True)
    df.to_csv(out_dir / "fp_alpha_summary.csv", index=False)

    pretty = df.copy()
    for c in ["alpha", "JS", "KS", "int_abs_J", "int_abs_S", "mean_abs_J", "mean_abs_S"]:
        pretty[c] = pretty[c].map(lambda v: f"{v:.6g}" if pd.notna(v) else "")
    pretty.to_csv(out_dir / "fp_alpha_summary_pretty.csv", index=False)


def write_report(datasets: List[Dict], out_dir: Path, args, ylims: Dict[str, Optional[Tuple[float, float]]]) -> None:
    lines = [
        "FP alpha diagnostics report",
        "===========================",
        "",
        f"Input pairs_csv: {args.pairs_csv}",
        f"Output dir: {args.out_dir}",
        f"Format: {args.fmt}",
        "",
        "Global settings",
        f"  x_min = {args.x_min}",
        f"  x_max = {args.x_max}",
        f"  ratio_smooth_sigma_pts = {args.ratio_smooth_sigma_pts}",
        f"  flux_smooth_sigma_pts = {args.flux_smooth_sigma_pts}",
        f"  ratio_ymin/ymax (user) = {args.ratio_ymin}, {args.ratio_ymax}",
        f"  flux_ymin/ymax (user) = {args.flux_ymin}, {args.flux_ymax}",
        f"  source_ymin/ymax (user) = {args.source_ymin}, {args.source_ymax}",
        "",
        "Resolved common y-limits",
        f"  pdf    = {ylims['pdf']}",
        f"  ratio  = {ylims['ratio']}",
        f"  flux   = {ylims['flux']}",
        f"  source = {ylims['source']}",
        "",
        "Datasets",
    ]

    for d in datasets:
        lines.append(
            f"  alpha={d['alpha']:.6g}  tag={d['tag']}  n_grid={len(d['x'])}  JS={d['JS']:.6g}  KS={d['KS']:.6g}"
        )
        lines.append(f"    grid:    {d['grid_csv']}")
        lines.append(f"    summary: {d['summary_json']}")

    (out_dir / "fp_alpha_report.txt").write_text("\n".join(lines), encoding="utf-8")


# =========================================================
# Plotting
# =========================================================
def plot_metrics(datasets: List[Dict], out_dir: Path, fmt: str, dpi: int) -> None:
    alphas = [d["alpha"] for d in datasets]
    JS = [d["JS"] for d in datasets]
    KS = [d["KS"] for d in datasets]

    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    ax.plot(alphas, JS, "-o", label="JS")
    ax.plot(alphas, KS, "-o", label="KS")
    ax.set_xscale("log")
    ax.set_xlabel("alpha")
    ax.set_ylabel("metric")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_metrics.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_pdf(
    d: Dict,
    out_dir: Path,
    fmt: str,
    dpi: int,
    w: float,
    h: float,
    pdf_ylim: Optional[Tuple[float, float]],
) -> None:
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(d["x"], d["p_emp"], label="empirical")
    ax.plot(d["x"], d["p_closed"], label="stationary")
    ax.set_xlabel("log-frequency")
    ax.set_ylabel("density")
    ax.set_title(d["alpha_label"])
    maybe_set_ylim(ax, pdf_ylim)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_pdf_{d['alpha_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_ratio(
    d: Dict,
    out_dir: Path,
    fmt: str,
    dpi: int,
    w: float,
    h: float,
    ratio_ylim: Optional[Tuple[float, float]],
) -> None:
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(d["x"], d["ratio"])
    ax.axhline(1.0, ls="--", lw=1.0)
    ax.set_xlabel("log-frequency")
    ax.set_ylabel("emp/model")
    ax.set_title(d["alpha_label"])
    maybe_set_ylim(ax, ratio_ylim)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_ratio_{d['alpha_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_flux(
    d: Dict,
    out_dir: Path,
    fmt: str,
    dpi: int,
    w: float,
    h: float,
    flux_ylim: Optional[Tuple[float, float]],
) -> None:
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(d["x"], d["J"])
    ax.axhline(0.0, ls="--", lw=1.0)
    ax.set_xlabel("log-frequency")
    ax.set_ylabel("J(x)")
    ax.set_title(d["alpha_label"])
    maybe_set_ylim(ax, flux_ylim)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_flux_{d['alpha_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_single_source(
    d: Dict,
    out_dir: Path,
    fmt: str,
    dpi: int,
    w: float,
    h: float,
    source_ylim: Optional[Tuple[float, float]],
) -> None:
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot(d["x"], d["S"])
    ax.axhline(0.0, ls="--", lw=1.0)
    ax.set_xlabel("log-frequency")
    ax.set_ylabel("S(x)")
    ax.set_title(d["alpha_label"])
    maybe_set_ylim(ax, source_ylim)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_source_{d['alpha_file']}.{fmt}", dpi=dpi)
    plt.close(fig)


def plot_multipanel_4x4(
    datasets: List[Dict],
    out_dir: Path,
    fmt: str,
    dpi: int,
    fig_w: float,
    fig_h: float,
    pdf_ylim: Optional[Tuple[float, float]],
    ratio_ylim: Optional[Tuple[float, float]],
    flux_ylim: Optional[Tuple[float, float]],
    source_ylim: Optional[Tuple[float, float]],
) -> None:
    if len(datasets) != 4:
        raise ValueError(f"4x4 multipanel requires exactly 4 alpha values; found {len(datasets)}")

    fig, axes = plt.subplots(4, 4, figsize=(fig_w, fig_h), sharex=False)

    for col, d in enumerate(datasets):
        # PDFs
        ax = axes[0, col]
        ax.plot(d["x"], d["p_emp"], label="empirical")
        ax.plot(d["x"], d["p_closed"], label="stationary")
        ax.set_title(d["alpha_label"])
        if col == 0:
            ax.set_ylabel("density")
        maybe_set_ylim(ax, pdf_ylim)
        if col == 3:
            ax.legend(frameon=False, loc="upper right")

        # ratio
        ax = axes[1, col]
        ax.plot(d["x"], d["ratio"])
        ax.axhline(1.0, ls="--", lw=1.0)
        if col == 0:
            ax.set_ylabel("emp/model")
        maybe_set_ylim(ax, ratio_ylim)

        # flux
        ax = axes[2, col]
        ax.plot(d["x"], d["J"])
        ax.axhline(0.0, ls="--", lw=1.0)
        if col == 0:
            ax.set_ylabel("J(x)")
        maybe_set_ylim(ax, flux_ylim)

        # source
        ax = axes[3, col]
        ax.plot(d["x"], d["S"])
        ax.axhline(0.0, ls="--", lw=1.0)
        if col == 0:
            ax.set_ylabel("S(x)")
        ax.set_xlabel("log-frequency")
        maybe_set_ylim(ax, source_ylim)

    fig.tight_layout()
    fig.savefig(out_dir / f"fig_fp_alpha_multipanel_4x4.{fmt}", dpi=dpi)
    plt.close(fig)


# =========================================================
# Main
# =========================================================
def main() -> None:
    args = build_argparser().parse_args()
    setup_matplotlib()

    out_dir = ensure_dir(Path(args.out_dir).resolve())

    ratio_ylim_user = normalize_ylim(args.ratio_ymin, args.ratio_ymax)
    flux_ylim_user = normalize_ylim(args.flux_ymin, args.flux_ymax)
    source_ylim_user = normalize_ylim(args.source_ymin, args.source_ymax)

    pairs = pd.read_csv(Path(args.pairs_csv).resolve())
    required = {"alpha", "grid_csv", "summary_json"}
    missing = required - set(pairs.columns)
    if missing:
        raise ValueError(f"pairs_csv missing required columns: {sorted(missing)}")

    pairs = pairs.copy()
    pairs["alpha"] = pd.to_numeric(pairs["alpha"], errors="coerce")
    if pairs["alpha"].isna().any():
        raise ValueError("Column 'alpha' contains non-numeric values.")
    pairs = pairs.sort_values("alpha").reset_index(drop=True)

    datasets: List[Dict] = []
    for _, row in pairs.iterrows():
        alpha = float(row["alpha"])
        grid_csv = Path(str(row["grid_csv"])).expanduser().resolve()
        summary_json = Path(str(row["summary_json"])).expanduser().resolve()

        if not grid_csv.exists():
            raise FileNotFoundError(f"Grid CSV not found: {grid_csv}")
        if not summary_json.exists():
            raise FileNotFoundError(f"Summary JSON not found: {summary_json}")

        datasets.append(load_dataset(alpha, grid_csv, summary_json, args))

    ylims = resolve_global_ylims(
        datasets=datasets,
        ratio_ylim_user=ratio_ylim_user,
        flux_ylim_user=flux_ylim_user,
        source_ylim_user=source_ylim_user,
        auto_qlo=args.auto_qlo,
        auto_qhi=args.auto_qhi,
        auto_pad_frac=args.auto_pad_frac,
    )

    save_summary_table(datasets, out_dir)
    write_report(datasets, out_dir, args, ylims)
    plot_metrics(datasets, out_dir, args.fmt, args.dpi)

    for d in datasets:
        plot_single_pdf(
            d, out_dir, args.fmt, args.dpi,
            args.single_width, args.single_height,
            ylims["pdf"]
        )
        plot_single_ratio(
            d, out_dir, args.fmt, args.dpi,
            args.single_width, args.single_height,
            ylims["ratio"]
        )
        plot_single_flux(
            d, out_dir, args.fmt, args.dpi,
            args.single_width, args.single_height,
            ylims["flux"]
        )
        plot_single_source(
            d, out_dir, args.fmt, args.dpi,
            args.single_width, args.single_height,
            ylims["source"]
        )

    if len(datasets) == 4:
        plot_multipanel_4x4(
            datasets=datasets,
            out_dir=out_dir,
            fmt=args.fmt,
            dpi=args.dpi,
            fig_w=args.multi_width,
            fig_h=args.multi_height,
            pdf_ylim=ylims["pdf"],
            ratio_ylim=ylims["ratio"],
            flux_ylim=ylims["flux"],
            source_ylim=ylims["source"],
        )
    else:
        print(f"[INFO] Multipanel 4x4 not generated because number of alphas is {len(datasets)} (need exactly 4).")

    print(f"[DONE] Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()