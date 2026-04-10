#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_target_trajectory_and_dual_calibration.py

Read one or more per_clone_denoised_subject.csv files, reconstruct the trajectory
of a user-specified target clonotype (aaSeqCDR3), display it together with a
representative set of non-target clonotypes, and estimate two calibration functions
linking molecular frequencies (freq_geo_x) to cellular target frequencies:

1) calibration using all target time points
2) calibration using only target time points with observable == TRUE

Non-target clonotypes are selected so as to be distributed across the full
frequency range, rather than being concentrated only among the most abundant ones.

Main outputs
------------
1) target_plus_background.png
   Target trajectory + representative non-target clonotypes

2) target_calibration_all_points.png
   Calibration using all target time points

3) target_calibration_observable_only.png
   Calibration using only target time points with observable == TRUE

4) CSV tables:
   - selected_clones.csv
   - extracted_trajectories.csv
   - observability_matrix.csv
   - target_calibration_points_all.csv
   - target_calibration_fit_all.csv
   - target_calibration_curve_all.csv
   - target_calibration_points_observable.csv
   - target_calibration_fit_observable.csv
   - target_calibration_curve_observable.csv
"""

import argparse
import glob
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# =========================
# Global matplotlib style
# =========================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 10,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})

DEFAULT_TARGET = "CASSFSTCSANYGYTF"


# =========================
# Utilities
# =========================
def as_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    ss = s.astype(str).str.strip().str.lower()
    return ss.isin(["true", "t", "1", "yes", "y"])


def collect_input_files(inputs: List[str]) -> List[Path]:
    files: List[Path] = []

    for a in inputs:
        p = Path(a)
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            found = sorted(
                Path(x) for x in glob.glob(
                    str(p / "**" / "per_clone_denoised_subject.csv"),
                    recursive=True
                )
            )
            if not found:
                found = sorted(
                    Path(x) for x in glob.glob(
                        str(p / "**" / "per_clone_denoised_subject*"),
                        recursive=True
                    )
                )
            if not found:
                raise RuntimeError(f"No per_clone_denoised_subject.csv files found in directory: {p}")
            files.extend(found)
        else:
            raise RuntimeError(f"Input path not found: {p}")

    seen = set()
    uniq: List[Path] = []
    for f in files:
        if f not in seen:
            uniq.append(f)
            seen.add(f)

    if not uniq:
        raise RuntimeError("No input files found.")
    return uniq


def get_frequency_column(df: pd.DataFrame, srcfile: Path) -> pd.Series:
    """
    Return a linear-scale frequency series.

    Priority:
      1) freq_geo_x
      2) frequency
      3) freq
      4) f_denoised
      5) posterior_mean
      6) freq_denoised
      7) f_post
      8) f_hat
      9) mean_freq
      10) log_freq / log10_freq -> 10**(...)
      11) ln_freq / loge_freq -> exp(...)
    """
    linear_candidates = [
        "freq_geo_x",
        "frequency",
        "freq",
        "f_denoised",
        "posterior_mean",
        "freq_denoised",
        "f_post",
        "f_hat",
        "mean_freq",
    ]

    for c in linear_candidates:
        if c in df.columns:
            return pd.to_numeric(df[c], errors="coerce")

    log10_candidates = ["log_freq", "log10_freq"]
    for c in log10_candidates:
        if c in df.columns:
            return 10 ** pd.to_numeric(df[c], errors="coerce")

    ln_candidates = ["ln_freq", "loge_freq"]
    for c in ln_candidates:
        if c in df.columns:
            return np.exp(pd.to_numeric(df[c], errors="coerce"))

    raise RuntimeError(
        f"Missing frequency-like column in {srcfile}\n"
        f"Available columns are:\n{list(df.columns)}"
    )


def infer_time_column(df: pd.DataFrame, srcfile: Path) -> str:
    candidates = ["time", "t", "week", "timepoint"]
    for c in candidates:
        if c in df.columns:
            return c
    raise RuntimeError(
        f"Missing time column in {srcfile}. Tried candidates: {candidates}\n"
        f"Available columns: {list(df.columns)}"
    )


def infer_subject_column(df: pd.DataFrame) -> Optional[str]:
    for c in ["subject", "donor", "id"]:
        if c in df.columns:
            return c
    return None


def infer_observable_column(df: pd.DataFrame) -> Optional[str]:
    for c in ["observable", "detectable", "is_detectable", "obs"]:
        if c in df.columns:
            return c
    return None


def parse_subject_filter(subject_arg: Optional[str]):
    if subject_arg is None:
        return None
    return str(subject_arg).strip()


def parse_target_cell_freqs(s: str) -> np.ndarray:
    vals = [x.strip() for x in s.split(",") if x.strip() != ""]
    if len(vals) == 0:
        raise ValueError("--target-cell-freqs is empty.")
    arr = np.array([float(x) for x in vals], dtype=float)
    if np.any(arr <= 0):
        raise ValueError("All target cellular frequencies must be > 0.")
    return arr


# =========================
# Non-target selection
# =========================
def select_non_target_clones_stratified(
    df: pd.DataFrame,
    target_norm: str,
    n_non_target: int,
    freq_min: Optional[float],
    freq_max: Optional[float],
    time_col: str,
    n_bins: int = 5,
) -> List[str]:
    """
    Select non-target clonotypes approximately evenly across the frequency range.

    Rules:
      - clonotype must be present at all time points
      - optional filters on mean frequency
      - candidates are stratified across log10(mean frequency) bins
      - approximately equal numbers are drawn from each bin
      - within each bin, picks are spread across the local rank
    """
    times = sorted(df[time_col].dropna().unique().tolist())
    n_times = len(times)
    if n_times == 0:
        return []

    # clonotypes present at all time points
    nunique_time = df.groupby("aaSeqCDR3_norm")[time_col].nunique()
    candidates = nunique_time[nunique_time == n_times].index.tolist()

    mean_freq = (
        df.groupby("aaSeqCDR3_norm")["frequency"]
        .mean()
        .rename("mean_freq")
        .reset_index()
    )

    cand_df = pd.DataFrame({"aaSeqCDR3_norm": candidates}).merge(
        mean_freq, on="aaSeqCDR3_norm", how="left"
    )
    cand_df["mean_freq"] = cand_df["mean_freq"].fillna(0.0)

    # remove target
    cand_df = cand_df[cand_df["aaSeqCDR3_norm"] != target_norm].copy()

    # optional range filter
    if freq_min is not None:
        cand_df = cand_df[cand_df["mean_freq"] >= float(freq_min)]
    if freq_max is not None:
        cand_df = cand_df[cand_df["mean_freq"] <= float(freq_max)]

    # positive frequencies only for log binning
    cand_df = cand_df[cand_df["mean_freq"] > 0].copy()
    if cand_df.empty:
        return []

    # if very few candidates, return them all
    if len(cand_df) <= n_non_target:
        return cand_df.sort_values("mean_freq", ascending=False)["aaSeqCDR3_norm"].tolist()

    cand_df["log10_mean_freq"] = np.log10(cand_df["mean_freq"])

    qbins = min(n_bins, cand_df["aaSeqCDR3_norm"].nunique())
    cand_df["freq_bin"] = pd.qcut(
        cand_df["log10_mean_freq"],
        q=qbins,
        duplicates="drop"
    )

    grouped_bins = [g.copy() for _, g in cand_df.groupby("freq_bin", observed=True)]
    n_actual_bins = len(grouped_bins)
    if n_actual_bins == 0:
        return []

    base = n_non_target // n_actual_bins
    remainder = n_non_target % n_actual_bins

    selected = []

    for i, g in enumerate(grouped_bins):
        take = base + (1 if i < remainder else 0)
        g = g.sort_values("mean_freq", ascending=False).reset_index(drop=True)

        if len(g) <= take:
            chosen = g
        else:
            idx = np.linspace(0, len(g) - 1, num=take, dtype=int)
            chosen = g.iloc[idx]

        selected.extend(chosen["aaSeqCDR3_norm"].tolist())

    # de-duplicate preserving order
    seen = set()
    out = []
    for c in selected:
        if c not in seen:
            out.append(c)
            seen.add(c)

    return out[:n_non_target]


# =========================
# Calibration
# =========================
def fit_loglog_calibration(freq_mol: np.ndarray, freq_cell: np.ndarray) -> Tuple[float, float, float]:
    """
    Fit:
        log10(freq_cell) = a + b * log10(freq_mol)
    """
    x = np.log10(freq_mol)
    y = np.log10(freq_cell)

    b, a = np.polyfit(x, y, 1)
    yhat = a + b * x

    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = np.nan if ss_tot == 0 else 1 - ss_res / ss_tot

    return a, b, r2


def calibration_function(freq_mol: np.ndarray, a: float, b: float) -> np.ndarray:
    return (10 ** a) * (np.asarray(freq_mol) ** b)


def build_calibration_summary(
    target_df: pd.DataFrame,
    target_cell_freqs: np.ndarray,
    observable_only: bool = False
) -> pd.DataFrame:
    full_times = sorted(target_df["time"].dropna().unique().tolist())

    if len(full_times) != len(target_cell_freqs):
        raise RuntimeError(
            f"Number of provided target cellular frequencies ({len(target_cell_freqs)}) "
            f"does not match the number of sorted target time points ({len(full_times)}).\n"
            f"Target times are: {full_times}"
        )

    time_to_cell = {t: f for t, f in zip(full_times, target_cell_freqs)}

    work = target_df[target_df["observable_norm"] == True].copy() if observable_only else target_df.copy()

    if work.empty:
        raise RuntimeError("No target points available for this calibration mode.")

    summary = (
        work.groupby("time", as_index=False)
        .agg(
            freq_mol=("frequency", "mean"),
            observable_any=("observable_norm", "max"),
            n_rows=("frequency", "size"),
        )
        .sort_values("time")
    )

    summary["freq_cell"] = summary["time"].map(time_to_cell)
    summary = summary[(summary["freq_mol"] > 0) & (summary["freq_cell"] > 0)].copy()

    if len(summary) < 2:
        raise RuntimeError("Need at least 2 positive calibration points to fit the calibration curve.")

    return summary


# =========================
# Plot helpers
# =========================
def make_figure(fig_width: float, fig_height: float, dpi: int):
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    return fig, ax


def place_legend_below(ax, legend_ncol: int, legend_yoffset: float):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h

    ax.legend(
        uniq.values(),
        uniq.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, legend_yoffset),
        ncol=legend_ncol,
        frameon=True,
        borderaxespad=0.0,
        handlelength=2.5,
    )


def finalize_figure(fig, out_plot: Path, dpi: int, no_show: bool, rect=(0, 0.14, 1, 1)):
    fig.tight_layout(rect=rect)
    fig.savefig(out_plot, dpi=dpi, bbox_inches="tight")
    if not no_show:
        plt.show()
    plt.close(fig)


def short_trajectory_title(subject_filter, freq_min, freq_max, n_non_target):
    parts = [f"Target + {n_non_target} representative non-target clonotypes"]
    if subject_filter is not None:
        parts.append(f"subject={subject_filter}")
    return " | ".join(parts)


def save_calibration_outputs(
    summary_df: pd.DataFrame,
    out_points: Path,
    out_fit: Path,
    out_curve: Path,
    out_plot: Path,
    fig_width: float,
    fig_height: float,
    dpi: int,
    point_color: str,
    line_color: str,
    title: str,
    no_show: bool,
    legend_ncol: int,
    legend_yoffset: float,
):
    a, b, r2 = fit_loglog_calibration(
        freq_mol=summary_df["freq_mol"].to_numpy(),
        freq_cell=summary_df["freq_cell"].to_numpy(),
    )

    fit_df = pd.DataFrame([{
        "model": "log10(freq_cell) = a + b * log10(freq_mol)",
        "intercept_a": a,
        "slope_b": b,
        "r2_log10": r2,
        "n_points": len(summary_df),
    }])

    x_curve = np.geomspace(
        summary_df["freq_mol"].min() * 0.8,
        summary_df["freq_mol"].max() * 1.2,
        300
    )
    y_curve = calibration_function(x_curve, a=a, b=b)

    curve_df = pd.DataFrame({
        "freq_mol_grid": x_curve,
        "freq_cell_fit": y_curve,
    })

    summary_df.to_csv(out_points, index=False)
    fit_df.to_csv(out_fit, index=False)
    curve_df.to_csv(out_curve, index=False)

    fig, ax = make_figure(fig_width, fig_height, dpi)

    ax.scatter(
        summary_df["freq_mol"],
        summary_df["freq_cell"],
        color=point_color,
        s=35,
        label="Target calibration points"
    )
    ax.plot(
        x_curve,
        y_curve,
        color=line_color,
        linewidth=2.0,
        label=f"Fit: log10(y) = {a:.3f} + {b:.3f} log10(x)\n$R^2$ = {r2:.3f}"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Target molecular frequency (freq_geo_x)", labelpad=10)
    ax.set_ylabel("Target cellular frequency", labelpad=10)
    ax.set_title(title, pad=10)

    place_legend_below(ax, legend_ncol=legend_ncol, legend_yoffset=legend_yoffset)
    finalize_figure(fig, out_plot, dpi=dpi, no_show=no_show)

    return fit_df, curve_df


# =========================
# Main
# =========================
def main():
    parser = argparse.ArgumentParser(
        description="Plot target clonotype trajectory and estimate two calibration curves: all points and observable-only."
    )
    parser.add_argument(
        "--input", "-i",
        nargs="+",
        required=True,
        help="One or more input CSV files and/or directories containing per_clone_denoised_subject.csv"
    )
    parser.add_argument("--outdir", "-o", required=True, help="Output directory")
    parser.add_argument("--target", default=DEFAULT_TARGET, help="Target aaSeqCDR3 sequence")
    parser.add_argument(
        "--target-cell-freqs",
        required=True,
        help="Comma-separated target cellular frequencies, one per sorted time point"
    )
    parser.add_argument("--subject", default=None, help="Optional subject/donor filter")
    parser.add_argument("--n-non-target", type=int, default=30, help="Number of non-target clonotypes")
    parser.add_argument("--freq-min", type=float, default=None, help="Minimum mean linear frequency for non-target selection")
    parser.add_argument("--freq-max", type=float, default=None, help="Maximum mean linear frequency for non-target selection")
    parser.add_argument("--n-freq-bins", type=int, default=5, help="Number of log-frequency bins for stratified non-target selection")

    # aesthetics
    parser.add_argument("--fig-width", type=float, default=7.0, help="Figure width in inches")
    parser.add_argument("--fig-height", type=float, default=5.0, help="Figure height in inches")
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    parser.add_argument("--yscale", choices=["linear", "log"], default="linear", help="Y scale for trajectory plot")
    parser.add_argument("--xmin", type=float, default=None, help="Minimum x for trajectory plot")
    parser.add_argument("--xmax", type=float, default=None, help="Maximum x for trajectory plot")
    parser.add_argument("--ymin", type=float, default=None, help="Minimum y for trajectory plot")
    parser.add_argument("--ymax", type=float, default=None, help="Maximum y for trajectory plot")
    parser.add_argument("--target-color", default="black", help="Color for target trajectory")
    parser.add_argument("--background-color", default="gray", help="Color for non-target trajectories")
    parser.add_argument("--alpha-non-target", type=float, default=0.25, help="Alpha for non-target trajectories")
    parser.add_argument("--linewidth-target", type=float, default=2.5, help="Line width for target")
    parser.add_argument("--linewidth-non-target", type=float, default=1.0, help="Line width for non-target")
    parser.add_argument("--marker-size", type=float, default=4.0, help="Marker size")
    parser.add_argument("--calib-color-all", default="tab:red", help="Color for all-points calibration fit")
    parser.add_argument("--calib-color-obs", default="tab:green", help="Color for observable-only calibration fit")
    parser.add_argument("--calib-point-color-all", default="tab:blue", help="Color for all-points calibration points")
    parser.add_argument("--calib-point-color-obs", default="tab:orange", help="Color for observable-only calibration points")
    parser.add_argument("--legend-ncol", type=int, default=1, help="Number of columns in legends")
    parser.add_argument("--legend-yoffset", type=float, default=-0.20, help="Vertical offset for legends placed below the x-axis")
    parser.add_argument("--no-show", action="store_true", help="Do not open interactive plot windows")

    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    input_files = collect_input_files(args.input)
    subject_filter = parse_subject_filter(args.subject)
    target_norm = args.target.strip().upper()
    target_cell_freqs = parse_target_cell_freqs(args.target_cell_freqs)

    print("Input files:")
    for f in input_files:
        print(" ", f)

    all_rows = []
    selection_rows = []

    for f in input_files:
        df = pd.read_csv(f, sep=None, engine="python")
        df.columns = [c.strip() for c in df.columns]

        if "aaSeqCDR3" not in df.columns:
            raise RuntimeError(f"Missing 'aaSeqCDR3' in {f}")

        time_col = infer_time_column(df, f)
        subject_col = infer_subject_column(df)
        observable_col = infer_observable_column(df)

        df["aaSeqCDR3_norm"] = df["aaSeqCDR3"].astype(str).str.strip().str.upper()
        df[time_col] = pd.to_numeric(df[time_col], errors="coerce")

        if df[time_col].isna().any():
            bad = df[df[time_col].isna()].head(5)
            raise RuntimeError(f"Non-numeric time values found in {f}. Example rows:\n{bad}")

        df["frequency"] = get_frequency_column(df, f)

        if observable_col is not None:
            df["observable_norm"] = as_bool_series(df[observable_col])
        else:
            df["observable_norm"] = df["frequency"].fillna(0.0) > 0

        if subject_filter is not None:
            if subject_col is None:
                raise RuntimeError(
                    f"--subject was provided, but no subject column was found in {f}\n"
                    f"Available columns: {list(df.columns)}"
                )
            keep = df[subject_col].astype(str).str.strip() == subject_filter
            df = df[keep].copy()

        if df.empty:
            print(f"\n[WARNING] File {f.name}: no rows left after filtering. Skipping.")
            continue

        non_target = select_non_target_clones_stratified(
            df=df,
            target_norm=target_norm,
            n_non_target=args.n_non_target,
            freq_min=args.freq_min,
            freq_max=args.freq_max,
            time_col=time_col,
            n_bins=args.n_freq_bins,
        )

        selected_norm: List[str] = []
        if (df["aaSeqCDR3_norm"] == target_norm).any():
            selected_norm.append(target_norm)
        else:
            print(f"\n[WARNING] Target clonotype not found in {f.name}")

        selected_norm.extend(non_target)

        if not selected_norm:
            print(f"\n[WARNING] No clonotypes selected in {f.name}. Skipping.")
            continue

        print(f"\n[file={f.name}] selected non-target clonotypes across frequency bins: {len(non_target)}")

        for c in selected_norm:
            selection_rows.append({
                "source_file": f.name,
                "subject": subject_filter if subject_filter is not None else "",
                "aaSeqCDR3_norm": c,
                "is_target": (c == target_norm),
            })

        df_sel = df[df["aaSeqCDR3_norm"].isin(selected_norm)].copy()
        df_sel["source_file"] = f.name
        df_sel["time"] = df_sel[time_col]

        keep_cols = []
        if subject_col is not None:
            keep_cols.append(subject_col)
        keep_cols += [
            "source_file",
            "aaSeqCDR3",
            "aaSeqCDR3_norm",
            "time",
            "frequency",
            "observable_norm",
        ]
        for extra in ["p_emp_low", "p_emp", "posterior_mean", "f_denoised", "freq_geo_x", "log_freq", "observable"]:
            if extra in df_sel.columns and extra not in keep_cols:
                keep_cols.append(extra)

        all_rows.append(df_sel[keep_cols].copy())

    if not all_rows:
        raise RuntimeError("No data extracted from any input file.")

    traj = pd.concat(all_rows, ignore_index=True).sort_values(["source_file", "aaSeqCDR3_norm", "time"])

    sel_df = (
        pd.DataFrame(selection_rows)
        .drop_duplicates()
        .sort_values(["source_file", "is_target"], ascending=[True, False])
    )

    obs_table = (
        traj.pivot_table(
            index=["source_file", "aaSeqCDR3_norm", "time"],
            values="observable_norm",
            aggfunc="max",
        )
        .fillna(False)
        .astype(bool)
    )

    target_df = traj[traj["aaSeqCDR3_norm"] == target_norm].copy()
    if target_df.empty:
        raise RuntimeError("Target clonotype was not found in the extracted trajectories, so calibration cannot be computed.")

    calib_all = build_calibration_summary(target_df=target_df, target_cell_freqs=target_cell_freqs, observable_only=False)
    calib_obs = build_calibration_summary(target_df=target_df, target_cell_freqs=target_cell_freqs, observable_only=True)

    p_sel = outdir / "selected_clones.csv"
    p_traj = outdir / "extracted_trajectories.csv"
    p_obs = outdir / "observability_matrix.csv"
    p_png_traj = outdir / "target_plus_background.png"

    p_calib_pts_all = outdir / "target_calibration_points_all.csv"
    p_calib_fit_all = outdir / "target_calibration_fit_all.csv"
    p_calib_curve_all = outdir / "target_calibration_curve_all.csv"
    p_png_calib_all = outdir / "target_calibration_all_points.png"

    p_calib_pts_obs = outdir / "target_calibration_points_observable.csv"
    p_calib_fit_obs = outdir / "target_calibration_fit_observable.csv"
    p_calib_curve_obs = outdir / "target_calibration_curve_observable.csv"
    p_png_calib_obs = outdir / "target_calibration_observable_only.png"

    sel_df.to_csv(p_sel, index=False)
    traj.to_csv(p_traj, index=False)
    obs_table.to_csv(p_obs, index=True)

    # Plot 1: trajectory
    fig, ax = make_figure(args.fig_width, args.fig_height, args.dpi)

    for (_, clone), g in traj.groupby(["source_file", "aaSeqCDR3_norm"]):
        is_target = clone == target_norm
        ax.plot(
            g["time"],
            g["frequency"],
            marker="o",
            markersize=args.marker_size,
            linewidth=args.linewidth_target if is_target else args.linewidth_non_target,
            alpha=1.0 if is_target else args.alpha_non_target,
            color=args.target_color if is_target else args.background_color,
            label=("TARGET: CASSFSTCSANYGYTF" if is_target else None),
        )

    ax.set_xlabel("Time point", labelpad=10)
    ax.set_ylabel("Clonotype molecular frequency", labelpad=10)
    ax.set_yscale(args.yscale)

    if args.xmin is not None or args.xmax is not None:
        ax.set_xlim(left=args.xmin, right=args.xmax)
    if args.ymin is not None or args.ymax is not None:
        ax.set_ylim(bottom=args.ymin, top=args.ymax)

    ax.set_title(short_trajectory_title(subject_filter, args.freq_min, args.freq_max, args.n_non_target), pad=10)

    place_legend_below(ax, legend_ncol=args.legend_ncol, legend_yoffset=args.legend_yoffset)
    finalize_figure(fig, p_png_traj, dpi=args.dpi, no_show=args.no_show)

    # Plot 2: calibration all points
    fit_all_df, curve_all_df = save_calibration_outputs(
        summary_df=calib_all,
        out_points=p_calib_pts_all,
        out_fit=p_calib_fit_all,
        out_curve=p_calib_curve_all,
        out_plot=p_png_calib_all,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
        dpi=args.dpi,
        point_color=args.calib_point_color_all,
        line_color=args.calib_color_all,
        title="Spike-in calibration (all target time points)",
        no_show=args.no_show,
        legend_ncol=args.legend_ncol,
        legend_yoffset=args.legend_yoffset,
    )

    # Plot 3: calibration observable only
    fit_obs_df, curve_obs_df = save_calibration_outputs(
        summary_df=calib_obs,
        out_points=p_calib_pts_obs,
        out_fit=p_calib_fit_obs,
        out_curve=p_calib_curve_obs,
        out_plot=p_png_calib_obs,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
        dpi=args.dpi,
        point_color=args.calib_point_color_obs,
        line_color=args.calib_color_obs,
        title="Spike-in calibration (target points with observable = TRUE)",
        no_show=args.no_show,
        legend_ncol=args.legend_ncol,
        legend_yoffset=args.legend_yoffset,
    )

    print("\nSaved outputs:")
    print(f" - {p_sel}")
    print(f" - {p_traj}")
    print(f" - {p_obs}")
    print(f" - {p_png_traj}")
    print(f" - {p_calib_pts_all}")
    print(f" - {p_calib_fit_all}")
    print(f" - {p_calib_curve_all}")
    print(f" - {p_png_calib_all}")
    print(f" - {p_calib_pts_obs}")
    print(f" - {p_calib_fit_obs}")
    print(f" - {p_calib_curve_obs}")
    print(f" - {p_png_calib_obs}")

    print("\nCalibration summary (all points):")
    print(fit_all_df.to_string(index=False))

    print("\nCalibration summary (observable-only):")
    print(fit_obs_df.to_string(index=False))


if __name__ == "__main__":
    main()