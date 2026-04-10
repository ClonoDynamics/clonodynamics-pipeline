#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
build_multipanel.py

Independent multipanel builder from existing single-plot images.

Now supports:
- raster inputs: PNG, JPG, JPEG, TIFF, BMP, WEBP
- PDF inputs: first page is rasterized and used as a panel
- mixed inputs: e.g. panel1.pdf, panel2.png, panel3.jpg

Goal: preserve "original size" of singles in the multipanel as much as possible.

Definitions of "original":
  1) Pixel-identical: do NOT resample raster images (no resizing).
  2) Inch-identical: the image occupies the same physical size (inches) as when
     the single was saved.

Modes:
  --mode preserve_inches   (DEFAULT)
      - keeps raster pixels intact (no resampling)
      - keeps inches intact by using the SINGLE DPI as the output DPI
      - requires knowing the single DPI:
          * tries to read embedded DPI metadata for raster images
          * if missing, uses --fallback_single_dpi
      - for PDF inputs, the page is rasterized at --pdf_dpi
      - if --out_dpi differs from the inferred single DPI, we FORCE out_dpi
        to match the reference DPI in preserve_inches mode

  --mode preserve_pixels
      - keeps raster pixels intact
      - ignores physical inches (depends on --out_dpi)
      - for PDF inputs, rasterization still occurs at --pdf_dpi

Inputs:
  Use either:
    --files "a.png,b.pdf,c.jpg"
  or:
    --array "msd,median,ppos" --pattern "{name}.png" --in_dir <dir>

Layout:
  --layout 2x2
  --fill row|col controls placement order

Labels:
  auto: --labels upper|lower|numeric|roman|none
  or: --label_list "A,B,C,D"
  formatting: --label_style "({})"  (must include {})

Output:
  Provide --out as a FILE path (png or pdf recommended).
  Or provide --out_dir + --out_name.

Examples
--------
python build_multipanel.py \
    --layout 1x2 --fill col \
    --files "./fig1.png,./fig2.pdf" \
    --mode preserve_inches \
    --single_dpi auto \
    --fallback_single_dpi 200 \
    --pdf_dpi 300 \
    --out ./multipanel.png \
    --labels upper \
    --label_style "{})" \
    --font Arial \
    --font_size 12

    Fig. 2
    python ./code/utils/buil_multipanel.py \
    --mode preserve_inches \
    --single_dpi auto \
    --pdf_dpi 300 \
    --labels upper \
    --label_style "{})" \
    --font Arial \
    --font_size 12 \
    --out ./figures/main/figure_2_prova.png \
    --files "./figures/5-compare_conditions_dt1/fig_2A-C/compare_median_dx_dt1.png, \
    ./figures/5-compare_conditions_dt1/fig_2A-C/compare_ppos_dt1.png, \
    ./figures/5-compare_conditions_dt1/fig_2A-C/compare_msd_dt1.png, \
    ./figures/5-finite_time_diagnostics/prova3/figA_drift_overlay.png, \
    ./figures/5-finite_time_diagnostics/prova3/ppos.png, \
    ./figures/5-finite_time_diagnostics/prova3/figB_diffusion_overlay.png, \
    ./figures/5-finite_time_diagnostics/prova3/figC_xstar_vs_dt.png, \
    ./figures/5-finite_time_diagnostics/prova3/figE_slope_vs_dt.png, \
    ./figures/5-finite_time_diagnostics/prova3/figG_D_pooled_vs_dt.png" \
    --layout 3x3 --fill row

    python ./code/utils/buil_multipanel.py \
    --mode preserve_inches \
    --single_dpi auto \
    --pdf_dpi 300 \
    --labels upper \
    --label_style "{})" \
    --font Arial \
    --font_size 12 \
    --out ./figures/supplementary/figure_SX.png \
    --files "./figures/5-compare_conditions_dt1/figure_SX.png/compare_median_dx_dt1.png, \
    ./figures/5-compare_conditions_dt1/figure_SX.png/compare_ppos_dt1.png, \
    ./figures/5-compare_conditions_dt1/figure_SX.png/compare_msd_dt1.png, \
    ./figures/5-compare_conditions_dt1/figure_SY.png/compare_xstar_vs_dt.png, \
    ./figures/5-compare_conditions_dt1/figure_SY.png/compare_slope_vs_dt.png, \
    ./figures/5-compare_conditions_dt1/figure_SY.png/compare_Dmad_pooled_vs_dt.png" \
    --layout 2x3 --fill row

    
    
    Fig.3
python ./code/utils/buil_multipanel.py \
    --mode preserve_inches \
    --single_dpi auto \
    --fallback_single_dpi 200 \
    --pdf_dpi 300 \
    --out ./figures/main/prova.png \
    --labels upper \
    --label_style "{})" \
    --font Arial \
    --font_size 12 \
--files "./figures/8-fp_stationary/TT_p05_w_clone/fp.grid.p_emp_vs_p.png, \
./figures/8-fp_stationary/TT_p05_w_clone/fp.grid.ratio.png, \
./figures/8-fp_stationary/TT_p05_w_clone/fp.grid.diffusion.png, \
./figures/8-fp_stationary/TT_p02_w_clone/fp.grid.p_emp_vs_p.png, \
./figures/8-fp_stationary/TT_p02_w_clone/fp.grid.ratio.png, \
./figures/8-fp_stationary/TT_p02_w_clone/fp.grid.diffusion.png, \
./figures/8-fp_stationary/TT_p005_w_clone/fp.grid.p_emp_vs_p.png, \
./figures/8-fp_stationary/TT_p005_w_clone/fp.grid.ratio.png, \
./figures/8-fp_stationary/TT_p005_w_clone/fp.grid.diffusion.png" \
--layout 3x3  --fill col


Notes
-----
- Images are embedded without resizing. If cells are bigger than an image,
  the image is centered and the remaining area is blank.
- If images in the same column/row have different sizes, column width/row height
  becomes the max in that column/row.
- PDF inputs are rasterized from page 1 at --pdf_dpi before placement.

Figure 1
python ./code/utils/buil_multipanel.py \
--files "./figures/1-powerlaw_qc/figA_ccdf_pooled.png,./figures/1-powerlaw_qc/figB_model_compare_hist.png,./figures/1-powerlaw_qc/figC_threshold_stability.png,./figures/1-powerlaw_qc/figD_gamma_scatter.png,./figures/2-diagnostics_qc/figT_delta_gamma_subject_time.png,./figures/2-diagnostics_qc/figW_qc_map.png,./figures/2-diagnostics_qc/figS_gamma_by_subject.png,./figures/2-diagnostics_qc/figP_xmin_by_subject.png,./figures/2-diagnostics_qc/figQ_ks_by_subject.png" \
    --mode preserve_inches \
    --single_dpi auto \
    --fallback_single_dpi 200 \
    --pdf_dpi 300 \
    --out ./figures/main/figure_1b.pdf \
    --labels upper \
    --label_style "{})" \
    --font Arial \
    --font_size 12 --layout 3x3
python ./code/utils/buil_multipanel.py \
--files "./figures/9-fp_boot_full_pipeline/p005_w_clone/fig_drift_bootstrap.png, \ 
./figures/9-fp_boot_full_pipeline/p005_w_clone/fig_diffusion_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p005_w_clone/fig_flux_J_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p005_w_clone/fig_source_S_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p005_w_clone/fig_tau_distribution.png, \
./figures/9-fp_boot_full_pipeline/p02_w_clone/fig_drift_bootstrap.png, \ 
./figures/9-fp_boot_full_pipeline/p02_w_clone/fig_diffusion_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p02_w_clone/fig_flux_J_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p02_w_clone/fig_source_S_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p02_w_clone/fig_tau_distribution.png, \
./figures/9-fp_boot_full_pipeline/p05_w_clone/fig_drift_bootstrap.png, \ 
./figures/9-fp_boot_full_pipeline/p05_w_clone/fig_diffusion_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p05_w_clone/fig_flux_J_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p05_w_clone/fig_source_S_bootstrap.png, \
./figures/9-fp_boot_full_pipeline/p05_w_clone/fig_tau_distribution.png" \
--mode preserve_inches \
--single_dpi auto \
--fallback_single_dpi 200 \
--pdf_dpi 300 \
--out ./figures/main/prova.png \
--labels upper \
--label_style "{})" \
--font Arial \
--font_size 12 --layout 5x3 --fill row
    
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

try:
    from pdf2image import convert_from_path
except ImportError:
    convert_from_path = None


# -------------------------
# Parsing helpers
# -------------------------
def parse_layout(s: str) -> Tuple[int, int]:
    s = s.strip().lower().replace("×", "x")
    if "x" not in s:
        raise ValueError(f"Bad --layout '{s}'. Use like '2x2'.")
    r, c = s.split("x", 1)
    return int(r), int(c)


def parse_list_str(s: str) -> List[str]:
    s = (s or "").strip()
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def parse_4(s: str) -> Tuple[float, float, float, float]:
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 4:
        raise ValueError(f"Expected 'l,r,t,b' but got '{s}'")
    return float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])


def excel_letters(n: int, upper: bool = True) -> str:
    if n <= 0:
        return str(n)
    chars = []
    x = n
    base = 65 if upper else 97
    while x > 0:
        x, rem = divmod(x - 1, 26)
        chars.append(chr(base + rem))
    return "".join(reversed(chars))


def roman(n: int) -> str:
    if n <= 0:
        return str(n)
    vals = [
        (1000, "m"), (900, "cm"), (500, "d"), (400, "cd"),
        (100, "c"), (90, "xc"), (50, "l"), (40, "xl"),
        (10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")
    ]
    out = []
    x = n
    for v, sym in vals:
        while x >= v:
            out.append(sym)
            x -= v
    return "".join(out)


def make_labels(kind: str, n: int, start: int) -> List[str]:
    labels = []
    for i in range(n):
        k = start + i
        if kind == "none":
            labels.append("")
        elif kind == "numeric":
            labels.append(str(k))
        elif kind == "roman":
            labels.append(roman(k))
        elif kind == "upper":
            labels.append(excel_letters(k, upper=True))
        elif kind == "lower":
            labels.append(excel_letters(k, upper=False))
        else:
            raise ValueError(kind)
    return labels


def label_xy(anchor: str, dx: float, dy: float, pos: str) -> Tuple[float, float, str, str]:
    if anchor == "tl":
        x0, y0, ha, va, sx, sy = 0.0, 1.0, "left", "top", +1, -1
    elif anchor == "tr":
        x0, y0, ha, va, sx, sy = 1.0, 1.0, "right", "top", -1, -1
    elif anchor == "bl":
        x0, y0, ha, va, sx, sy = 0.0, 0.0, "left", "bottom", +1, +1
    elif anchor == "br":
        x0, y0, ha, va, sx, sy = 1.0, 0.0, "right", "bottom", -1, +1
    else:
        raise ValueError(anchor)

    if pos == "inside":
        x = x0 + sx * dx
        y = y0 + sy * dy
    else:
        x = x0 - sx * dx
        y = y0 - sy * dy
    return x, y, ha, va


# -------------------------
# Input / image helpers
# -------------------------
RASTER_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
PDF_SUFFIXES = {".pdf"}


def is_pdf(path: Path) -> bool:
    return path.suffix.lower() in PDF_SUFFIXES


def is_raster(path: Path) -> bool:
    return path.suffix.lower() in RASTER_SUFFIXES


def require_pdf_support() -> None:
    if convert_from_path is None:
        raise ImportError(
            "PDF input requested, but pdf2image is not installed.\n"
            "Install with: pip install pdf2image\n"
            "Also install poppler on your system."
        )


def render_pdf_first_page(path: Path, pdf_dpi: int) -> Image.Image:
    require_pdf_support()
    pages = convert_from_path(str(path), dpi=pdf_dpi, first_page=1, last_page=1)
    if not pages:
        raise ValueError(f"Could not rasterize PDF: {path}")
    return pages[0].convert("RGB")


def read_image_rgb(path: Path, pdf_dpi: int) -> np.ndarray:
    suffix = path.suffix.lower()

    if suffix in PDF_SUFFIXES:
        return np.asarray(render_pdf_first_page(path, pdf_dpi))

    if suffix in RASTER_SUFFIXES:
        return np.asarray(Image.open(path).convert("RGB"))

    raise ValueError(f"Unsupported input format: {path}")


def image_px(path: Path, pdf_dpi: int) -> Tuple[int, int]:
    suffix = path.suffix.lower()

    if suffix in PDF_SUFFIXES:
        im = render_pdf_first_page(path, pdf_dpi)
        w, h = im.size
        return int(w), int(h)

    if suffix in RASTER_SUFFIXES:
        with Image.open(path) as im:
            w, h = im.size
        return int(w), int(h)

    raise ValueError(f"Unsupported input format: {path}")


def infer_dpi(path: Path, fallback_single_dpi: float, pdf_dpi: int) -> float:
    """
    Infer DPI for preserve_inches logic.
    - Raster images: read embedded DPI if present, else fallback.
    - PDF: use pdf_dpi as effective DPI because the page is rasterized at that DPI.
    """
    suffix = path.suffix.lower()

    if suffix in PDF_SUFFIXES:
        return float(pdf_dpi)

    try:
        with Image.open(path) as im:
            info = getattr(im, "info", {}) or {}
            if "dpi" in info and isinstance(info["dpi"], tuple) and len(info["dpi"]) >= 2:
                xdpi, ydpi = info["dpi"][0], info["dpi"][1]
                if xdpi and ydpi and xdpi > 1 and ydpi > 1:
                    return float(xdpi)
    except Exception:
        pass

    return float(fallback_single_dpi)


def compute_grid_cell_sizes_px(
    files: List[Path],
    nrows: int,
    ncols: int,
    fill: str,
    pdf_dpi: int
) -> Tuple[List[int], List[int], List[Tuple[int, int]]]:
    def idx_to_rc(k: int) -> Tuple[int, int]:
        if fill == "row":
            return k // ncols, k % ncols
        return k % nrows, k // nrows

    placements = [idx_to_rc(k) for k in range(len(files))]
    col_w = [0] * ncols
    row_h = [0] * nrows

    for k, p in enumerate(files):
        r, c = placements[k]
        w, h = image_px(p, pdf_dpi=pdf_dpi)
        col_w[c] = max(col_w[c], w)
        row_h[r] = max(row_h[r], h)

    return col_w, row_h, placements


def total_canvas_px(
    col_w: List[int],
    row_h: List[int],
    wspace_px: int,
    hspace_px: int,
    margins_px: Tuple[int, int, int, int]
) -> Tuple[int, int]:
    left, right, top, bottom = margins_px
    W = left + sum(col_w) + right + (len(col_w) - 1) * wspace_px
    H = bottom + sum(row_h) + top + (len(row_h) - 1) * hspace_px
    return int(W), int(H)


# -------------------------
# CLI
# -------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Build a multipanel from existing images and PDFs.")

    ap.add_argument("--in_dir", default=".", help="Base directory for resolving --array/--pattern or relative --files.")

    ap.add_argument(
        "--files",
        default="",
        help="Comma-separated list of input file paths. Supports png,jpg,jpeg,tif,tiff,bmp,webp,pdf."
    )
    ap.add_argument(
        "--array",
        default="",
        help="Comma-separated panel names in order (used if --files is empty)."
    )
    ap.add_argument(
        "--pattern",
        default="{name}.png",
        help="Pattern to resolve --array names under --in_dir. Change this if using PDFs, e.g. '{name}.pdf'."
    )

    ap.add_argument("--layout", default="2x2", help="Grid layout like '2x2'.")
    ap.add_argument("--fill", choices=["row", "col"], default="row", help="Placement order.")
    ap.add_argument("--empty", choices=["off", "keep"], default="off", help="If grid has extra slots.")
    ap.add_argument("--bg", default="white", help="Background color.")

    ap.add_argument(
        "--mode",
        choices=["preserve_inches", "preserve_pixels"],
        default="preserve_inches",
        help="Preserve inches or preserve pixels."
    )

    ap.add_argument(
        "--single_dpi",
        default="auto",
        help="Single-plot DPI. 'auto' reads raster metadata when possible; PDFs use --pdf_dpi."
    )
    ap.add_argument(
        "--fallback_single_dpi",
        type=float,
        default=200.0,
        help="Fallback DPI when metadata is missing."
    )
    ap.add_argument(
        "--out_dpi",
        type=float,
        default=200.0,
        help="Output DPI. In preserve_inches mode may be forced to match reference DPI."
    )
    ap.add_argument(
        "--pdf_dpi",
        type=int,
        default=300,
        help="DPI used to rasterize PDF inputs."
    )

    ap.add_argument("--pad_x_in", type=float, default=0.15, help="Horizontal gap between panels in inches.")
    ap.add_argument("--pad_y_in", type=float, default=0.18, help="Vertical gap between panels in inches.")
    ap.add_argument(
        "--margins_in",
        default="0.2,0.2,0.2,0.2",
        help="Margins (left,right,top,bottom) in inches."
    )

    ap.add_argument("--labels", choices=["none", "upper", "lower", "numeric", "roman"], default="upper")
    ap.add_argument("--label_start", type=int, default=1)
    ap.add_argument("--label_list", default="", help="Explicit comma-separated labels overriding --labels.")
    ap.add_argument("--label_style", default="({})", help="Label formatter, must include {}.")
    ap.add_argument("--label_pos", choices=["inside", "outside"], default="inside")
    ap.add_argument("--label_anchor", choices=["tl", "tr", "bl", "br"], default="tl")
    ap.add_argument("--label_dx", type=float, default=0.02)
    ap.add_argument("--label_dy", type=float, default=0.02)
    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font_size", type=float, default=9.0)

    ap.add_argument("--out", default="", help="Full output file path (.png or .pdf).")
    ap.add_argument("--out_dir", default="", help="Output directory if --out is not used.")
    ap.add_argument("--out_name", default="multipanel", help="Base output filename if --out is not used.")
    ap.add_argument("--fmt", choices=["png", "pdf"], default="png", help="Output format if --out is not used.")
    ap.add_argument("--prefix", default="")
    ap.add_argument("--suffix", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--auto_tag", action="store_true")

    return ap


# -------------------------
# Main
# -------------------------
def main() -> None:
    args = build_argparser().parse_args()

    in_dir = Path(args.in_dir).expanduser().resolve()

    files: List[Path] = []
    if args.files.strip():
        for item in parse_list_str(args.files):
            p = Path(item).expanduser()
            if not p.is_absolute():
                p = in_dir / p
            p = p.resolve()
            if not p.exists():
                raise FileNotFoundError(f"Missing file: {p}")
            files.append(p)
    else:
        names = parse_list_str(args.array)
        if not names:
            raise ValueError("Provide either --files or --array.")
        for nm in names:
            rel = args.pattern.format(name=nm)
            p = (in_dir / rel).resolve()
            if not p.exists():
                raise FileNotFoundError(f"Missing input for '{nm}': {p}")
            files.append(p)

    n_panels = len(files)
    nrows, ncols = parse_layout(args.layout)

    if nrows * ncols < n_panels:
        raise ValueError(f"Layout {args.layout} has {nrows*ncols} slots but you have {n_panels} panels.")

    # Validate input types early
    for p in files:
        suff = p.suffix.lower()
        if suff not in RASTER_SUFFIXES and suff not in PDF_SUFFIXES:
            raise ValueError(f"Unsupported input format: {p}")

    # DPI handling
    if args.single_dpi.strip().lower() == "auto":
        panel_dpis = [
            infer_dpi(p, fallback_single_dpi=args.fallback_single_dpi, pdf_dpi=args.pdf_dpi)
            for p in files
        ]
    else:
        forced = float(args.single_dpi)
        panel_dpis = [forced] * n_panels

    single_dpi_ref = float(np.median(panel_dpis)) if panel_dpis else float(args.fallback_single_dpi)

    out_dpi = float(args.out_dpi)
    if args.mode == "preserve_inches":
        if abs(out_dpi - single_dpi_ref) > 1e-6:
            out_dpi = single_dpi_ref

    col_w_px, row_h_px, placements = compute_grid_cell_sizes_px(
        files=files,
        nrows=nrows,
        ncols=ncols,
        fill=args.fill,
        pdf_dpi=args.pdf_dpi,
    )

    margins_in = parse_4(args.margins_in)
    margin_px = (
        int(round(margins_in[0] * out_dpi)),
        int(round(margins_in[1] * out_dpi)),
        int(round(margins_in[2] * out_dpi)),
        int(round(margins_in[3] * out_dpi)),
    )
    wspace_px = int(round(float(args.pad_x_in) * out_dpi))
    hspace_px = int(round(float(args.pad_y_in) * out_dpi))

    canvas_w_px, canvas_h_px = total_canvas_px(col_w_px, row_h_px, wspace_px, hspace_px, margin_px)
    figsize = (canvas_w_px / out_dpi, canvas_h_px / out_dpi)

    if args.label_list.strip():
        labels = parse_list_str(args.label_list)
        if len(labels) != n_panels:
            raise ValueError(f"--label_list has {len(labels)} labels but there are {n_panels} panels.")
    else:
        labels = make_labels(args.labels, n_panels, int(args.label_start))

    if args.out.strip():
        out_path = Path(args.out).expanduser()
    else:
        out_dir = Path(args.out_dir).expanduser() if args.out_dir.strip() else in_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        auto = f"_{in_dir.name}" if args.auto_tag else ""
        tag = f"_{args.tag}" if args.tag else ""
        final_tag = f"{auto}{tag}"
        out_name = f"{args.prefix}{args.out_name}{final_tag}{args.suffix}.{args.fmt}"
        out_path = out_dir / out_name

    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.family"] = args.font
    plt.rcParams["font.size"] = float(args.font_size)

    fig = plt.figure(figsize=figsize, facecolor=args.bg, dpi=out_dpi)

    left_px, right_px, top_px, bottom_px = margin_px
    left = left_px / canvas_w_px
    right = 1.0 - (right_px / canvas_w_px)
    bottom = bottom_px / canvas_h_px
    top = 1.0 - (top_px / canvas_h_px)

    mean_col_px = max(1.0, float(np.mean([x for x in col_w_px if x > 0] or [1.0])))
    mean_row_px = max(1.0, float(np.mean([x for x in row_h_px if x > 0] or [1.0])))
    wspace = wspace_px / mean_col_px
    hspace = hspace_px / mean_row_px

    gs = fig.add_gridspec(
        nrows=nrows,
        ncols=ncols,
        width_ratios=[max(1, x) for x in col_w_px],
        height_ratios=[max(1, x) for x in row_h_px],
        left=left, right=right, bottom=bottom, top=top,
        wspace=wspace, hspace=hspace
    )

    def idx_to_rc(k: int) -> Tuple[int, int]:
        if args.fill == "row":
            return k // ncols, k % ncols
        return k % nrows, k // ncols

    for k in range(nrows * ncols):
        if args.fill == "row":
            r, c = k // ncols, k % ncols
        else:
            r, c = k % nrows, k // nrows

        ax = fig.add_subplot(gs[r, c])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

        if k < n_panels:
            img = read_image_rgb(files[k], pdf_dpi=args.pdf_dpi)
            ih, iw = img.shape[0], img.shape[1]
            cell_w = col_w_px[c]
            cell_h = row_h_px[r]

            ax.set_xlim(0, cell_w)
            ax.set_ylim(cell_h, 0)

            x0 = (cell_w - iw) / 2.0
            y0 = (cell_h - ih) / 2.0
            ax.imshow(img, extent=[x0, x0 + iw, y0 + ih, y0], interpolation="none")

            lab = labels[k]
            if lab:
                txt = args.label_style.format(lab)
                x, y, ha, va = label_xy(
                    args.label_anchor,
                    float(args.label_dx),
                    float(args.label_dy),
                    args.label_pos
                )
                ax.text(
                    x, y, txt,
                    transform=ax.transAxes,
                    ha=ha, va=va,
                    fontsize=float(args.font_size),
                    fontfamily=args.font,
                    color="black",
                )
        else:
            if args.empty == "off":
                ax.axis("off")

    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path, dpi=out_dpi, facecolor=args.bg)
    else:
        fig.savefig(out_path, facecolor=args.bg)

    plt.close(fig)

    print(f"[SAVED] {out_path.resolve()}")
    print(f"[MODE] {args.mode}")
    print(f"[SINGLE_DPI_REF] {single_dpi_ref:.3f}")
    print(f"[OUT_DPI_USED] {out_dpi:.3f}")
    print(f"[PDF_DPI] {args.pdf_dpi}")
    print(f"[CANVAS_PX] {canvas_w_px} x {canvas_h_px}")
    print(f"[FIGSIZE_IN] {figsize[0]:.3f} x {figsize[1]:.3f}")


if __name__ == "__main__":
    main()