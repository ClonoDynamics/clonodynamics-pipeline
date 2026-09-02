#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clonodynamics_orchestrator_plotting.py
======================================

Interactive plotting orchestrator for ClonoDynamics.

The orchestrator separates visualization from scientific analysis, resolves
an exact Step output directory from either that directory or a parent results
root, and writes figures to <step-output>/figures.

Step 16 is coordinated with its plotter-specific CLI:

    python3 ./code/plotting/16-plot_observation_threshold_robustness.py \
        --input-dir <longitudinal>/16-observation_threshold_robustness \
        --figdir <longitudinal>/16-observation_threshold_robustness/figures \
        --formats pdf \
        --font-family Arial \
        --font-size 8 \
        --ci-mode reference

The P16 defaults above are applied automatically by this orchestrator. Extra
arguments forwarded after `--` are appended last and may override non-protected
plot settings.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


VERSION = "v7-plotting-orchestrator-step16-2026-09-01"


@dataclass(frozen=True)
class Plotter:
    key: str
    step: int
    title: str
    script: str
    analysis_dir_name: str
    input_option: str = "--analysis-dir"
    output_option: Optional[str] = "--figures-dir"
    default_args: Tuple[str, ...] = ()


PLOTTER_LIST: Tuple[Plotter, ...] = (
    Plotter(
        "step1", 1, "Repertoire characterization",
        "1-plot_repertoire_characterization.py",
        "1-repertoire_characterization",
    ),
    Plotter(
        "step3", 3, "Latent posterior-quality characterization",
        "3-plot_latent_posterior_quality_characterization.py",
        "3-latent_posterior_quality_characterization",
    ),
    Plotter(
        "step6", 6, "Latent-transition dataset characterization",
        "6-plot_latent_transition_dataset_characterization.py",
        "6-latent_transition_dataset_characterization",
    ),
    Plotter(
        "step7", 7, "Pseudo-longitudinal conditioning benchmark",
        "7-plot_pseudo_longitudinal_conditioning_benchmark.py",
        "7-pseudo_longitudinal_conditioning_benchmark",
    ),
    Plotter(
        "step8", 8, "Pseudo-longitudinal representation assessment",
        "8-plot_pseudo_longitudinal_representation_assessment.py",
        "8-pseudo_longitudinal_representation_assessment",
    ),
    Plotter(
        "step9", 9, "Replicate-decoupled forward drift",
        "9-plot_replicate_decoupled_forward_drift.py",
        "9-replicate_decoupled_forward_drift",
    ),
    Plotter(
        "step10", 10, "Longitudinal versus pseudo forward drift",
        "10-plot_longitudinal_vs_pseudo_forward_drift.py",
        "10-longitudinal_vs_pseudo_forward_drift",
    ),
    Plotter(
        "step11", 11, "Longitudinal fluctuation dynamics",
        "11-plot_longitudinal_fluctuation_dynamics.py",
        "11-longitudinal_fluctuation_dynamics",
    ),
    Plotter(
        "step12", 12, "Longitudinal versus pseudo fluctuation baseline",
        "12-plot_longitudinal_vs_pseudo_fluctuation_baseline.py",
        "12-longitudinal_vs_pseudo_fluctuation_baseline",
    ),
    Plotter(
        "step13", 13, "Temporal fluctuation scaling",
        "13-plot_temporal_fluctuation_scaling.py",
        "13-temporal_fluctuation_scaling",
    ),
    Plotter(
        "step14", 14, "Interval-position structure",
        "14-plot_interval_position_structure.py",
        "14-interval_position_structure",
    ),
    Plotter(
        "step15", 15, "Observation-boundary sensitivity",
        "15-plot_detectability_boundary_sensitivity.py",
        "15-detectability_boundary_sensitivity",
    ),
    Plotter(
        "step16", 16, "Operational observation-threshold robustness",
        "16-plot_observation_threshold_robustness.py",
        "16-observation_threshold_robustness",
        input_option="--input-dir",
        output_option="--figdir",
        default_args=(
            "--formats", "pdf",
            "--font-family", "Arial",
            "--font-size", "8",
            "--ci-mode", "reference",
        ),
    ),
)

PLOTTERS: Dict[str, Plotter] = {item.key: item for item in PLOTTER_LIST}
STEP_TO_KEY: Dict[int, str] = {item.step: item.key for item in PLOTTER_LIST}
MENU_TO_KEY: Dict[int, str] = {
    index: item.key for index, item in enumerate(PLOTTER_LIST, start=1)
}

PROTECTED_ARGS = {
    "--analysis-dir",
    "--benchmark-dir",
    "--input-dir",
    "--figures-dir",
    "--figdir",
    "--outdir",
}


# =============================================================================
# Generic utilities
# =============================================================================


def shell_join(command: Sequence[str]) -> str:
    try:
        return shlex.join([str(value) for value in command])
    except AttributeError:
        return " ".join(shlex.quote(str(value)) for value in command)


def strip_matching_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def dedupe_paths(paths: Sequence[Path]) -> List[Path]:
    output: List[Path] = []
    seen = set()
    for path in paths:
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path.absolute()
        token = str(resolved)
        if token not in seen:
            seen.add(token)
            output.append(resolved)
    return output


def nearest_existing_parent(path: Path) -> Optional[Path]:
    current = path
    while True:
        if current.exists() and current.is_dir():
            return current.resolve()
        if current.parent == current:
            return None
        current = current.parent


# =============================================================================
# Menu and selector normalization
# =============================================================================


def banner() -> None:
    print()
    print("=" * 80)
    print("ClonoDynamics | Figures")
    print("Plot already-generated ClonoDynamics analysis outputs")
    print("Plotting orchestrator | {}".format(VERSION))
    print("=" * 80)
    print()
    print("Choose a plotting module. Enter either the exact Step analysis directory")
    print("or a parent ClonoDynamics results directory. The Step directory will be")
    print("resolved automatically and figures will be written to its 'figures' folder.")
    print(flush=True)


def print_menu() -> None:
    print("Available plotters\n")
    for menu_number, plotter in enumerate(PLOTTER_LIST, start=1):
        print(" [{:>2}] {}".format(menu_number, plotter.script))
    print(" [ 0] Exit")
    print(flush=True)


def prompt_plot() -> Optional[str]:
    print_menu()
    maximum = len(PLOTTER_LIST)
    while True:
        raw = input("\nSelection [0-{}]: ".format(maximum)).strip()
        if raw.lower() in {"0", "q", "quit", "exit"}:
            return None
        try:
            choice = int(raw)
        except ValueError:
            print("Please enter a number from 0 to {}.".format(maximum), flush=True)
            continue
        if choice in MENU_TO_KEY:
            key = MENU_TO_KEY[choice]
            selected = PLOTTERS[key]
            print("\nSelected: {}".format(selected.script), flush=True)
            return key
        print("Please enter a number from 0 to {}.".format(maximum), flush=True)


def normalize_plot(value: str) -> str:
    raw = str(value).strip().lower().replace("_", "-")
    compact = raw.replace("-", "")
    if compact.startswith("step"):
        suffix = compact[4:]
        if suffix.isdigit() and int(suffix) in STEP_TO_KEY:
            return STEP_TO_KEY[int(suffix)]
    if raw.isdigit() and int(raw) in STEP_TO_KEY:
        return STEP_TO_KEY[int(raw)]

    aliases = {
        "repertoire": "step1",
        "posterior": "step3",
        "transition-dataset": "step6",
        "transition-characterization": "step6",
        "conditioning": "step7",
        "representation": "step8",
        "forward": "step9",
        "forward-drift": "step9",
        "forward-compare": "step10",
        "fluctuation": "step11",
        "fluctuation-baseline": "step12",
        "temporal-scaling": "step13",
        "interval-position": "step14",
        "detectability": "step15",
        "observation-boundary": "step15",
        "threshold": "step16",
        "observation-threshold": "step16",
        "operational-threshold": "step16",
    }
    if raw in aliases:
        return aliases[raw]
    allowed = ", ".join("step{}".format(item.step) for item in PLOTTER_LIST)
    raise ValueError("Unknown plot selection. Use one of: {}.".format(allowed))


# =============================================================================
# Plotter discovery and validation
# =============================================================================


def candidate_plotting_dirs(
    explicit: Optional[Path],
    code_dir: Path,
) -> List[Path]:
    if explicit is not None:
        candidates = [explicit.expanduser()]
    else:
        candidates = [
            code_dir.parent / "plotting",  # canonical code/dynamics -> code/plotting
            code_dir / "plotting",
            code_dir,
            code_dir.parent,
        ]
    return dedupe_paths(candidates)


def compatibility_candidates(directory: Path, canonical_script: str) -> List[Path]:
    stem = Path(canonical_script).stem
    matches: List[Path] = []
    for pattern in ("{}*.py".format(stem), "{}_v*.py".format(stem)):
        for candidate in sorted(directory.glob(pattern)):
            if candidate.name == canonical_script:
                continue
            if candidate.is_file() and candidate not in matches:
                matches.append(candidate)
    return matches


def validate_plotter_source(path: Path, plotter: Plotter) -> Tuple[bool, str]:
    try:
        source = path.read_text(encoding="utf-8", errors="replace").lower()
    except OSError as exc:
        return False, "cannot read source: {}".format(exc)

    has_backend = (
        "import matplotlib" in source
        or "from matplotlib" in source
        or "import plotly" in source
    )
    if not has_backend:
        return False, "no plotting-library import found"
    if plotter.input_option.lower() not in source:
        return False, "missing declared input option {}".format(plotter.input_option)
    if plotter.output_option is not None and plotter.output_option.lower() not in source:
        return False, "missing declared output option {}".format(plotter.output_option)
    return True, "ok"


def resolve_plotter(
    plotter: Plotter,
    plotting_dirs: Sequence[Path],
) -> Tuple[Path, bool]:
    checked: List[str] = []
    invalid: List[str] = []
    for directory in plotting_dirs:
        canonical = directory / plotter.script
        checked.append(str(canonical))
        if canonical.is_file():
            valid, reason = validate_plotter_source(canonical, plotter)
            if valid:
                return canonical.resolve(), False
            invalid.append("{}: {}".format(canonical, reason))

        for candidate in compatibility_candidates(directory, plotter.script):
            checked.append(str(candidate))
            valid, reason = validate_plotter_source(candidate, plotter)
            if valid:
                return candidate.resolve(), True
            invalid.append("{}: {}".format(candidate, reason))

    message = (
        "No valid plotter found: {}.\nChecked:\n{}".format(
            plotter.script,
            "\n".join("  - {}".format(item) for item in checked),
        )
    )
    if invalid:
        message += "\nRejected:\n" + "\n".join(
            "  - {}".format(item) for item in invalid
        )
    raise FileNotFoundError(message)


# =============================================================================
# Input-directory resolution
# =============================================================================


def discover_analysis_dir(root: Path, plotter: Plotter) -> List[Path]:
    name = plotter.analysis_dir_name
    direct = [
        root / name,
        root / "longitudinal" / name,
        root / "pseudo" / name,
        root / "shared" / name,
        root / "comparison" / name,
        root / "compare" / name,
    ]
    matches = [path for path in direct if path.is_dir()]
    if root.is_dir() and root.name == name:
        matches.insert(0, root)
    if root.is_dir():
        try:
            for candidate in root.rglob(name):
                if candidate.is_dir():
                    matches.append(candidate)
                    if len(matches) >= 20:
                        break
        except OSError:
            pass
    return dedupe_paths(matches)


def resolve_input_dir(raw: str, plotter: Plotter) -> Path:
    entered = Path(strip_matching_quotes(raw)).expanduser()
    if entered.exists() and entered.is_dir() and entered.name == plotter.analysis_dir_name:
        return entered.resolve()

    roots: List[Path] = []
    if entered.exists():
        if not entered.is_dir():
            raise NotADirectoryError(
                "Input path is not a directory: {}".format(entered.resolve())
            )
        roots.append(entered.resolve())
    else:
        parent = entered.parent
        if parent.exists() and parent.is_dir():
            roots.append(parent.resolve())
        else:
            nearest = nearest_existing_parent(entered)
            if nearest is not None:
                roots.append(nearest)

    if not entered.exists() and entered.parent.name == entered.name and entered.parent.is_dir():
        roots.insert(0, entered.parent.resolve())

    found: List[Path] = []
    for root in dedupe_paths(roots):
        found.extend(discover_analysis_dir(root, plotter))
    found = dedupe_paths(found)

    if len(found) == 1:
        print(
            "[AUTO] Resolved Step {} analysis directory:\n       {}".format(
                plotter.step, found[0]
            ),
            flush=True,
        )
        return found[0]

    if len(found) > 1:
        print(
            "[INFO] Multiple candidate directories found for Step {}:".format(
                plotter.step
            ),
            flush=True,
        )
        for index, candidate in enumerate(found, 1):
            print("  [{}] {}".format(index, candidate), flush=True)
        while True:
            raw_choice = input("Choose directory [1-{}]: ".format(len(found))).strip()
            try:
                choice = int(raw_choice)
            except ValueError:
                print("Please enter a valid number.", flush=True)
                continue
            if 1 <= choice <= len(found):
                return found[choice - 1]

    raise FileNotFoundError(
        "Could not locate the selected Step analysis directory.\n"
        "Entered: {}\nExpected directory name: {}\n"
        "Enter either the exact directory or the parent results root.".format(
            entered.absolute(), plotter.analysis_dir_name
        )
    )


def prompt_input_dir(plotter: Plotter) -> Path:
    while True:
        raw = input("Input analysis directory or results root: ").strip()
        if not raw:
            print("A directory is required.", flush=True)
            continue
        try:
            return resolve_input_dir(raw, plotter)
        except (FileNotFoundError, NotADirectoryError) as exc:
            print("[ERROR] {}".format(exc), flush=True)


# =============================================================================
# Forwarded arguments and dispatch
# =============================================================================


def forwarded_has_option(args: Sequence[str], option_names: Sequence[str]) -> bool:
    names = set(option_names)
    forwarded = list(args)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    return any(token.split("=", 1)[0] in names for token in forwarded)


def prompt_optional_repertoire_dir() -> Optional[Path]:
    print()
    print("Step 1 panel A can also use the original repertoire files.")
    print("Press Enter to skip panel A, or provide the raw repertoire directory.")
    while True:
        raw = input("Original repertoire directory [optional]: ").strip()
        if not raw:
            return None
        path = Path(strip_matching_quotes(raw)).expanduser()
        if path.exists() and path.is_dir():
            return path.resolve()
        print("[ERROR] Repertoire directory does not exist: {}".format(path), flush=True)


def check_forwarded_args(args: Sequence[str]) -> List[str]:
    forwarded = list(args)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    for token in forwarded:
        name = token.split("=", 1)[0]
        if name in PROTECTED_ARGS:
            raise ValueError(
                "{} is controlled by the plotting orchestrator. Choose the input "
                "directory with the orchestrator's --input-dir; output is always "
                "<step-output>/figures.".format(name)
            )
    return forwarded


def build_plotter_args(
    plotter: Plotter,
    input_dir: Path,
    extra_args: Sequence[str],
) -> List[str]:
    figures_dir = input_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    command_args = [plotter.input_option, str(input_dir)]
    if plotter.output_option is not None:
        command_args.extend([plotter.output_option, str(figures_dir)])
    command_args.extend(plotter.default_args)
    command_args.extend(check_forwarded_args(extra_args))
    return command_args


def dispatch(
    plotter_key: str,
    plotting_dirs: Sequence[Path],
    input_dir: Path,
    extra_args: Sequence[str],
    dry_dispatch: bool,
) -> int:
    plotter = PLOTTERS[plotter_key]
    plotter_path, used_compat = resolve_plotter(plotter, plotting_dirs)
    plotter_args = build_plotter_args(plotter, input_dir, extra_args)
    figures_dir = input_dir / "figures"
    command = [
        str(Path(sys.executable).resolve()),
        "-u",
        str(plotter_path),
        *plotter_args,
    ]

    print("\nReady to launch")
    print("---------------")
    print("Plotter : {}".format(plotter_path.name))
    print("Title   : {}".format(plotter.title))
    print("Script  : {}".format(plotter_path))
    if used_compat:
        print("Mode    : compatible filename")
    print("Input   : {}".format(input_dir))
    print("Output  : {}".format(figures_dir))
    print("\nCommand:")
    print("$ " + shell_join(command), flush=True)

    if dry_dispatch:
        print("\n[DRY DISPATCH] Plotter resolved; nothing executed.")
        return 0

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    try:
        completed = subprocess.run(
            command,
            cwd=str(plotter_path.parent),
            env=env,
            check=False,
        )
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    return int(completed.returncode)


# =============================================================================
# CLI and main
# =============================================================================


def parse_args(
    argv: Optional[Sequence[str]] = None,
) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser(
        description="Launch ClonoDynamics plotting modules.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--plot",
        default=None,
        help="Direct selector, e.g. step9, step15, step16.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Step analysis directory or parent results root.",
    )
    parser.add_argument("--plotting-dir", type=Path, default=None)
    parser.add_argument("--dry-dispatch", action="store_true")
    parser.add_argument("--list", action="store_true")
    return parser.parse_known_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    banner()
    args, extra_args = parse_args(argv)
    code_dir = Path(__file__).resolve().parent
    plotting_dirs = candidate_plotting_dirs(args.plotting_dir, code_dir)

    if args.list:
        print_menu()
        return 0

    if args.plot is None:
        if extra_args:
            raise ValueError(
                "Extra plotter arguments were supplied before selecting a plotter. "
                "Use --plot in direct mode or run interactively without extra arguments."
            )
        plotter_key = prompt_plot()
        if plotter_key is None:
            print("\nExit Figures.")
            return 0
    else:
        plotter_key = normalize_plot(args.plot)

    selected = PLOTTERS[plotter_key]
    resolve_plotter(selected, plotting_dirs)

    if args.input_dir is None:
        input_dir = prompt_input_dir(selected)
    else:
        input_dir = resolve_input_dir(str(args.input_dir), selected)

    effective_extra = list(extra_args)
    if (
        plotter_key == "step1"
        and args.plot is None
        and not forwarded_has_option(
            effective_extra, ("--repertoires-dir", "--repertoires_dir")
        )
    ):
        repertoire_dir = prompt_optional_repertoire_dir()
        if repertoire_dir is not None:
            effective_extra.extend(["--repertoires-dir", str(repertoire_dir)])

    return dispatch(
        plotter_key,
        plotting_dirs,
        input_dir,
        effective_extra,
        bool(args.dry_dispatch),
    )


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        print("\n[ERROR] {}".format(exc), file=sys.stderr)
        raise SystemExit(2)
