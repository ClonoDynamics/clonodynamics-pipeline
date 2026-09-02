#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clonodynamics.py
================

Top-level launcher for ClonoDynamics.

The interface separates:

    ANALYSIS
        scientific workflows and finalized analysis outputs;

    FIGURES
        plot-only scripts that visualize existing outputs.

Analysis architecture
---------------------
    Steps 1-6             clonodynamics_orchestrator_steps_1_6.py
    Steps 7-8             clonodynamics_orchestrator_pseudo_reference.py
    Steps 9,11,13-16      clonodynamics_orchestrator_longitudinal_dynamics.py
    Steps 10,12           clonodynamics_orchestrator_longitudinal_vs_pseudo.py

Step 16 performs operational observation-threshold robustness at
alpha = 0.10, 0.05, 0.025 and 0.01, reusing P09/P11/P13.

Figure architecture
-------------------
    clonodynamics_orchestrator_plotting.py

The plotting orchestrator includes the dedicated Step-16 plotter and applies
its production defaults: PDF, Arial 8, and reference-alpha confidence bands.

Interactive use
---------------
    python3 clonodynamics.py

Direct Step-16 analysis through the master launcher
---------------------------------------------------
    python3 clonodynamics.py --workflow longitudinal -- \
        --workflow 16 \
        --results-root ./dataset_longitudinal_clonodynamics_results/longitudinal \
        --dataset-label healthy \
        --yes \
        --non-interactive

Direct Step-16 plotting through the master launcher
---------------------------------------------------
    python3 clonodynamics.py \
        --mode figures \
        --plot step16 \
        --input-dir ./dataset_longitudinal_clonodynamics_results/longitudinal/16-observation_threshold_robustness

Arguments following `--` are forwarded to the selected orchestrator.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


VERSION = "v7-analysis-figures-master-step16-2026-09-01"
PLOTTING_ORCHESTRATOR = "clonodynamics_orchestrator_plotting.py"
PLOTTING_ORCHESTRATOR_ALIASES = (
    "clonodynamics_orchestrator_plotting_v6_1_step6.py",
    "clonodynamics_orchestrator_plotting_v5_fixed.py",
    "clonodynamics_plotting_orchestrator.py",
    "orchestrator_plotting.py",
)


@dataclass(frozen=True)
class Workflow:
    key: str
    menu_number: int
    title: str
    step_summary: str
    canonical_script: str
    aliases: Tuple[str, ...] = ()


WORKFLOWS: Dict[str, Workflow] = {
    "shared": Workflow(
        key="shared",
        menu_number=1,
        title="Repertoire to latent transitions",
        step_summary="Steps 1-6",
        canonical_script="clonodynamics_orchestrator_steps_1_6.py",
        aliases=("clonodynamics_orchestrator_steps_1_6_v2.py",),
    ),
    "pseudo": Workflow(
        key="pseudo",
        menu_number=2,
        title="Build pseudo-longitudinal reference",
        step_summary="Steps 7-8",
        canonical_script="clonodynamics_orchestrator_pseudo_reference.py",
    ),
    "longitudinal": Workflow(
        key="longitudinal",
        menu_number=3,
        title="Analyze genuine longitudinal dynamics",
        step_summary="Steps 9,11,13,14,15,16",
        canonical_script="clonodynamics_orchestrator_longitudinal_dynamics.py",
        aliases=("clonodynamics_orchestrator_longitudinal_dynamics_v2.py",),
    ),
    "compare": Workflow(
        key="compare",
        menu_number=4,
        title="Compare longitudinal versus pseudo",
        step_summary="Steps 10,12",
        canonical_script="clonodynamics_orchestrator_longitudinal_vs_pseudo.py",
        aliases=(
            "orchestrator_longitudinal_vs_pseudo.py",
            "orchestrator_longitudinal_vs_pseudo_v2.py",
        ),
    ),
}

MENU_TO_KEY = {
    workflow.menu_number: workflow.key for workflow in WORKFLOWS.values()
}


# =============================================================================
# Utilities and menus
# =============================================================================


def shell_join(command: Sequence[str]) -> str:
    try:
        return shlex.join([str(value) for value in command])
    except AttributeError:
        return " ".join(shlex.quote(str(value)) for value in command)


def sanitize_forwarded_args(args: Sequence[str]) -> List[str]:
    forwarded = list(args)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    return forwarded


def master_banner() -> None:
    print()
    print("=" * 78)
    print("ClonoDynamics")
    print("Noise-aware longitudinal clonotype dynamics")
    print("Master interface | {}".format(VERSION))
    print("=" * 78)
    print()
    print("  ANALYSIS  run scientific workflows and generate result tables.")
    print("  FIGURES   render figures from finalized analysis outputs.")


def print_mode_menu() -> None:
    print("\nWhat do you want to do?\n")
    print(" [1] Analysis")
    print(" [2] Figures")
    print(" [3] Exit")


def prompt_mode() -> Optional[str]:
    print_mode_menu()
    aliases = {
        "1": "analysis",
        "analysis": "analysis",
        "a": "analysis",
        "2": "figures",
        "figure": "figures",
        "figures": "figures",
        "plot": "figures",
        "plots": "figures",
        "f": "figures",
        "3": None,
        "exit": None,
        "quit": None,
        "q": None,
    }
    while True:
        raw = input("\nSelection [1-3]: ").strip().lower()
        if raw in aliases:
            return aliases[raw]
        print("Please enter 1, 2, or 3.")


def normalize_mode(value: str) -> str:
    raw = str(value).strip().lower().replace("_", "-")
    aliases = {
        "analysis": "analysis",
        "analyse": "analysis",
        "analyze": "analysis",
        "a": "analysis",
        "1": "analysis",
        "figures": "figures",
        "figure": "figures",
        "plots": "figures",
        "plot": "figures",
        "f": "figures",
        "2": "figures",
    }
    if raw not in aliases:
        raise ValueError("Unknown mode. Use 'analysis' or 'figures'.")
    return aliases[raw]


def analysis_banner() -> None:
    print()
    print("=" * 74)
    print("ClonoDynamics | Analysis workflows")
    print("=" * 74)


def print_analysis_menu() -> None:
    print("\nWhat do you want to run?\n")
    for workflow in sorted(WORKFLOWS.values(), key=lambda item: item.menu_number):
        print(" [{}] {}".format(workflow.menu_number, workflow.title))
        print("     {}".format(workflow.step_summary))
        print()
    print(" [5] Exit")


def prompt_workflow() -> Optional[str]:
    print_analysis_menu()
    while True:
        raw = input("\nSelection [1-5]: ").strip()
        if raw == "5":
            return None
        try:
            number = int(raw)
        except ValueError:
            print("Please enter 1, 2, 3, 4, or 5.")
            continue
        if number in MENU_TO_KEY:
            return MENU_TO_KEY[number]
        print("Please enter 1, 2, 3, 4, or 5.")


def normalize_workflow(value: str) -> str:
    raw = str(value).strip().lower().replace("_", "-")
    aliases = {
        "shared": "shared",
        "1": "shared",
        "1-6": "shared",
        "steps1-6": "shared",
        "steps-1-6": "shared",
        "pseudo": "pseudo",
        "2": "pseudo",
        "7-8": "pseudo",
        "steps7-8": "pseudo",
        "steps-7-8": "pseudo",
        "pseudo-reference": "pseudo",
        "longitudinal": "longitudinal",
        "3": "longitudinal",
        "long": "longitudinal",
        "genuine": "longitudinal",
        "genuine-longitudinal": "longitudinal",
        "compare": "compare",
        "comparison": "compare",
        "4": "compare",
        "10,12": "compare",
        "10-12": "compare",
        "longitudinal-vs-pseudo": "compare",
    }
    if raw not in aliases:
        raise ValueError(
            "Unknown workflow. Use shared, pseudo, longitudinal, or compare."
        )
    return aliases[raw]


# =============================================================================
# Script resolution and dispatch
# =============================================================================


def resolve_code_dir(raw: Optional[Path]) -> Path:
    path = (
        Path(__file__).resolve().parent
        if raw is None
        else raw.expanduser().resolve(strict=True)
    )
    if not path.is_dir():
        raise NotADirectoryError(path)
    return path


def resolve_script(
    code_dir: Path,
    canonical_script: str,
    aliases: Sequence[str],
    label: str,
) -> Tuple[Path, bool]:
    canonical = code_dir / canonical_script
    if canonical.is_file():
        return canonical.resolve(), False
    for alias in aliases:
        candidate = code_dir / alias
        if candidate.is_file():
            return candidate.resolve(), True
    raise FileNotFoundError(
        "Required {} was not found.\nCode directory: {}\nChecked:\n{}".format(
            label,
            code_dir,
            "\n".join(
                "  - {}".format(name) for name in (canonical_script, *aliases)
            ),
        )
    )


def resolve_analysis_orchestrator(
    workflow: Workflow,
    code_dir: Path,
) -> Tuple[Path, bool]:
    return resolve_script(
        code_dir,
        workflow.canonical_script,
        workflow.aliases,
        "orchestrator for '{}'".format(workflow.title),
    )


def resolve_plotting_orchestrator(code_dir: Path) -> Tuple[Path, bool]:
    return resolve_script(
        code_dir,
        PLOTTING_ORCHESTRATOR,
        PLOTTING_ORCHESTRATOR_ALIASES,
        "plotting orchestrator",
    )


def build_command(script: Path, forwarded_args: Sequence[str]) -> List[str]:
    return [
        str(Path(sys.executable).resolve()),
        "-u",
        str(script),
        *sanitize_forwarded_args(forwarded_args),
    ]


def run_child(command: Sequence[str], cwd: Path, dry_dispatch: bool) -> int:
    print("\nDelegating:")
    print("$ " + shell_join(command))
    if dry_dispatch:
        print("\n[DRY DISPATCH] Child command resolved but not launched.")
        return 0
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    try:
        completed = subprocess.run(
            list(command), cwd=str(cwd), env=env, check=False
        )
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    return int(completed.returncode)


def dispatch_analysis(
    workflow_key: str,
    code_dir: Path,
    forwarded_args: Sequence[str],
    dry_dispatch: bool,
) -> int:
    workflow = WORKFLOWS[workflow_key]
    orchestrator, used_alias = resolve_analysis_orchestrator(workflow, code_dir)
    if used_alias:
        print("\n[compat] Using orchestrator alias: {}".format(orchestrator.name))
    command = build_command(orchestrator, forwarded_args)
    print("\nSelected Analysis workflow")
    print("--------------------------")
    print("Analysis : {}".format(workflow.title))
    print("Steps    : {}".format(workflow.step_summary))
    print("Launcher : {}".format(orchestrator))
    print("Python   : {}".format(Path(sys.executable).resolve()))
    if forwarded_args:
        print("Forwarded: " + shell_join(sanitize_forwarded_args(forwarded_args)))
    return run_child(command, code_dir, dry_dispatch)


def dispatch_figures(
    code_dir: Path,
    plot: Optional[str],
    input_dir: Optional[Path],
    plotting_dir: Optional[Path],
    forwarded_args: Sequence[str],
    dry_dispatch: bool,
) -> int:
    """Run the plotting orchestrator in-process for reliable interactive input."""
    orchestrator, used_alias = resolve_plotting_orchestrator(code_dir)
    if used_alias:
        print("\n[compat] Using plotting orchestrator alias: {}".format(orchestrator.name))

    child_args: List[str] = []
    if plot is not None:
        child_args.extend(["--plot", str(plot)])
    if input_dir is not None:
        child_args.extend(["--input-dir", str(input_dir)])
    if plotting_dir is not None:
        child_args.extend(["--plotting-dir", str(plotting_dir)])
    if dry_dispatch:
        child_args.append("--dry-dispatch")
    forwarded = sanitize_forwarded_args(forwarded_args)
    if forwarded:
        child_args.append("--")
        child_args.extend(forwarded)

    print("\nSelected layer")
    print("--------------")
    print("Mode     : Figures")
    print("Launcher : {}".format(orchestrator))
    if plot is not None:
        print("Plot     : {}".format(plot))
    if input_dir is not None:
        print("Input    : {}".format(input_dir))
    if plotting_dir is not None:
        print("Plot dir : {}".format(plotting_dir))
    print("[INFO] Figures interface running in the current process.", flush=True)

    module_name = "_clonodynamics_plotting_orchestrator_runtime"
    spec = importlib.util.spec_from_file_location(module_name, orchestrator)
    if spec is None or spec.loader is None:
        raise ImportError("Could not load plotting orchestrator: {}".format(orchestrator))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    main_fn = getattr(module, "main", None)
    if not callable(main_fn):
        raise AttributeError(
            "Plotting orchestrator does not expose main(): {}".format(orchestrator)
        )
    result = main_fn(child_args)
    return int(result or 0)


# =============================================================================
# CLI and main
# =============================================================================


def parse_args(
    argv: Optional[Sequence[str]] = None,
) -> Tuple[argparse.Namespace, List[str]]:
    parser = argparse.ArgumentParser(
        description="Top-level Analysis/Figures launcher for ClonoDynamics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mode", default=None)
    parser.add_argument("--workflow", default=None)
    parser.add_argument(
        "--plot",
        default=None,
        help="Direct Figure selection, e.g. step9, step15, step16.",
    )
    parser.add_argument("--code-dir", type=Path, default=None)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Analysis directory or results root for direct Figure dispatch.",
    )
    parser.add_argument("--plotting-dir", type=Path, default=None)
    parser.add_argument("--dry-dispatch", action="store_true")
    return parser.parse_known_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    master_banner()
    args, forwarded = parse_args(argv)
    code_dir = resolve_code_dir(args.code_dir)

    if args.mode is None:
        if args.workflow is not None and args.plot is not None:
            raise ValueError("Use either --workflow or --plot, not both.")
        if args.workflow is not None:
            mode = "analysis"
        elif args.plot is not None:
            mode = "figures"
        else:
            if forwarded:
                raise ValueError(
                    "Forwarded arguments were supplied before choosing Analysis "
                    "or Figures. Use --mode/--workflow/--plot or run interactively."
                )
            mode = prompt_mode()
            if mode is None:
                print("\nExit.")
                return 0
    else:
        mode = normalize_mode(args.mode)

    if mode == "analysis":
        if args.plot is not None:
            raise ValueError("--plot can only be used with --mode figures.")
        if args.workflow is None:
            if forwarded:
                raise ValueError(
                    "Analysis-child arguments were provided without --workflow."
                )
            analysis_banner()
            workflow_key = prompt_workflow()
            if workflow_key is None:
                print("\nExit Analysis.")
                return 0
        else:
            workflow_key = normalize_workflow(args.workflow)
        return dispatch_analysis(
            workflow_key,
            code_dir,
            forwarded,
            bool(args.dry_dispatch),
        )

    if args.workflow is not None:
        raise ValueError("--workflow can only be used with --mode analysis.")
    return dispatch_figures(
        code_dir,
        args.plot,
        args.input_dir,
        args.plotting_dir,
        forwarded,
        bool(args.dry_dispatch),
    )


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        print("\n[ERROR] {}".format(exc), file=sys.stderr)
        raise SystemExit(2)
