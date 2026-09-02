#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clonodynamics_orchestrator_longitudinal_dynamics.py
====================================================

Dependency-aware orchestrator for the genuine-longitudinal ClonoDynamics
analysis branch.

Managed public Steps
--------------------
    9   replicate-decoupled forward drift
    11  longitudinal fluctuation dynamics
    13  temporal fluctuation scaling
    14  interval-position and anchoring sensitivity
    15  observation-boundary transition-domain sensitivity
    16  operational observation-threshold robustness

Step 16 reclassifies endpoint states at alpha = 0.10, 0.05, 0.025 and 0.01
without refitting latent abundance, and reuses the current P09/P11/P13
implementations. Its upstream inputs are Step-2 outputs, Step-4 trajectories
and Step-5 transitions.

Steps 10 and 12 remain in the longitudinal-versus-pseudo comparison
orchestrator.

Typical interactive use
-----------------------
    python3 clonodynamics_orchestrator_longitudinal_dynamics.py

Run Step 16 only
----------------
    python3 clonodynamics_orchestrator_longitudinal_dynamics.py \
        --workflow 16 \
        --results-root ./dataset_longitudinal_clonodynamics_results/longitudinal \
        --dataset-label healthy \
        --yes \
        --non-interactive

Complete branch
---------------
    python3 clonodynamics_orchestrator_longitudinal_dynamics.py \
        --workflow complete \
        --results-root ./dataset_longitudinal_clonodynamics_results/longitudinal \
        --dataset-label healthy \
        --yes \
        --non-interactive

Fresh/resume policy for Step 16
-------------------------------
A complete Step-16 result is skipped in resume mode. An incomplete Step-16
result is restarted at public-Step granularity, because P16 coordinates four
alpha-specific P09/P11/P13 runs and their common support/core must remain
internally coherent.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


VERSION = "v2-longitudinal-dynamics-step16-2026-09-01"

MANAGED_STEPS = (9, 11, 13, 14, 15, 16)

SCRIPT_NAMES = {
    9: "9-replicate_decoupled_forward_drift.py",
    11: "11-longitudinal_fluctuation_dynamics.py",
    13: "13-temporal_fluctuation_scaling.py",
    14: "14-interval_position_structure.py",
    15: "15-detectability_boundary_sensitivity.py",
    16: "16-analyze_observation_threshold_robustness.py",
}

OUTDIR_NAMES = {
    9: "9-replicate_decoupled_forward_drift",
    11: "11-longitudinal_fluctuation_dynamics",
    13: "13-temporal_fluctuation_scaling",
    14: "14-interval_position_structure",
    15: "15-detectability_boundary_sensitivity",
    16: "16-observation_threshold_robustness",
}

DIR2 = "2-noise_aware_latent_inference"
DIR4 = "4-latent_trajectory_construction"
DIR5 = "5-latent_transition_construction"
TRAJECTORY_BASENAME = "latent_trajectories_long.parquet"
TRANSITION_BASENAME = "latent_transitions.parquet"
STEP2_PRIMARY_MARKER = "per_clone_latent_subject.parquet"

DEFAULT_DT_VALUES = (1, 2, 3, 4, 5)
P16_ALPHAS = (0.10, 0.05, 0.025, 0.01)
P16_REFERENCE_ALPHA = 0.05

COMPLETION_MARKERS = {
    9: (
        "00_run_config.json",
        "forward_drift_by_bin.csv",
        "forward_longitudinal_report.json",
    ),
    11: (
        "00_analysis_metadata.json",
        "09_subject_bin_dynamics.parquet",
        "10_binned_dynamics_long.parquet",
    ),
    13: (
        "00_run_config.json",
        "05_model_comparison_by_bin.csv",
        "07_interpretation_aid.csv",
    ),
    14: (
        "00_run_config.json",
        "06_interval_position_permutation_tests.csv",
        "10_analysis_summary.csv",
        "DEBUG_REPORT.txt",
    ),
    15: (
        "00_run_config.json",
        "15_fluctuation_mode_robustness_summary.csv",
        "19_temporal_scaling_sensitivity_summary.csv",
        "20_step15_summary.csv",
        "DEBUG_REPORT.txt",
    ),
    16: (
        "00_run_config.json",
        "03_class_composition_by_alpha_dt.csv",
        "05_common_xstar_support.csv",
        "06_tt_retention_by_alpha.csv",
        "09_combined_p11_binned_dynamics.csv",
        "12_combined_p13_model_comparison_by_bin.csv",
        "14_cross_alpha_common_core.csv",
        "15_combined_p09_forward_drift_by_bin.csv",
        "DEBUG_REPORT.txt",
    ),
}


# =============================================================================
# Generic utilities
# =============================================================================


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def shell_join(command: Sequence[str]) -> str:
    try:
        return shlex.join([str(x) for x in command])
    except AttributeError:
        return " ".join(shlex.quote(str(x)) for x in command)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def outdir(root: Path, step: int) -> Path:
    return root / OUTDIR_NAMES[int(step)]


def canonical_step2_dir(root: Path) -> Path:
    return root / DIR2


def canonical_trajectories(root: Path) -> Path:
    return root / DIR4 / TRAJECTORY_BASENAME


def canonical_transitions(root: Path) -> Path:
    return root / DIR5 / TRANSITION_BASENAME


def dt_csv(values: Sequence[int]) -> str:
    return ",".join(str(int(v)) for v in values)


def alpha_text(values: Sequence[float]) -> str:
    return ",".join("{:g}".format(float(v)) for v in values)


def normalize_steps(values: Sequence[int]) -> Tuple[int, ...]:
    steps = tuple(sorted(set(int(v) for v in values)))
    if not steps:
        raise ValueError("At least one longitudinal-dynamics Step is required.")
    invalid = [v for v in steps if v not in MANAGED_STEPS]
    if invalid:
        raise ValueError(
            "Unsupported Steps: {}. Allowed: {}".format(invalid, MANAGED_STEPS)
        )
    return steps


def parse_step_csv(text: str) -> Tuple[int, ...]:
    tokens = [
        token.strip()
        for token in str(text).replace(";", ",").split(",")
        if token.strip()
    ]
    try:
        return normalize_steps([int(token) for token in tokens])
    except Exception as exc:
        raise ValueError(
            "Custom Steps must be comma-separated integers from "
            "9,11,13,14,15,16."
        ) from exc


def parse_dt_values(text: str) -> Tuple[int, ...]:
    values = tuple(
        sorted(
            set(
                int(token.strip())
                for token in str(text).replace(";", ",").split(",")
                if token.strip()
            )
        )
    )
    if not values or any(v <= 0 for v in values):
        raise ValueError("dt values must be positive integers.")
    return values


# =============================================================================
# Workflow selection
# =============================================================================


def workflow_from_name(value: str) -> Tuple[int, ...]:
    key = str(value).strip().lower().replace(" ", "")
    mapping = {
        "complete": (9, 11, 13, 14, 15, 16),
        "all": (9, 11, 13, 14, 15, 16),
        "1": (9, 11, 13, 14, 15, 16),
        "core": (9, 11, 13),
        "2": (9, 11, 13),
        "sensitivity": (15, 16),
        "sensitivities": (15, 16),
        "15-16": (15, 16),
        "15,16": (15, 16),
        "fluctuation": (11, 13, 14),
        "11-14": (11, 13, 14),
        "11,13,14": (11, 13, 14),
        "9": (9,),
        "11": (11,),
        "13": (13,),
        "14": (14,),
        "15": (15,),
        "16": (16,),
        "observation-threshold": (16,),
        "threshold": (16,),
    }
    if key in mapping:
        return mapping[key]
    return parse_step_csv(key)


def steps_label(steps: Sequence[int]) -> str:
    return "Steps " + ",".join(str(int(v)) for v in steps)


def prompt_one_step() -> Tuple[int, ...]:
    options = (
        (1, 9, "replicate-decoupled forward drift"),
        (2, 11, "fluctuation dynamics"),
        (3, 13, "temporal scaling"),
        (4, 14, "interval-position structure"),
        (5, 15, "observation-boundary transition-domain sensitivity"),
        (6, 16, "operational observation-threshold robustness"),
    )
    print("\nOne longitudinal-dynamics Step:")
    for menu, step, title in options:
        print("  [{}] Step {:>2} — {}".format(menu, step, title))
    mapping = {str(menu): step for menu, step, _title in options}
    while True:
        raw = input("Selection: ").strip()
        if raw in mapping:
            return (mapping[raw],)
        if raw in {str(step) for _menu, step, _title in options}:
            return (int(raw),)
        print("Please select 1-6 or enter the Step number.")


def prompt_workflow() -> Tuple[int, ...]:
    print("\nLongitudinal-dynamics workflow:")
    print("  [1] Complete       — Steps 9,11,13,14,15,16")
    print("  [2] Core dynamics  — Steps 9,11,13")
    print("  [3] One Step")
    print("  [4] Fluctuation     — Steps 11,13,14")
    print("  [5] Sensitivities   — Steps 15,16")
    print("  [6] Custom Step selection")
    while True:
        raw = input("Selection [1]: ").strip()
        if raw in {"", "1"}:
            return (9, 11, 13, 14, 15, 16)
        if raw == "2":
            return (9, 11, 13)
        if raw == "3":
            return prompt_one_step()
        if raw == "4":
            return (11, 13, 14)
        if raw == "5":
            return (15, 16)
        if raw == "6":
            while True:
                custom = input(
                    "Steps (comma-separated subset of 9,11,13,14,15,16): "
                ).strip()
                try:
                    return parse_step_csv(custom)
                except Exception as exc:
                    print("Invalid selection: {}".format(exc))
        print("Please enter 1-6.")


# =============================================================================
# Results-root and provenance
# =============================================================================


def looks_like_longitudinal_root(path: Path) -> bool:
    indicators = (
        "00_orchestrator_config.json",
        DIR2,
        DIR4,
        DIR5,
        OUTDIR_NAMES[9],
        OUTDIR_NAMES[11],
        OUTDIR_NAMES[13],
        OUTDIR_NAMES[14],
        OUTDIR_NAMES[15],
        OUTDIR_NAMES[16],
    )
    return any((path / item).exists() for item in indicators)


def normalize_results_root(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)
    if looks_like_longitudinal_root(path):
        return path
    child = path / "longitudinal"
    if child.is_dir() and looks_like_longitudinal_root(child):
        print("[root] Using longitudinal child:\n  {}".format(child))
        return child.resolve()
    return path


def prompt_results_root() -> Path:
    while True:
        raw = input("\nExisting longitudinal ClonoDynamics results root: ").strip()
        if not raw:
            print("A results-root path is required.")
            continue
        try:
            return normalize_results_root(Path(raw))
        except Exception as exc:
            print("Invalid results root: {}".format(exc))


def upstream_config(root: Path) -> Optional[Dict[str, object]]:
    path = root / "00_orchestrator_config.json"
    if not path.exists():
        return None
    try:
        return load_json(path)
    except Exception as exc:
        print("[warning] Could not read {}: {}".format(path, exc))
        return None


def validate_longitudinal_provenance(
    root: Path,
    config: Optional[Dict[str, object]],
) -> Dict[str, object]:
    config_path = root / "00_orchestrator_config.json"
    if config is None:
        print(
            "[warning] No readable Steps-1-6 orchestrator config was found; "
            "longitudinal provenance cannot be verified automatically."
        )
        return {
            "config_found": False,
            "config_path": str(config_path),
            "dataset_type": None,
        }
    dataset_type = config.get("dataset_type")
    if dataset_type is not None:
        normalized = str(dataset_type).strip().lower().replace("_", "-")
        if normalized != "longitudinal":
            raise ValueError(
                "This orchestrator accepts genuine longitudinal result trees only. "
                "Existing config reports dataset_type={!r}:\n  {}".format(
                    dataset_type, config_path
                )
            )
    return {
        "config_found": True,
        "config_path": str(config_path),
        "dataset_type": dataset_type,
    }


def inherited_alpha(config: Optional[Dict[str, object]]) -> Optional[float]:
    if not config:
        return None
    try:
        value = float(config.get("observability_alpha"))
    except Exception:
        return None
    return value if 0 < value < 1 else None


def previous_longitudinal_config(root: Path) -> Optional[Dict[str, object]]:
    path = root / "00_orchestrator_longitudinal_dynamics_config.json"
    if not path.exists():
        return None
    try:
        return load_json(path)
    except Exception:
        return None


def resolve_dataset_label(
    explicit: Optional[str],
    previous: Optional[Dict[str, object]],
    interactive: bool,
) -> str:
    if explicit:
        return str(explicit).strip()
    if previous and previous.get("dataset_label"):
        return str(previous["dataset_label"])
    if interactive:
        return input("\nDataset label [healthy]: ").strip() or "healthy"
    return "healthy"


def resolve_step15_alpha(
    selected_steps: Sequence[int],
    explicit: Optional[float],
    upstream: Optional[float],
    interactive: bool,
) -> Optional[float]:
    if 15 not in selected_steps:
        return upstream
    if explicit is not None:
        value = float(explicit)
        if not 0 < value < 1:
            raise ValueError("--detectability-alpha must satisfy 0 < alpha < 1")
        if upstream is not None and abs(value - upstream) > 1e-12:
            raise ValueError(
                "--detectability-alpha conflicts with upstream observability_alpha "
                "({} vs {}). Step 15 records the existing threshold; it does not "
                "redefine it.".format(value, upstream)
            )
        return value
    if upstream is not None:
        print("[config] Step-15 alpha inherited: {:g}".format(upstream))
        return upstream
    if interactive:
        while True:
            raw = input("\nUpstream operational alpha [0.05]: ").strip()
            if not raw:
                return 0.05
            try:
                value = float(raw)
            except ValueError:
                print("Alpha must be numeric.")
                continue
            if 0 < value < 1:
                return value
            print("Alpha must satisfy 0 < alpha < 1.")
    return 0.05


# =============================================================================
# Script and upstream artifact resolution
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


def required_script_steps(selected_steps: Sequence[int]) -> Tuple[int, ...]:
    required = set(int(v) for v in selected_steps)
    if 15 in required or 16 in required:
        required.update({9, 11, 13})
    return tuple(sorted(required))


def validate_scripts(
    code_dir: Path,
    selected_steps: Sequence[int],
) -> Dict[int, Path]:
    scripts: Dict[int, Path] = {}
    missing: List[Path] = []
    for step in required_script_steps(selected_steps):
        path = code_dir / SCRIPT_NAMES[step]
        scripts[step] = path
        if not path.is_file():
            missing.append(path)
    if missing:
        raise FileNotFoundError(
            "Required longitudinal-dynamics scripts are missing:\n"
            + "\n".join("  - {}".format(path) for path in missing)
        )
    return scripts


def validate_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError("{} missing/empty:\n  {}".format(label, path))
    return path


def validate_step_dir(
    path: Path,
    label: str,
    required_markers: Sequence[str],
) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError("{}: {}".format(label, path))
    missing = [marker for marker in required_markers if not (path / marker).exists()]
    if missing:
        raise FileNotFoundError(
            "{} is incomplete. Missing markers: {}\n  {}".format(
                label, missing, path
            )
        )
    return path


def file_candidates(root: Path, basename: str) -> List[Path]:
    return sorted(
        {
            path.resolve()
            for path in root.rglob(basename)
            if path.is_file() and path.stat().st_size > 0
        },
        key=lambda path: (len(path.parts), str(path)),
    )


def dir_candidates_from_markers(
    root: Path,
    markers: Sequence[str],
) -> List[Path]:
    parents = None
    for marker in markers:
        current = {
            path.parent.resolve()
            for path in root.rglob(marker)
            if path.is_file() and path.stat().st_size > 0
        }
        parents = current if parents is None else parents & current
    return sorted(parents or set(), key=lambda path: (len(path.parts), str(path)))


def choose_file_interactively(label: str, candidates: Sequence[Path]) -> Path:
    print("\nRequired upstream input: {}".format(label))
    if candidates:
        for index, path in enumerate(candidates, 1):
            print("  [{}] {}".format(index, path))
        print("  [0] enter another path")
        while True:
            raw = input("Selection: ").strip()
            try:
                choice = int(raw)
            except ValueError:
                print("Enter a candidate number.")
                continue
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return candidates[choice - 1]
            print("Invalid candidate number.")
    while True:
        raw = input("Path for {}: ".format(label)).strip()
        try:
            return validate_file(Path(raw), label)
        except Exception as exc:
            print(exc)


def choose_dir_interactively(
    label: str,
    candidates: Sequence[Path],
    markers: Sequence[str],
) -> Path:
    print("\nRequired upstream directory: {}".format(label))
    if candidates:
        for index, path in enumerate(candidates, 1):
            print("  [{}] {}".format(index, path))
        print("  [0] enter another path")
        while True:
            raw = input("Selection: ").strip()
            try:
                choice = int(raw)
            except ValueError:
                print("Enter a candidate number.")
                continue
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return candidates[choice - 1]
            print("Invalid candidate number.")
    while True:
        raw = input("Path for {}: ".format(label)).strip()
        try:
            return validate_step_dir(Path(raw), label, markers)
        except Exception as exc:
            print(exc)


def resolve_file(
    root: Path,
    canonical: Path,
    basename: str,
    label: str,
    override: Optional[Path],
    interactive: bool,
) -> Path:
    if override is not None:
        return validate_file(override, label)
    if canonical.exists():
        return validate_file(canonical, label)
    candidates = file_candidates(root, basename)
    if len(candidates) == 1:
        print("[input] {}:\n  {}".format(label, candidates[0]))
        return candidates[0]
    if interactive:
        return choose_file_interactively(label, candidates)
    if not candidates:
        raise FileNotFoundError(
            "Could not find {}. Canonical path:\n  {}".format(label, canonical)
        )
    raise RuntimeError(
        "Multiple candidates found for {}; use an explicit override:\n{}".format(
            label, "\n".join("  - {}".format(path) for path in candidates)
        )
    )


def resolve_directory(
    root: Path,
    canonical: Path,
    label: str,
    markers: Sequence[str],
    override: Optional[Path],
    interactive: bool,
) -> Path:
    if override is not None:
        return validate_step_dir(override, label, markers)
    if canonical.is_dir():
        try:
            return validate_step_dir(canonical, label, markers)
        except Exception:
            pass
    candidates = dir_candidates_from_markers(root, markers)
    if len(candidates) == 1:
        print("[input] {}:\n  {}".format(label, candidates[0]))
        return candidates[0]
    if interactive:
        return choose_dir_interactively(label, candidates, markers)
    if not candidates:
        raise FileNotFoundError(
            "Could not find a complete {}. Canonical path:\n  {}".format(
                label, canonical
            )
        )
    raise RuntimeError(
        "Multiple candidates found for {}; use an explicit override:\n{}".format(
            label, "\n".join("  - {}".format(path) for path in candidates)
        )
    )


def resolve_preexisting_inputs(
    root: Path,
    selected_steps: Sequence[int],
    trajectories_override: Optional[Path],
    transitions_override: Optional[Path],
    step2_override: Optional[Path],
    step11_override: Optional[Path],
    step13_override: Optional[Path],
    interactive: bool,
) -> Dict[str, Path]:
    selected = set(int(v) for v in selected_steps)
    inputs: Dict[str, Path] = {}

    if selected & {9, 15, 16}:
        inputs["trajectories"] = resolve_file(
            root,
            canonical_trajectories(root),
            TRAJECTORY_BASENAME,
            "Step-4 latent_trajectories_long.parquet",
            trajectories_override,
            interactive,
        )

    if selected & {11, 14, 15, 16}:
        inputs["transitions"] = resolve_file(
            root,
            canonical_transitions(root),
            TRANSITION_BASENAME,
            "Step-5 latent_transitions.parquet",
            transitions_override,
            interactive,
        )

    if 16 in selected:
        inputs["step2_dir"] = resolve_directory(
            root,
            canonical_step2_dir(root),
            "Step-2 noise-aware latent-inference directory",
            (STEP2_PRIMARY_MARKER,),
            step2_override,
            interactive,
        )

    if 13 in selected and 11 not in selected:
        inputs["step11_dir"] = resolve_directory(
            root,
            outdir(root, 11),
            "Step-11 longitudinal fluctuation-dynamics directory",
            ("00_analysis_metadata.json", "10_binned_dynamics_long.parquet"),
            step11_override,
            interactive,
        )

    if 14 in selected and 11 not in selected and "step11_dir" not in inputs:
        inputs["step11_dir"] = resolve_directory(
            root,
            outdir(root, 11),
            "Step-11 longitudinal fluctuation-dynamics directory",
            ("00_analysis_metadata.json", "10_binned_dynamics_long.parquet"),
            step11_override,
            interactive,
        )

    if 14 in selected and 13 not in selected:
        inputs["step13_dir"] = resolve_directory(
            root,
            outdir(root, 13),
            "Step-13 temporal-scaling directory",
            ("00_run_config.json", "07_interpretation_aid.csv"),
            step13_override,
            interactive,
        )

    return inputs


def effective_inputs(
    root: Path,
    selected_steps: Sequence[int],
    resolved: Dict[str, Path],
) -> Dict[str, Path]:
    selected = set(int(v) for v in selected_steps)
    output = dict(resolved)
    if 11 in selected:
        output["step11_dir"] = outdir(root, 11)
    if 13 in selected:
        output["step13_dir"] = outdir(root, 13)
    return output


# =============================================================================
# Completion, cleanup, and run mode
# =============================================================================


def required_outputs(root: Path, step: int) -> List[Path]:
    return [outdir(root, step) / marker for marker in COMPLETION_MARKERS[step]]


def step_complete(root: Path, step: int) -> bool:
    for path in required_outputs(root, step):
        if not path.exists():
            return False
        if path.is_file() and path.stat().st_size <= 0:
            return False
    return True


def selected_outputs_exist(root: Path, selected_steps: Sequence[int]) -> bool:
    return any(
        outdir(root, step).is_dir() and any(outdir(root, step).iterdir())
        for step in selected_steps
    )


def clear_step_output(root: Path, step: int) -> None:
    path = outdir(root, step)
    if path.exists():
        shutil.rmtree(path)
    log = root / "logs" / "step{}.log".format(step)
    if log.exists():
        log.unlink()


def clear_selected_outputs(root: Path, selected_steps: Sequence[int]) -> None:
    for step in selected_steps:
        clear_step_output(root, step)


def choose_run_mode(
    root: Path,
    selected_steps: Sequence[int],
    restart: bool,
    resume: bool,
    interactive: bool,
) -> str:
    if restart and resume:
        raise ValueError("Use only one of --restart or --resume.")
    if not selected_outputs_exist(root, selected_steps):
        return "fresh"
    if restart:
        return "fresh"
    if resume:
        return "resume"
    if not interactive:
        raise RuntimeError(
            "Selected Step outputs already exist. Use --restart or --resume."
        )
    print("\nExisting outputs detected for {}:".format(steps_label(selected_steps)))
    print("  [1] fresh   — delete only selected Step outputs and rerun")
    print("  [2] resume  — skip complete Steps; restart incomplete Step 16 if needed")
    print("  [3] abort")
    while True:
        raw = input("Selection [2]: ").strip().lower()
        if raw in {"", "2", "resume", "r"}:
            return "resume"
        if raw in {"1", "fresh", "f"}:
            return "fresh"
        if raw in {"3", "abort", "a", "q", "quit"}:
            return "abort"
        print("Please enter 1, 2, or 3.")


# =============================================================================
# Command construction
# =============================================================================


def build_step9_command(
    script: Path,
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    n_bootstrap: int,
    seed: int,
    fresh: bool,
) -> List[str]:
    command = [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--trajectories", str(inputs["trajectories"]),
        "--dataset-label", dataset_label,
        "--outdir", str(outdir(root, 9)),
        "--dt", "1",
        "--n-bins", "20",
        "--binning", "quantile",
        "--q-low", "0.005",
        "--q-high", "0.995",
        "--n-bootstrap", str(int(n_bootstrap)),
        "--bootstrap-seed", str(int(seed)),
        "--parquet-compression", "zstd",
        "--parquet-row-group-size", "100000",
    ]
    if fresh:
        command.append("--restart")
    return command


def build_step11_command(
    script: Path,
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
) -> List[str]:
    return [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--transitions", str(inputs["transitions"]),
        "--dataset-label", dataset_label,
        "--outdir", str(outdir(root, 11)),
        "--representation", "latent",
        "--conditioning", "xstar",
        "--dt", *[str(int(v)) for v in dt_values],
        "--mode", "TT",
        "--shared-range", "intersection",
        "--n-bins", "30",
        "--min-n", "30",
        "--min-subjects", "2",
        "--n-bootstrap", str(int(n_bootstrap)),
        "--ci", "95",
        "--seed", str(int(seed)),
        "--streaming",
    ]


def build_step13_command(
    script: Path,
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    dt_values: Sequence[int],
    min_dt_points: int,
    n_bootstrap: int,
    seed: int,
) -> List[str]:
    return [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--input-dir", str(inputs["step11_dir"]),
        "--outdir", str(outdir(root, 13)),
        "--dataset-label", dataset_label,
        "--dt-values", dt_csv(dt_values),
        "--metrics", "var_dx,msd",
        "--estimators", "transition_weighted,equal_subject_weighted",
        "--bin-scope", "core",
        "--core-min-subject-fraction", "1.0",
        "--min-dt-points", str(int(min_dt_points)),
        "--n-bootstrap", str(int(n_bootstrap)),
        "--seed", str(int(seed)),
        "--ci", "95",
    ]


def build_step14_command(
    script: Path,
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
) -> List[str]:
    return [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--transitions", str(inputs["transitions"]),
        "--step11-dir", str(inputs["step11_dir"]),
        "--step13-dir", str(inputs["step13_dir"]),
        "--outdir", str(outdir(root, 14)),
        "--dataset-label", dataset_label,
        "--dt-values", dt_csv(dt_values),
        "--metrics", "var_dx,msd",
        "--n-perm", str(int(n_bootstrap)),
        "--n-boot-anchor", str(int(n_bootstrap)),
        "--min-interval-positions", "3",
        "--min-subjects-per-interval", "3",
        "--min-anchor-dt-points", "4",
        "--min-matched-anchor-subjects", "3",
        "--seed", str(int(seed)),
        "--streaming",
    ]


def build_step15_command(
    script: Path,
    scripts: Dict[int, Path],
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    detectability_alpha: float,
    dt_values: Sequence[int],
    min_dt_points: int,
    n_bootstrap: int,
    seed: int,
    fresh: bool,
    resume: bool,
) -> List[str]:
    command = [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--trajectories", str(inputs["trajectories"]),
        "--transitions", str(inputs["transitions"]),
        "--outdir", str(outdir(root, 15)),
        "--dataset-label", dataset_label,
        "--detectability-alpha", repr(float(detectability_alpha)),
        "--step9-script", str(scripts[9]),
        "--step11-script", str(scripts[11]),
        "--step13-script", str(scripts[13]),
        "--python-executable", str(Path(sys.executable).resolve()),
        "--dt-values", dt_csv(dt_values),
        "--forward-dt", "1",
        "--forward-n-bins", "20",
        "--fluctuation-n-bins", "30",
        "--min-n", "30",
        "--min-subjects", "2",
        "--n-bootstrap", str(int(n_bootstrap)),
        "--ci", "95",
        "--core-min-subject-fraction", "1.0",
        "--min-dt-points", str(int(min_dt_points)),
        "--seed", str(int(seed)),
        "--run-step11",
        "--run-step13",
        "--streaming",
    ]
    if fresh:
        command.append("--restart-forward")
    elif resume:
        command.append("--reuse-existing")
    return command


def build_step16_command(
    script: Path,
    scripts: Dict[int, Path],
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    dt_values: Sequence[int],
    n_bootstrap: int,
    seed: int,
) -> List[str]:
    """Build the validated P16 production command.

    The core arguments intentionally match the manually validated command.
    Internal P09/P11/P13 paths and global bootstrap/seed values are explicit to
    remove dependence on the current working directory and keep orchestrator
    provenance complete.
    """
    return [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--trajectories", str(inputs["trajectories"]),
        "--transitions", str(inputs["transitions"]),
        "--step2-dir", str(inputs["step2_dir"]),
        "--out-dir", str(outdir(root, 16)),
        "--dataset-label", dataset_label,
        "--alphas", "0.10", "0.05", "0.025", "0.01",
        "--reference-alpha", "0.05",
        "--dt", *[str(int(v)) for v in dt_values],
        "--seed", str(int(seed)),
        "--p09-script", str(scripts[9]),
        "--p11-script", str(scripts[11]),
        "--p13-script", str(scripts[13]),
        "--python-executable", str(Path(sys.executable).resolve()),
        "--p09-n-bootstrap", str(int(n_bootstrap)),
        "--p11-n-bootstrap", str(int(n_bootstrap)),
        "--p13-n-bootstrap", str(int(n_bootstrap)),
        "--run-all",
    ]


def build_command(
    step: int,
    scripts: Dict[int, Path],
    root: Path,
    inputs: Dict[str, Path],
    dataset_label: str,
    detectability_alpha: Optional[float],
    dt_values: Sequence[int],
    min_dt_points: int,
    n_bootstrap: int,
    seed: int,
    fresh: bool,
    resume: bool,
) -> List[str]:
    if step == 9:
        return build_step9_command(
            scripts[9], root, inputs, dataset_label, n_bootstrap, seed, fresh
        )
    if step == 11:
        return build_step11_command(
            scripts[11], root, inputs, dataset_label, dt_values, n_bootstrap, seed
        )
    if step == 13:
        return build_step13_command(
            scripts[13], root, inputs, dataset_label, dt_values,
            min_dt_points, n_bootstrap, seed
        )
    if step == 14:
        return build_step14_command(
            scripts[14], root, inputs, dataset_label, dt_values, n_bootstrap, seed
        )
    if step == 15:
        if detectability_alpha is None:
            raise RuntimeError("Step 15 requires an upstream operational alpha.")
        return build_step15_command(
            scripts[15], scripts, root, inputs, dataset_label,
            detectability_alpha, dt_values, min_dt_points,
            n_bootstrap, seed, fresh, resume
        )
    if step == 16:
        return build_step16_command(
            scripts[16], scripts, root, inputs, dataset_label,
            dt_values, n_bootstrap, seed
        )
    raise ValueError(step)


# =============================================================================
# Logging and display
# =============================================================================


def run_and_tee(
    command: Sequence[str],
    log_path: Path,
    cwd: Path,
    env: Dict[str, str],
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("\n$ " + shell_join(command))
    print("[log] {}".format(log_path))
    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        process = subprocess.Popen(
            [str(value) for value in command],
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log.write(line)
        return int(process.wait())


def write_manifest(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = [
        "step", "script", "status", "started_at", "finished_at",
        "returncode", "log_path", "command", "message",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def print_preflight(
    selected_steps: Sequence[int],
    root: Path,
    dataset_label: str,
    inputs: Dict[str, Path],
    detectability_alpha: Optional[float],
    dt_values: Sequence[int],
    min_dt_points: int,
    n_bootstrap: int,
    seed: int,
) -> None:
    print("\n" + "=" * 84)
    print("ClonoDynamics longitudinal-dynamics orchestrator — {}".format(
        steps_label(selected_steps)
    ))
    print("=" * 84)
    print("Dataset type            : longitudinal")
    print("Dataset label           : {}".format(dataset_label))
    print("Results root            : {}".format(root))
    print("Selected Steps          : {}".format(
        ",".join(str(v) for v in selected_steps)
    ))
    print("Cross-dataset Steps     : 10 and 12 excluded")
    print("Temporal lags           : {}".format(dt_csv(dt_values)))
    print("Bootstrap/permutation n : {}".format(n_bootstrap))
    print("Random seed             : {}".format(seed))
    if 13 in selected_steps or 15 in selected_steps or 16 in selected_steps:
        print("Min dt points           : {}".format(min_dt_points))
    if 15 in selected_steps:
        print("Step-15 fixed alpha     : {}".format(detectability_alpha))
    if 16 in selected_steps:
        print("Step-16 alpha sweep     : {}".format(alpha_text(P16_ALPHAS)))
        print("Step-16 reference alpha : {:g}".format(P16_REFERENCE_ALPHA))

    print("\nResolved / effective inputs:")
    for key in ("step2_dir", "trajectories", "transitions", "step11_dir", "step13_dir"):
        if key in inputs:
            print("  {:18s}: {}".format(key, inputs[key]))

    print("\nDependency notes:")
    print("  Step 9  <- Step 4 trajectories")
    print("  Step 11 <- Step 5 transitions")
    print("  Step 13 <- Step 11 output")
    print("  Step 14 <- Step 5 + Step 11 + Step 13")
    print("  Step 15 <- Step 4 + Step 5 + current scripts 9/11/13")
    print("  Step 16 <- Step 2 + Step 4 + Step 5 + current scripts 9/11/13")
    print("=" * 84)


# =============================================================================
# CLI and main
# =============================================================================


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Dependency-aware orchestrator for genuine-longitudinal "
            "ClonoDynamics Steps 9,11,13,14,15,16."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--workflow",
        default=None,
        help=(
            "complete, core, sensitivity, 9, 11, 13, 14, 15, 16, "
            "or a comma-separated custom subset."
        ),
    )
    parser.add_argument("--results-root", type=Path, default=None)
    parser.add_argument("--code-dir", type=Path, default=None)
    parser.add_argument("--dataset-label", default=None)

    parser.add_argument("--step2-dir", type=Path, default=None)
    parser.add_argument("--trajectories", type=Path, default=None)
    parser.add_argument("--transitions", type=Path, default=None)
    parser.add_argument("--step11-dir", type=Path, default=None)
    parser.add_argument("--step13-dir", type=Path, default=None)

    parser.add_argument(
        "--detectability-alpha",
        type=float,
        default=None,
        help=(
            "Upstream alpha already represented by Step-15 labels. Step 15 "
            "records but does not re-threshold it."
        ),
    )
    parser.add_argument(
        "--dt-values",
        default="1,2,3,4,5",
        help="Comma-separated longitudinal lags.",
    )
    parser.add_argument("--min-dt-points", type=int, default=5)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=123)

    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--yes", action="store_true")
    parser.add_argument("--non-interactive", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    print("ClonoDynamics longitudinal-dynamics orchestrator | {}".format(VERSION))
    args = parse_args(argv)
    interactive = not bool(args.non_interactive)

    if args.workflow:
        selected_steps = normalize_steps(workflow_from_name(args.workflow))
    elif interactive:
        selected_steps = prompt_workflow()
    else:
        raise ValueError("--workflow is required in non-interactive mode.")

    if args.results_root is not None:
        root = normalize_results_root(args.results_root)
    elif interactive:
        root = prompt_results_root()
    else:
        raise ValueError("--results-root is required in non-interactive mode.")

    upstream = upstream_config(root)
    provenance = validate_longitudinal_provenance(root, upstream)
    previous = previous_longitudinal_config(root)
    dataset_label = resolve_dataset_label(
        args.dataset_label, previous, interactive
    )

    dt_values = parse_dt_values(args.dt_values)
    if int(args.min_dt_points) < 2:
        raise ValueError("--min-dt-points must be >=2")
    if int(args.min_dt_points) > len(dt_values) and (
        13 in selected_steps or 15 in selected_steps or 16 in selected_steps
    ):
        raise ValueError(
            "--min-dt-points cannot exceed the number of selected dt values."
        )
    if int(args.n_bootstrap) < 1:
        raise ValueError("--n-bootstrap must be >=1")

    step15_alpha = resolve_step15_alpha(
        selected_steps,
        args.detectability_alpha,
        inherited_alpha(upstream),
        interactive,
    )

    code_dir = resolve_code_dir(args.code_dir)
    scripts = validate_scripts(code_dir, selected_steps)
    resolved = resolve_preexisting_inputs(
        root,
        selected_steps,
        args.trajectories,
        args.transitions,
        args.step2_dir,
        args.step11_dir,
        args.step13_dir,
        interactive,
    )
    inputs = effective_inputs(root, selected_steps, resolved)

    print_preflight(
        selected_steps,
        root,
        dataset_label,
        inputs,
        step15_alpha,
        dt_values,
        int(args.min_dt_points),
        int(args.n_bootstrap),
        int(args.seed),
    )

    mode = choose_run_mode(
        root,
        selected_steps,
        bool(args.restart),
        bool(args.resume),
        interactive,
    )
    if mode == "abort":
        print("Aborted.")
        return 0
    if mode == "fresh":
        clear_selected_outputs(root, selected_steps)

    (root / "logs").mkdir(parents=True, exist_ok=True)
    run_dir = root / "orchestrator_runs_longitudinal_dynamics"
    run_dir.mkdir(parents=True, exist_ok=True)
    inputs = effective_inputs(root, selected_steps, resolved)

    commands: Dict[int, List[str]] = {}
    for step in selected_steps:
        commands[step] = build_command(
            step,
            scripts,
            root,
            inputs,
            dataset_label,
            step15_alpha,
            dt_values,
            int(args.min_dt_points),
            int(args.n_bootstrap),
            int(args.seed),
            mode == "fresh",
            mode == "resume",
        )

    config = {
        "orchestrator": Path(__file__).name,
        "orchestrator_version": VERSION,
        "generated_at": now_iso(),
        "dataset_type": "longitudinal",
        "dataset_label": dataset_label,
        "selected_steps": list(selected_steps),
        "results_root": str(root),
        "code_dir": str(code_dir),
        "dt_values": list(dt_values),
        "min_dt_points": int(args.min_dt_points),
        "n_bootstrap": int(args.n_bootstrap),
        "seed": int(args.seed),
        "step15_fixed_alpha": step15_alpha,
        "step16_alphas": list(P16_ALPHAS),
        "step16_reference_alpha": P16_REFERENCE_ALPHA,
        "step16_latent_states_refitted": False,
        "steps_1_6_provenance": provenance,
        "resolved_preexisting_inputs": {
            key: str(value) for key, value in resolved.items()
        },
        "effective_inputs": {
            key: str(value) for key, value in inputs.items()
        },
        "script_sha256": {
            str(step): sha256_file(path) for step, path in scripts.items()
        },
        "dependency_graph": {
            "step9": ["Step-4 trajectories"],
            "step11": ["Step-5 transitions"],
            "step13": ["Step-11 output"],
            "step14": ["Step-5 transitions", "Step-11 output", "Step-13 output"],
            "step15": [
                "Step-4 trajectories", "Step-5 transitions",
                "current Step-9/11/13 scripts",
            ],
            "step16": [
                "Step-2 empirical endpoint p-values",
                "Step-4 trajectories",
                "Step-5 transitions",
                "current Step-9/11/13 scripts",
            ],
            "step10_managed_here": False,
            "step12_managed_here": False,
        },
        "primary_settings": {
            "step9": {"dt": 1, "n_bins": 20, "q_low": 0.005, "q_high": 0.995},
            "step11": {"mode": "TT", "conditioning": "xstar", "n_bins": 30},
            "step13": {"metrics": ["var_dx", "msd"], "bin_scope": "core"},
            "step15": {"forward": "TT_vs_T0PLUS", "fluctuation": "TT_vs_NON_FF"},
            "step16": {
                "alphas": list(P16_ALPHAS),
                "reference_alpha": P16_REFERENCE_ALPHA,
                "reuses": ["P09", "P11", "P13"],
                "run_all": True,
            },
        },
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
    }
    config_text = json.dumps(config, indent=2, ensure_ascii=False)
    config_path = root / "00_orchestrator_longitudinal_dynamics_config.json"
    config_path.write_text(config_text, encoding="utf-8")
    run_config = run_dir / (
        "{}_{}_config.json".format(
            timestamp_slug(), "-".join(str(v) for v in selected_steps)
        )
    )
    run_config.write_text(config_text, encoding="utf-8")

    print("\nExecution plan:")
    for step in selected_steps:
        state = "SKIP (complete)" if mode == "resume" and step_complete(root, step) else "RUN"
        print("\nStep {}: {}".format(step, state))
        print("  " + shell_join(commands[step]))

    if args.dry_run:
        print("\n[DRY RUN] Dependencies resolved and commands built. No Step was executed.")
        return 0

    if not args.yes:
        answer = input("\nRun {} now? [Y/n]: ".format(steps_label(selected_steps))).strip().lower()
        if answer not in {"", "y", "yes", "s", "si", "sì"}:
            print("Aborted.")
            return 0

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    old_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(code_dir) if not old_pythonpath else str(code_dir) + os.pathsep + old_pythonpath

    manifest_path = root / "00_orchestrator_longitudinal_dynamics_manifest.csv"
    run_manifest_path = run_dir / (
        "{}_{}_manifest.csv".format(
            timestamp_slug(), "-".join(str(v) for v in selected_steps)
        )
    )
    rows: List[Dict[str, object]] = []

    for step in selected_steps:
        script_name = SCRIPT_NAMES[step]
        log_path = root / "logs" / "step{}.log".format(step)

        if mode == "resume" and step_complete(root, step):
            row = {
                "step": step,
                "script": script_name,
                "status": "skipped_complete",
                "started_at": now_iso(),
                "finished_at": now_iso(),
                "returncode": 0,
                "log_path": str(log_path),
                "command": shell_join(commands[step]),
                "message": "Required completion markers already present.",
            }
            rows.append(row)
            write_manifest(manifest_path, rows)
            write_manifest(run_manifest_path, rows)
            print("\n[SKIP] Step {}: complete.".format(step))
            continue

        # P16 is resumed only at public-Step granularity. Remove an incomplete
        # P16 tree before launch so all alpha-specific internal analyses are
        # regenerated coherently.
        if step == 16 and mode == "resume" and outdir(root, 16).exists():
            print("\n[resume] Step 16 is incomplete; restarting its output tree.")
            clear_step_output(root, 16)
            commands[16] = build_command(
                16, scripts, root, inputs, dataset_label, step15_alpha,
                dt_values, int(args.min_dt_points), int(args.n_bootstrap),
                int(args.seed), True, False,
            )

        if step == 13 and 11 in selected_steps:
            inputs["step11_dir"] = validate_step_dir(
                outdir(root, 11),
                "Current-run Step-11 output",
                ("00_analysis_metadata.json", "10_binned_dynamics_long.parquet"),
            )
            commands[13] = build_command(
                13, scripts, root, inputs, dataset_label, step15_alpha,
                dt_values, int(args.min_dt_points), int(args.n_bootstrap),
                int(args.seed), mode == "fresh", mode == "resume",
            )

        if step == 14:
            if 11 in selected_steps:
                inputs["step11_dir"] = validate_step_dir(
                    outdir(root, 11),
                    "Current-run Step-11 output",
                    ("00_analysis_metadata.json", "10_binned_dynamics_long.parquet"),
                )
            if 13 in selected_steps:
                inputs["step13_dir"] = validate_step_dir(
                    outdir(root, 13),
                    "Current-run Step-13 output",
                    ("00_run_config.json", "07_interpretation_aid.csv"),
                )
            commands[14] = build_command(
                14, scripts, root, inputs, dataset_label, step15_alpha,
                dt_values, int(args.min_dt_points), int(args.n_bootstrap),
                int(args.seed), mode == "fresh", mode == "resume",
            )

        outdir(root, step).mkdir(parents=True, exist_ok=True)
        started = now_iso()
        print("\n" + "=" * 84)
        print("STEP {} — {}".format(step, script_name))
        print("=" * 84)
        returncode = run_and_tee(commands[step], log_path, code_dir, env)
        finished = now_iso()
        status = "completed" if returncode == 0 else "failed"
        message = ""

        if returncode == 0 and not step_complete(root, step):
            status = "failed_output_check"
            missing = [
                str(path)
                for path in required_outputs(root, step)
                if not path.exists() or (path.is_file() and path.stat().st_size <= 0)
            ]
            message = "Process exited with code 0 but required markers are missing/empty: " + "; ".join(missing)

        row = {
            "step": step,
            "script": script_name,
            "status": status,
            "started_at": started,
            "finished_at": finished,
            "returncode": returncode,
            "log_path": str(log_path),
            "command": shell_join(commands[step]),
            "message": message,
        }
        rows.append(row)
        write_manifest(manifest_path, rows)
        write_manifest(run_manifest_path, rows)

        if status != "completed":
            print("\n[FAIL] Step {} stopped the workflow.".format(step))
            print("See log:\n  {}".format(log_path))
            if message:
                print(message)
            return returncode if returncode != 0 else 2
        print("[OK] Step {} completed and required outputs were verified.".format(step))

    print("\n" + "=" * 84)
    print("ClonoDynamics longitudinal-dynamics workflow completed successfully.")
    print("Selected : {}".format(steps_label(selected_steps)))
    print("Results  : {}".format(root))
    print("Config   : {}".format(config_path))
    print("Manifest : {}".format(manifest_path))
    print("=" * 84)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, NotADirectoryError, RuntimeError) as exc:
        print("\n[ERROR] {}".format(exc), file=sys.stderr)
        raise SystemExit(2)
