#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clonodynamics_orchestrator_steps_1_6.py
=======================================

Interactive, dependency-aware orchestrator for ClonoDynamics Steps 1-6.

Version
-------
v2 adds selective execution:

    - all Steps 1-6;
    - one single Step;
    - any contiguous Step range, e.g. 2-4, 3-6, 5-6.

The orchestrator resolves the inputs required by the FIRST selected Step and
then uses outputs generated inside the selected range whenever they become
available.

Examples
--------
Interactive:

    python3 clonodynamics_orchestrator_steps_1_6.py

Single Step:

    python3 clonodynamics_orchestrator_steps_1_6.py --step 6

Contiguous range:

    python3 clonodynamics_orchestrator_steps_1_6.py \
        --start-step 3 \
        --end-step 5

Full pipeline:

    python3 clonodynamics_orchestrator_steps_1_6.py \
        --start-step 1 \
        --end-step 6


Pipeline dependency graph
-------------------------
The numbered workflow is not a simple linear chain:

    RAW DATA
      |
      +----------------------> Step 1
      |
      +--> Step 2
             |
             +---------------> Step 3
             |
             +---------------> Step 4
                                |
                                +-------> Step 5 -------> Step 6
             |                             ^
             |                             |
             +-- full posterior NPZ -------+

Exact dependencies:

    Step 1:
        raw repertoire directory

    Step 2:
        raw repertoire directory
        + noiseK_latent.py

    Step 3:
        Step-2 per_clone_latent_subject.parquet

    Step 4:
        Step-2 per_clone_latent_subject.parquet

    Step 5:
        Step-4 latent_trajectories_long.parquet
        +
        Step-2 per_clone_latent_logP/*.posterior.npz

    Step 6:
        Step-5 latent_transitions.parquet


Selective-execution policy
--------------------------
If the selected range starts with Step 1 or Step 2, the orchestrator asks for
the raw dataset path because those Steps consume raw repertoire files.

If the selected range starts at Step 3 or later, the orchestrator asks for an
existing ClonoDynamics result root and searches for the required upstream
artifacts.

Canonical upstream paths are tried first. If a canonical path is absent, the
result tree is searched recursively by exact filename/directory name. If the
search is ambiguous or empty, interactive mode asks the user to choose/provide
the required path.

Examples:

    Step 3 only
        needs existing:
            Step 2 / per_clone_latent_subject.parquet

    Step 4-6
        starts from existing:
            Step 2 / per_clone_latent_subject.parquet

        Step 4 then generates:
            latent_trajectories_long.parquet

        Step 5 uses:
            current Step-4 trajectories
            + existing Step-2 posterior NPZ directory

        Step 6 uses:
            current Step-5 transitions

    Step 5 only
        needs existing:
            Step-4 latent_trajectories_long.parquet
            Step-2 per_clone_latent_logP/

    Step 6 only
        needs existing:
            Step-5 latent_transitions.parquet


Dataset type
------------
The supported dataset labels are:

    longitudinal
    pseudo-longitudinal.

Steps 1-6 use the same computational algorithms for both dataset types.
Dataset type is used for provenance and default result-tree organization.


Step-2 observability alpha
--------------------------
The orchestrator asks for alpha only when Step 2 is part of the selected
execution range.

This alpha controls the empirical Step-2 observability call:

    observable = p_value < alpha.

It is not a threshold on p_detect_state.


Default result tree
-------------------
When raw data are used and --results-root is not supplied:

    /path/to/DATASET

produces:

    /path/to/DATASET_clonodynamics_results/<dataset-type>/

When execution begins at Step 3 or later, an existing result root must be
supplied or selected interactively.


Result directories
------------------
    1-repertoire_characterization/
    2-noise_aware_latent_inference/
    3-latent_posterior_quality_characterization/
    4-latent_trajectory_construction/
    5-latent_transition_construction/
    6-latent_transition_dataset_characterization/


Upstream-input discovery
------------------------
The following CLI overrides are available when canonical upstream outputs have
been moved or stored elsewhere:

    --step2-per-clone
    --step2-posterior-dir
    --step4-trajectories
    --step5-transitions
    --max-dt.

Otherwise, the orchestrator first checks the canonical location and then
searches recursively inside --results-root.


Fresh versus resume
-------------------
Only SELECTED Step outputs are affected.

Fresh:
    deletes only the selected Step directories and selected Step logs.
    Upstream dependency directories outside the selected range are preserved.

Resume:
    skips a selected Step only when its required outputs are already complete.
    Upstream inputs are still validated before execution starts.

This makes partial reruns safe, e.g. a fresh Step 5 rerun does not delete
Step-2 posterior files or Step-4 trajectories.


Fail-fast
---------
Each Step is executed as a subprocess with:

    the same sys.executable used for the orchestrator;

    the code directory prepended to PYTHONPATH;

    stdout/stderr streamed to terminal and Step-specific log.

A non-zero return code or missing required output stops the selected range.


Reproducibility
---------------
Each run records:

    00_orchestrator_config.json
    00_orchestrator_step_manifest.csv

and an immutable timestamped copy under:

    orchestrator_runs/

The configuration records:

    selected Step range;
    dataset type;
    raw-data path when applicable;
    observability alpha when applicable;
    resolved upstream inputs;
    derived temporal span;
    script SHA-256 hashes;
    commands and execution parameters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


ORCHESTRATOR_VERSION = "v2-selective-steps1-6-2026-08-24"

COMMON_PATTERN = (
    r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])"
    r"(?:\.(?:tsv|csv|parquet|pq|feather|arrow))?$"
)

STEP_SCRIPTS = {
    1: "1-repertoire_characterization.py",
    2: "2-noise_aware_latent_inference.py",
    3: "3-latent_posterior_quality_characterization.py",
    4: "4-latent_trajectory_construction.py",
    5: "5-latent_transition_construction.py",
    6: "6-latent_transition_dataset_characterization.py",
}

NOISE_MODULE = "noiseK_latent.py"

STEP_DIR_NAMES = {
    1: "1-repertoire_characterization",
    2: "2-noise_aware_latent_inference",
    3: "3-latent_posterior_quality_characterization",
    4: "4-latent_trajectory_construction",
    5: "5-latent_transition_construction",
    6: "6-latent_transition_dataset_characterization",
}


# =============================================================================
# Basic utilities
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
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def slug_dataset_type(value: str) -> str:
    value = str(value).strip().lower().replace("_", "-")
    aliases = {
        "1": "longitudinal",
        "l": "longitudinal",
        "long": "longitudinal",
        "longitudinal": "longitudinal",
        "2": "pseudo-longitudinal",
        "p": "pseudo-longitudinal",
        "pseudo": "pseudo-longitudinal",
        "pseudo-longitudinal": "pseudo-longitudinal",
        "pseudolongitudinal": "pseudo-longitudinal",
    }
    if value not in aliases:
        raise ValueError(
            "dataset type must be 'longitudinal' or 'pseudo-longitudinal'"
        )
    return aliases[value]


def step_range(start_step: int, end_step: int) -> List[int]:
    if not (1 <= int(start_step) <= 6):
        raise ValueError("start step must be between 1 and 6")
    if not (1 <= int(end_step) <= 6):
        raise ValueError("end step must be between 1 and 6")
    if int(end_step) < int(start_step):
        raise ValueError("end step must be >= start step")
    return list(range(int(start_step), int(end_step) + 1))


def step_paths(results_root: Path) -> Dict[int, Path]:
    return {
        step: results_root / dirname
        for step, dirname in STEP_DIR_NAMES.items()
    }


def default_results_root(data_path: Path, dataset_type: str) -> Path:
    return (
        data_path.parent
        / f"{data_path.name}_clonodynamics_results"
        / dataset_type
    )


# =============================================================================
# Interactive prompts
# =============================================================================


def prompt_execution_scope() -> Tuple[int, int]:
    print("\nExecution scope:")
    print("  [1] all Steps 1-6")
    print("  [2] one Step")
    print("  [3] contiguous Step range")

    while True:
        raw = input("Selection [1]: ").strip().lower()
        if raw in {"", "1", "all", "a"}:
            return 1, 6

        if raw in {"2", "one", "single", "s"}:
            while True:
                value = input("Step to run [1-6]: ").strip()
                try:
                    step = int(value)
                except ValueError:
                    print("Enter an integer between 1 and 6.")
                    continue
                if 1 <= step <= 6:
                    return step, step
                print("Enter an integer between 1 and 6.")

        if raw in {"3", "range", "r"}:
            while True:
                first = input("First Step [1-6]: ").strip()
                last = input("Last Step  [1-6]: ").strip()
                try:
                    first_i = int(first)
                    last_i = int(last)
                    step_range(first_i, last_i)
                    return first_i, last_i
                except Exception as exc:
                    print(f"Invalid range: {exc}")

        print("Please enter 1, 2, or 3.")


def prompt_dataset_type(default: Optional[str] = None) -> str:
    print("\nDataset type:")
    print("  [1] longitudinal")
    print("  [2] pseudo-longitudinal")
    suffix = f" [{default}]" if default else ""
    while True:
        raw = input(f"Selection{suffix}: ").strip()
        if not raw and default:
            return slug_dataset_type(default)
        try:
            return slug_dataset_type(raw)
        except ValueError:
            print("Please enter 1/2, longitudinal, or pseudo-longitudinal.")


def prompt_existing_directory(label: str) -> Path:
    while True:
        raw = input(f"\n{label}: ").strip()
        if not raw:
            print("A directory path is required.")
            continue
        path = Path(raw).expanduser()
        try:
            path = path.resolve(strict=True)
        except FileNotFoundError:
            print(f"Path not found: {path}")
            continue
        if not path.is_dir():
            print(f"Not a directory: {path}")
            continue
        return path


def prompt_data_path() -> Path:
    return prompt_existing_directory("Raw dataset path")


def prompt_alpha(default: float = 0.05) -> float:
    while True:
        raw = input(
            "\nObservability alpha "
            f"(Step 2 empirical logP threshold) [{default:g}]: "
        ).strip()

        if not raw:
            return float(default)

        try:
            value = float(raw)
        except ValueError:
            print("Alpha must be numeric.")
            continue

        if not (0.0 < value < 1.0):
            print("Alpha must satisfy 0 < alpha < 1.")
            continue

        return value


def prompt_path_choice(
    *,
    label: str,
    candidates: Sequence[Path],
    expect_dir: bool,
) -> Path:
    print(f"\nRequired upstream input: {label}")

    if candidates:
        print("Candidates found:")
        for i, path in enumerate(candidates, start=1):
            print(f"  [{i}] {path}")
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
        raw = input(f"Path for {label}: ").strip()
        if not raw:
            print("A path is required.")
            continue

        path = Path(raw).expanduser()
        try:
            path = path.resolve(strict=True)
        except FileNotFoundError:
            print(f"Path not found: {path}")
            continue

        if expect_dir and not path.is_dir():
            print(f"Expected a directory: {path}")
            continue
        if not expect_dir and not path.is_file():
            print(f"Expected a file: {path}")
            continue
        return path


# =============================================================================
# Raw dataset audit
# =============================================================================


@dataclass(frozen=True)
class RawDatasetAudit:
    n_matching_files: int
    n_subjects: int
    n_unique_subject_time_pairs: int
    n_complete_pairs: int
    n_incomplete_pairs: int
    complete_pair_ids: Tuple[str, ...]
    incomplete_pair_ids: Tuple[str, ...]
    max_dt_complete_pairs: int
    file_separator: str
    suffixes: Tuple[str, ...]


def audit_raw_dataset(data_path: Path, pattern_str: str) -> RawDatasetAudit:
    pattern = re.compile(pattern_str)

    keyed: Dict[Tuple[int, int, int], Path] = {}
    by_pair: Dict[Tuple[int, int], set] = {}
    suffixes: set = set()

    for path in sorted(data_path.iterdir()):
        if not path.is_file():
            continue

        match = pattern.match(path.name)
        if match is None:
            continue

        subject = int(match.group("subject"))
        time = int(match.group("time"))
        replica = int(match.group("replica"))
        key = (subject, time, replica)

        if key in keyed:
            raise ValueError(
                "Duplicate raw repertoire key detected for "
                f"subject={subject}, time={time}, replica={replica}:\n"
                f"  {keyed[key]}\n"
                f"  {path}"
            )

        keyed[key] = path
        by_pair.setdefault((subject, time), set()).add(replica)
        suffixes.add(path.suffix.lower())

    if not keyed:
        raise ValueError(
            "No repertoire files matched the expected filename pattern:\n"
            f"  {pattern_str}\n"
            f"in:\n  {data_path}"
        )

    complete = sorted(
        pair for pair, reps in by_pair.items()
        if {1, 2}.issubset(reps)
    )
    incomplete = sorted(
        pair for pair, reps in by_pair.items()
        if not {1, 2}.issubset(reps)
    )

    by_subject_complete: Dict[int, List[int]] = {}
    for subject, time in complete:
        by_subject_complete.setdefault(subject, []).append(time)

    spans: List[int] = []
    for times in by_subject_complete.values():
        if len(times) >= 2:
            spans.append(max(times) - min(times))

    max_dt = max(spans) if spans else 0

    csv_like = ".csv" in suffixes
    tab_like = bool({".tsv", ".txt"} & suffixes)

    if csv_like and tab_like:
        raise ValueError(
            "Mixed comma- and tab-delimited raw text files were detected. "
            "Steps 1 and 2 use one shared --file-sep. "
            "Standardize the raw text delimiter before running."
        )

    file_sep = "," if csv_like else "\t"

    return RawDatasetAudit(
        n_matching_files=len(keyed),
        n_subjects=len({s for s, _t in by_pair}),
        n_unique_subject_time_pairs=len(by_pair),
        n_complete_pairs=len(complete),
        n_incomplete_pairs=len(incomplete),
        complete_pair_ids=tuple(f"{s}_{t}" for s, t in complete),
        incomplete_pair_ids=tuple(f"{s}_{t}" for s, t in incomplete),
        max_dt_complete_pairs=int(max_dt),
        file_separator=file_sep,
        suffixes=tuple(sorted(suffixes)),
    )


# =============================================================================
# Script discovery
# =============================================================================


def ensure_selected_scripts(
    code_dir: Path,
    selected_steps: Sequence[int],
) -> Dict[str, Path]:
    resolved: Dict[str, Path] = {}
    missing: List[Path] = []

    for step in selected_steps:
        name = STEP_SCRIPTS[int(step)]
        path = code_dir / name
        resolved[f"step{step}"] = path
        if not path.is_file():
            missing.append(path)

    if 2 in selected_steps:
        noise_path = code_dir / NOISE_MODULE
        resolved["noise_model"] = noise_path
        if not noise_path.is_file():
            missing.append(noise_path)

    if missing:
        details = "\n".join(f"  - {p}" for p in missing)
        raise FileNotFoundError(
            "Required files for the SELECTED Steps are missing from --code-dir:\n"
            + details
        )

    return resolved


# =============================================================================
# Upstream artifact discovery
# =============================================================================


def canonical_upstream_paths(results_root: Path) -> Dict[str, Path]:
    s = step_paths(results_root)
    return {
        "step2_per_clone": s[2] / "per_clone_latent_subject.parquet",
        "step2_posterior_dir": s[2] / "per_clone_latent_logP",
        "step4_trajectories": s[4] / "latent_trajectories_long.parquet",
        "step5_transitions": s[5] / "latent_transitions.parquet",
    }


def validate_resolved_input(
    path: Path,
    *,
    label: str,
    expect_dir: bool,
) -> Path:
    path = path.expanduser().resolve(strict=True)

    if expect_dir and not path.is_dir():
        raise NotADirectoryError(f"{label}: expected directory: {path}")
    if not expect_dir and not path.is_file():
        raise FileNotFoundError(f"{label}: expected file: {path}")

    if expect_dir and label == "Step-2 full posterior directory":
        if not any(path.glob("*.posterior.npz")):
            raise FileNotFoundError(
                "Step-2 posterior directory contains no *.posterior.npz files:\n"
                f"  {path}"
            )

    return path


def recursive_exact_candidates(
    results_root: Path,
    *,
    basename: str,
    expect_dir: bool,
) -> List[Path]:
    candidates: List[Path] = []

    if not results_root.exists():
        return candidates

    for path in results_root.rglob(basename):
        if expect_dir and path.is_dir():
            candidates.append(path.resolve())
        elif not expect_dir and path.is_file():
            candidates.append(path.resolve())

    return sorted(set(candidates), key=lambda p: (len(p.parts), str(p)))


def resolve_existing_artifact(
    *,
    label: str,
    canonical: Path,
    basename: str,
    expect_dir: bool,
    override: Optional[Path],
    results_root: Path,
    interactive: bool,
) -> Path:
    if override is not None:
        return validate_resolved_input(
            override,
            label=label,
            expect_dir=expect_dir,
        )

    if canonical.exists():
        return validate_resolved_input(
            canonical,
            label=label,
            expect_dir=expect_dir,
        )

    candidates = recursive_exact_candidates(
        results_root,
        basename=basename,
        expect_dir=expect_dir,
    )

    valid: List[Path] = []
    for candidate in candidates:
        try:
            valid.append(
                validate_resolved_input(
                    candidate,
                    label=label,
                    expect_dir=expect_dir,
                )
            )
        except Exception:
            continue

    if len(valid) == 1:
        print(f"[input] {label}: {valid[0]}")
        return valid[0]

    if not interactive:
        if not valid:
            raise FileNotFoundError(
                f"Could not find required upstream input: {label}\n"
                f"Canonical path checked:\n  {canonical}"
            )
        raise RuntimeError(
            f"Multiple candidates found for {label}; provide an explicit CLI override:\n"
            + "\n".join(f"  - {p}" for p in valid)
        )

    chosen = prompt_path_choice(
        label=label,
        candidates=valid,
        expect_dir=expect_dir,
    )
    return validate_resolved_input(
        chosen,
        label=label,
        expect_dir=expect_dir,
    )


def resolve_upstream_inputs(
    *,
    selected_steps: Sequence[int],
    results_root: Path,
    step2_per_clone_override: Optional[Path],
    step2_posterior_override: Optional[Path],
    step4_trajectories_override: Optional[Path],
    step5_transitions_override: Optional[Path],
    interactive: bool,
) -> Dict[str, Path]:
    """
    Resolve only dependencies that must PRE-EXIST before the selected range starts.

    Dependencies generated earlier inside the selected range are not searched for.
    """
    selected = set(int(x) for x in selected_steps)
    canonical = canonical_upstream_paths(results_root)
    resolved: Dict[str, Path] = {}

    # Step 3 and Step 4 require Step-2 aggregate input.
    if (3 in selected or 4 in selected) and 2 not in selected:
        resolved["step2_per_clone"] = resolve_existing_artifact(
            label="Step-2 per_clone_latent_subject.parquet",
            canonical=canonical["step2_per_clone"],
            basename="per_clone_latent_subject.parquet",
            expect_dir=False,
            override=step2_per_clone_override,
            results_root=results_root,
            interactive=interactive,
        )

    # Step 5 requires Step-2 posterior NPZ directory unless Step 2 runs now.
    if 5 in selected and 2 not in selected:
        resolved["step2_posterior_dir"] = resolve_existing_artifact(
            label="Step-2 full posterior directory",
            canonical=canonical["step2_posterior_dir"],
            basename="per_clone_latent_logP",
            expect_dir=True,
            override=step2_posterior_override,
            results_root=results_root,
            interactive=interactive,
        )

    # Step 5 requires Step-4 trajectories unless Step 4 runs now.
    if 5 in selected and 4 not in selected:
        resolved["step4_trajectories"] = resolve_existing_artifact(
            label="Step-4 latent_trajectories_long.parquet",
            canonical=canonical["step4_trajectories"],
            basename="latent_trajectories_long.parquet",
            expect_dir=False,
            override=step4_trajectories_override,
            results_root=results_root,
            interactive=interactive,
        )

    # Step 6 requires Step-5 transitions unless Step 5 runs now.
    if 6 in selected and 5 not in selected:
        resolved["step5_transitions"] = resolve_existing_artifact(
            label="Step-5 latent_transitions.parquet",
            canonical=canonical["step5_transitions"],
            basename="latent_transitions.parquet",
            expect_dir=False,
            override=step5_transitions_override,
            results_root=results_root,
            interactive=interactive,
        )

    return resolved


def effective_input_paths(
    *,
    selected_steps: Sequence[int],
    results_root: Path,
    resolved_upstream: Dict[str, Path],
) -> Dict[str, Path]:
    """
    Return the paths commands should use.

    If an upstream Step is inside the selected range, its canonical CURRENT-RUN
    output path is used. Otherwise, use the discovered pre-existing artifact.
    """
    selected = set(int(x) for x in selected_steps)
    canonical = canonical_upstream_paths(results_root)

    out: Dict[str, Path] = {}

    if 2 in selected:
        out["step2_per_clone"] = canonical["step2_per_clone"]
        out["step2_posterior_dir"] = canonical["step2_posterior_dir"]
    else:
        if "step2_per_clone" in resolved_upstream:
            out["step2_per_clone"] = resolved_upstream["step2_per_clone"]
        elif canonical["step2_per_clone"].exists():
            out["step2_per_clone"] = canonical["step2_per_clone"]

        if "step2_posterior_dir" in resolved_upstream:
            out["step2_posterior_dir"] = resolved_upstream["step2_posterior_dir"]
        elif canonical["step2_posterior_dir"].exists():
            out["step2_posterior_dir"] = canonical["step2_posterior_dir"]

    if 4 in selected:
        out["step4_trajectories"] = canonical["step4_trajectories"]
    elif "step4_trajectories" in resolved_upstream:
        out["step4_trajectories"] = resolved_upstream["step4_trajectories"]
    elif canonical["step4_trajectories"].exists():
        out["step4_trajectories"] = canonical["step4_trajectories"]

    if 5 in selected:
        out["step5_transitions"] = canonical["step5_transitions"]
    elif "step5_transitions" in resolved_upstream:
        out["step5_transitions"] = resolved_upstream["step5_transitions"]
    elif canonical["step5_transitions"].exists():
        out["step5_transitions"] = canonical["step5_transitions"]

    return out


# =============================================================================
# Infer temporal span when raw data are not part of the selected run
# =============================================================================


def infer_max_dt_from_posterior_dir(path: Path) -> int:
    """
    Infer the maximum within-subject nominal time span from Step-2 posterior
    filenames.

    Canonical posterior names are:

        latent_per_clone_logP_<subject>_<time>.posterior.npz

    This avoids loading the large Step-2/Step-4 parquet tables merely to
    reconstruct the nominal temporal span.
    """
    pattern = re.compile(
        r"^latent_per_clone_logP_(?P<subject>\d+)_(?P<time>\d+)"
        r"\.posterior\.npz$"
    )

    by_subject: Dict[int, List[int]] = {}

    for file_path in sorted(path.glob("*.posterior.npz")):
        match = pattern.match(file_path.name)
        if match is None:
            continue
        subject = int(match.group("subject"))
        time = int(match.group("time"))
        by_subject.setdefault(subject, []).append(time)

    spans = [
        max(times) - min(times)
        for times in by_subject.values()
        if len(set(times)) >= 2
    ]

    if not spans:
        raise ValueError(
            "Could not infer a positive temporal span from canonical Step-2 "
            "posterior filenames in:\n"
            f"  {path}"
        )

    return int(max(spans))


def prompt_max_dt() -> int:
    while True:
        raw = input(
            "\nMaximum nominal dt for Step 5 "
            "(could not be inferred automatically): "
        ).strip()
        try:
            value = int(raw)
        except ValueError:
            print("Enter a positive integer.")
            continue
        if value >= 1:
            return value
        print("Enter a positive integer.")


def resolve_max_dt(
    *,
    selected_steps: Sequence[int],
    explicit_max_dt: Optional[int],
    audit: Optional[RawDatasetAudit],
    inputs: Dict[str, Path],
    existing_config: Optional[Dict[str, object]],
    interactive: bool,
) -> Optional[int]:
    if 5 not in selected_steps:
        return None

    if explicit_max_dt is not None:
        if int(explicit_max_dt) < 1:
            raise ValueError("--max-dt must be >= 1")
        return int(explicit_max_dt)

    if audit is not None and audit.max_dt_complete_pairs >= 1:
        return int(audit.max_dt_complete_pairs)

    if existing_config is not None:
        old = existing_config.get("derived_max_dt")
        try:
            old_i = int(old)
            if old_i >= 1:
                return old_i
        except Exception:
            pass

    posterior_dir = inputs.get("step2_posterior_dir")
    if posterior_dir is not None and posterior_dir.exists():
        try:
            return infer_max_dt_from_posterior_dir(posterior_dir)
        except Exception as exc:
            print(f"[warning] Automatic max-dt inference failed: {exc}")

    if interactive:
        return prompt_max_dt()

    raise RuntimeError(
        "Step 5 is selected but max dt could not be inferred automatically. "
        "Provide --max-dt explicitly."
    )


# =============================================================================
# Output checks / selective cleanup
# =============================================================================


def required_outputs(results_root: Path) -> Dict[int, List[Path]]:
    s = step_paths(results_root)
    return {
        1: [
            s[1] / "repertoire_characteristics.csv",
            s[1] / "repertoire_tail_model_scan.csv",
        ],
        2: [
            s[2] / "latent_params_by_pair.csv",
            s[2] / "per_clone_latent_subject.parquet",
            s[2] / "per_clone_latent_logP",
        ],
        3: [
            s[3] / "clonotype_latent_posterior_qc.parquet",
        ],
        4: [
            s[4] / "latent_trajectories_long.parquet",
        ],
        5: [
            s[5] / "latent_transitions.parquet",
            s[5] / "latent_transition_build_report.md",
        ],
        6: [
            s[6] / "01_transition_dataset_summary.csv",
        ],
    }


def output_is_complete(
    step: int,
    outputs: Dict[int, List[Path]],
) -> bool:
    for path in outputs[int(step)]:
        if not path.exists():
            return False

        if path.is_file() and path.stat().st_size <= 0:
            return False

        if path.is_dir():
            if int(step) == 2 and path.name == "per_clone_latent_logP":
                if not any(path.glob("*.posterior.npz")):
                    return False
            elif not any(path.iterdir()):
                return False

    return True


def selected_output_tree_exists(
    results_root: Path,
    selected_steps: Sequence[int],
) -> bool:
    s = step_paths(results_root)
    for step in selected_steps:
        path = s[int(step)]
        if path.exists() and (
            (path.is_dir() and any(path.iterdir()))
            or path.is_file()
        ):
            return True
    return False


def clear_selected_outputs(
    results_root: Path,
    selected_steps: Sequence[int],
) -> None:
    s = step_paths(results_root)

    for step in selected_steps:
        path = s[int(step)]
        if path.exists():
            shutil.rmtree(path)

        log_path = results_root / "logs" / f"step{step}.log"
        if log_path.exists():
            log_path.unlink()


# =============================================================================
# Existing run metadata
# =============================================================================


def load_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def maybe_load_existing_config(
    results_root: Path,
) -> Optional[Dict[str, object]]:
    path = results_root / "00_orchestrator_config.json"
    if not path.exists():
        return None
    try:
        return load_json(path)
    except Exception as exc:
        print(f"[warning] Could not read existing config {path}: {exc}")
        return None


def infer_dataset_type_from_config(
    config: Optional[Dict[str, object]],
) -> Optional[str]:
    if not config:
        return None
    value = config.get("dataset_type")
    if value is None:
        return None
    try:
        return slug_dataset_type(str(value))
    except Exception:
        return None


def infer_alpha_from_config(
    config: Optional[Dict[str, object]],
) -> Optional[float]:
    if not config:
        return None
    value = config.get("observability_alpha")
    try:
        value_f = float(value)
    except Exception:
        return None
    return value_f if 0.0 < value_f < 1.0 else None


def infer_raw_path_from_config(
    config: Optional[Dict[str, object]],
) -> Optional[Path]:
    if not config:
        return None
    value = config.get("raw_dataset_path")
    if not value:
        return None
    path = Path(str(value)).expanduser()
    if path.exists() and path.is_dir():
        return path.resolve()
    return None


# =============================================================================
# Command construction
# =============================================================================


def build_step_command(
    *,
    step: int,
    code_dir: Path,
    results_root: Path,
    inputs: Dict[str, Path],
    data_path: Optional[Path],
    alpha: Optional[float],
    file_sep: Optional[str],
    max_dt: Optional[int],
    n_jobs: int,
    n_jobs_profile: str,
    posterior_samples: int,
    seed: int,
    restart_step5: bool,
) -> List[str]:
    python = str(Path(sys.executable).resolve())
    s = step_paths(results_root)

    if int(step) == 1:
        if data_path is None or file_sep is None:
            raise RuntimeError("Step 1 requires raw data and file separator.")
        return [
            python, "-u", str(code_dir / STEP_SCRIPTS[1]),
            "--data_dir", str(data_path),
            "--out_dir", str(s[1]),
            "--file_sep", str(file_sep),
            "--pattern", COMMON_PATTERN,
        ]

    if int(step) == 2:
        if data_path is None or file_sep is None or alpha is None:
            raise RuntimeError(
                "Step 2 requires raw data, file separator and observability alpha."
            )
        return [
            python, "-u", str(code_dir / STEP_SCRIPTS[2]),
            "--data-dir", str(data_path),
            "--results-dir", str(s[2]),
            "--file-sep", str(file_sep),
            "--pattern", COMMON_PATTERN,
            "--alpha", repr(float(alpha)),
            "--null", "subject",
            "--tail", "low",
            "--intermediate-format", "parquet",
            "--parquet-compression", "snappy",
            "--n-jobs", str(int(n_jobs)),
            "--n-jobs-profile", str(n_jobs_profile),
        ]

    if int(step) == 3:
        per_clone = inputs.get("step2_per_clone")
        if per_clone is None:
            raise RuntimeError("Step 3 requires Step-2 per-clone aggregate.")
        return [
            python, "-u", str(code_dir / STEP_SCRIPTS[3]),
            "--input", str(per_clone),
            "--output-dir", str(s[3]),
            "--exact-bc-col", "replicate_posterior_bhattacharyya",
            "--require-exact-overlap",
            "--compression", "zstd",
        ]

    if int(step) == 4:
        per_clone = inputs.get("step2_per_clone")
        if per_clone is None:
            raise RuntimeError("Step 4 requires Step-2 per-clone aggregate.")
        return [
            python, "-u", str(code_dir / STEP_SCRIPTS[4]),
            "--per-clone", str(per_clone),
            "--outdir", str(s[4]),
            "--output-format", "parquet",
            "--parquet-compression", "zstd",
            "--freq-col", "f_latent_median",
            "--log-freq-col", "x_latent_median",
            "--observable-col", "observable",
            "--streaming",
        ]

    if int(step) == 5:
        trajectories = inputs.get("step4_trajectories")
        posterior_dir = inputs.get("step2_posterior_dir")

        if trajectories is None:
            raise RuntimeError("Step 5 requires Step-4 trajectories.")
        if posterior_dir is None:
            raise RuntimeError("Step 5 requires Step-2 posterior directory.")
        if max_dt is None or int(max_dt) < 1:
            raise RuntimeError("Step 5 requires a positive max_dt.")

        command = [
            python, "-u", str(code_dir / STEP_SCRIPTS[5]),
            "--trajectories", str(trajectories),
            "--posterior-dir", str(posterior_dir),
            "--out", str(s[5] / "latent_transitions.parquet"),
            "--output-format", "parquet",
            "--parquet-compression", "zstd",
            "--parquet-row-group-size", "100000",
            "--transition-mode", "all",
            "--min-dt", "1",
            "--max-dt", str(int(max_dt)),
            "--n-posterior-samples", str(int(posterior_samples)),
            "--posterior-seed", str(int(seed)),
            "--posterior-chunk-size", "5000",
            "--work-dir", str(s[5] / ".step5_work"),
            "--report-md", str(
                s[5] / "latent_transition_build_report.md"
            ),
        ]

        if restart_step5:
            command.append("--restart")

        return command

    if int(step) == 6:
        transitions = inputs.get("step5_transitions")
        if transitions is None:
            raise RuntimeError("Step 6 requires Step-5 transitions.")
        return [
            python, "-u", str(code_dir / STEP_SCRIPTS[6]),
            "--transitions", str(transitions),
            "--outdir", str(s[6]),
            "--streaming",
            "--write-parquet",
            "--compression", "zstd",
            "--figures",
            "--report-md", str(
                s[6] / "latent_transition_dataset_summary.md"
            ),
        ]

    raise ValueError(f"Unsupported Step: {step}")


# =============================================================================
# Run configuration and logging
# =============================================================================


def build_run_config(
    *,
    start_step: int,
    end_step: int,
    dataset_type: str,
    data_path: Optional[Path],
    results_root: Path,
    code_dir: Path,
    alpha: Optional[float],
    audit: Optional[RawDatasetAudit],
    selected_scripts: Dict[str, Path],
    resolved_upstream: Dict[str, Path],
    effective_inputs: Dict[str, Path],
    max_dt: Optional[int],
    n_jobs: int,
    n_jobs_profile: str,
    posterior_samples: int,
    seed: int,
) -> Dict[str, object]:
    return {
        "orchestrator": Path(__file__).name,
        "orchestrator_version": ORCHESTRATOR_VERSION,
        "generated_at": now_iso(),
        "selected_start_step": int(start_step),
        "selected_end_step": int(end_step),
        "selected_steps": step_range(start_step, end_step),
        "dataset_type": dataset_type,
        "raw_dataset_path": str(data_path) if data_path is not None else None,
        "results_root": str(results_root),
        "code_dir": str(code_dir),
        "observability_alpha": (
            float(alpha) if alpha is not None else None
        ),
        "observability_definition": (
            "Step-2 empirical logP p_value < alpha"
            if alpha is not None
            else "not redefined in this selected run"
        ),
        "step2_null": "subject",
        "step2_tail": "low",
        "common_filename_pattern": COMMON_PATTERN,
        "raw_file_separator": (
            audit.file_separator if audit is not None else None
        ),
        "raw_dataset_audit": (
            asdict(audit) if audit is not None else None
        ),
        "derived_max_dt": (
            int(max_dt) if max_dt is not None else None
        ),
        "resolved_preexisting_upstream_inputs": {
            key: str(value)
            for key, value in resolved_upstream.items()
        },
        "effective_pipeline_inputs": {
            key: str(value)
            for key, value in effective_inputs.items()
        },
        "step2_n_jobs": int(n_jobs),
        "step2_n_jobs_profile": str(n_jobs_profile),
        "step5_n_posterior_samples": int(posterior_samples),
        "random_seed": int(seed),
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "selected_script_sha256": {
            key: sha256_file(path)
            for key, path in selected_scripts.items()
        },
        "pipeline_dependencies": {
            "step1": ["raw_dataset"],
            "step2": ["raw_dataset", "noiseK_latent.py"],
            "step3": ["step2/per_clone_latent_subject.parquet"],
            "step4": ["step2/per_clone_latent_subject.parquet"],
            "step5": [
                "step4/latent_trajectories_long.parquet",
                "step2/per_clone_latent_logP/*.posterior.npz",
            ],
            "step6": ["step5/latent_transitions.parquet"],
        },
    }


def write_manifest(
    path: Path,
    rows: List[Dict[str, object]],
) -> None:
    fieldnames = [
        "step",
        "script",
        "status",
        "started_at",
        "finished_at",
        "returncode",
        "log_path",
        "command",
        "message",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: row.get(key, "") for key in fieldnames}
            )


def run_and_tee(
    command: Sequence[str],
    *,
    log_path: Path,
    cwd: Path,
    env: Dict[str, str],
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n$ " + shell_join(command))
    print(f"[log] {log_path}")

    with log_path.open(
        "w",
        encoding="utf-8",
        buffering=1,
    ) as log:
        process = subprocess.Popen(
            [str(x) for x in command],
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


# =============================================================================
# Existing selected-output mode
# =============================================================================


def choose_selected_run_mode(
    results_root: Path,
    selected_steps: Sequence[int],
    *,
    forced_restart: bool,
    forced_resume: bool,
    interactive: bool,
) -> str:
    if forced_restart and forced_resume:
        raise ValueError("Use only one of --restart or --resume.")

    selected_exist = selected_output_tree_exists(
        results_root,
        selected_steps,
    )

    if not selected_exist:
        return "fresh"

    if forced_restart:
        return "fresh"

    if forced_resume:
        return "resume"

    if not interactive:
        raise RuntimeError(
            "Selected Step outputs already exist. "
            "Use --restart or --resume."
        )

    label = (
        f"Step {selected_steps[0]}"
        if len(selected_steps) == 1
        else f"Steps {selected_steps[0]}-{selected_steps[-1]}"
    )

    print(f"\nExisting outputs detected for selected {label}:")
    print(f"  {results_root}")
    print("  [1] fresh   — delete only selected Step outputs and rerun")
    print("  [2] resume  — skip complete selected Steps and continue")
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
# Display
# =============================================================================


def print_preflight(
    *,
    start_step: int,
    end_step: int,
    dataset_type: str,
    data_path: Optional[Path],
    results_root: Path,
    alpha: Optional[float],
    audit: Optional[RawDatasetAudit],
    resolved_upstream: Dict[str, Path],
    max_dt: Optional[int],
) -> None:
    label = (
        f"Step {start_step}"
        if start_step == end_step
        else f"Steps {start_step}-{end_step}"
    )

    print("\n" + "=" * 76)
    print(f"ClonoDynamics — {label}")
    print("=" * 76)
    print(f"Dataset type          : {dataset_type}")
    print(
        f"Raw dataset           : "
        f"{data_path if data_path is not None else '(not required)'}"
    )
    print(f"Results root          : {results_root}")
    print(
        f"Observability alpha   : "
        f"{alpha:g}" if alpha is not None
        else "Observability alpha   : inherited/not used"
    )

    if audit is not None:
        sep_label = (
            "comma"
            if audit.file_separator == ","
            else "tab"
        )
        print(f"Matching raw files    : {audit.n_matching_files}")
        print(f"Subjects              : {audit.n_subjects}")
        print(f"Subject/time states   : {audit.n_unique_subject_time_pairs}")
        print(f"Complete rep. pairs   : {audit.n_complete_pairs}")
        print(f"Incomplete pairs      : {audit.n_incomplete_pairs}")
        print(f"Raw text separator    : {sep_label}")
        print(
            f"Detected suffixes     : "
            f"{', '.join(audit.suffixes) or '(none)'}"
        )

    if max_dt is not None:
        print(f"Derived max dt        : {max_dt}")

    if resolved_upstream:
        print("\nResolved PRE-EXISTING upstream inputs:")
        for key, path in resolved_upstream.items():
            print(f"  {key:24s}: {path}")

    print("=" * 76)


# =============================================================================
# CLI
# =============================================================================


def parse_args(
    argv: Optional[Sequence[str]] = None,
) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interactive dependency-aware ClonoDynamics orchestrator "
            "for Steps 1-6."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    scope = parser.add_mutually_exclusive_group()
    scope.add_argument(
        "--step",
        type=int,
        default=None,
        help="Run one Step only.",
    )
    scope.add_argument(
        "--start-step",
        type=int,
        default=None,
        help="First Step of a contiguous range.",
    )

    parser.add_argument(
        "--end-step",
        type=int,
        default=None,
        help="Last Step of a contiguous range. Requires --start-step.",
    )

    parser.add_argument(
        "--dataset-type",
        choices=["longitudinal", "pseudo-longitudinal"],
        default=None,
        help=(
            "Dataset provenance label. For Step >=3 it can be inherited "
            "from an existing orchestrator config."
        ),
    )

    parser.add_argument(
        "--data-path",
        type=Path,
        default=None,
        help="Raw dataset directory. Required only when Step 1 or 2 is selected.",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help=(
            "Step-2 empirical-logP observability alpha. "
            "Requested only when Step 2 is selected."
        ),
    )

    parser.add_argument(
        "--results-root",
        type=Path,
        default=None,
        help=(
            "Result root. Required/asked interactively when starting at Step 3+; "
            "otherwise automatically derived from raw dataset when omitted."
        ),
    )

    parser.add_argument(
        "--code-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing selected Step scripts. "
            "Default: directory containing this orchestrator."
        ),
    )

    # Explicit upstream-input overrides.
    parser.add_argument(
        "--step2-per-clone",
        type=Path,
        default=None,
        help="Override Step-2 per_clone_latent_subject.parquet.",
    )
    parser.add_argument(
        "--step2-posterior-dir",
        type=Path,
        default=None,
        help="Override Step-2 per_clone_latent_logP directory.",
    )
    parser.add_argument(
        "--step4-trajectories",
        type=Path,
        default=None,
        help="Override Step-4 latent_trajectories_long.parquet.",
    )
    parser.add_argument(
        "--step5-transitions",
        type=Path,
        default=None,
        help="Override Step-5 latent_transitions.parquet.",
    )

    parser.add_argument(
        "--n-jobs",
        type=int,
        default=0,
        help="Step-2 worker count; <=0 delegates to Step-2 automatic resolution.",
    )
    parser.add_argument(
        "--n-jobs-profile",
        choices=["memory", "balanced", "speed"],
        default="balanced",
        help="Step-2 multiprocessing profile.",
    )
    parser.add_argument(
        "--posterior-samples",
        type=int,
        default=1000,
        help="Step-5 endpoint-posterior Monte Carlo samples per transition.",
    )
    parser.add_argument(
        "--max-dt",
        type=int,
        default=None,
        help=(
            "Override Step-5 maximum nominal dt. Normally inferred from raw "
            "complete pairs, previous config, or Step-2 posterior filenames."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Step-5 posterior propagation seed.",
    )

    parser.add_argument(
        "--restart",
        action="store_true",
        help="Delete selected Step outputs and rerun the selected range.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip selected Steps whose required outputs are already complete.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve dependencies and print commands without executing.",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Skip final execution confirmation.",
    )
    parser.add_argument(
        "--non-interactive",
        action="store_true",
        help="Never prompt. Missing required choices/inputs cause an error.",
    )

    return parser.parse_args(argv)


# =============================================================================
# Main
# =============================================================================


def resolve_scope(
    args: argparse.Namespace,
    *,
    interactive: bool,
) -> Tuple[int, int]:
    if args.step is not None:
        step = int(args.step)
        step_range(step, step)
        if args.end_step is not None:
            raise ValueError("--end-step cannot be used with --step")
        return step, step

    if args.start_step is not None:
        start = int(args.start_step)
        end = (
            int(args.end_step)
            if args.end_step is not None
            else start
        )
        step_range(start, end)
        return start, end

    if args.end_step is not None:
        raise ValueError("--end-step requires --start-step")

    if not interactive:
        raise ValueError(
            "Non-interactive mode requires --step or --start-step/--end-step."
        )

    return prompt_execution_scope()


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    print(f"ClonoDynamics orchestrator | {ORCHESTRATOR_VERSION}")
    args = parse_args(argv)
    interactive = not bool(args.non_interactive)

    start_step, end_step = resolve_scope(
        args,
        interactive=interactive,
    )
    selected_steps = step_range(start_step, end_step)

    if args.posterior_samples <= 0:
        raise ValueError("--posterior-samples must be > 0")

    # -------------------------------------------------------------------------
    # Code directory: require only scripts needed by selected Steps.
    # -------------------------------------------------------------------------
    code_dir = (
        args.code_dir.expanduser().resolve(strict=True)
        if args.code_dir is not None
        else Path(__file__).resolve().parent
    )
    if not code_dir.is_dir():
        raise NotADirectoryError(code_dir)

    selected_scripts = ensure_selected_scripts(
        code_dir,
        selected_steps,
    )

    # -------------------------------------------------------------------------
    # Raw data are needed only if Step 1 or Step 2 is selected.
    # -------------------------------------------------------------------------
    raw_needed = 1 in selected_steps or 2 in selected_steps

    data_path: Optional[Path] = None
    audit: Optional[RawDatasetAudit] = None

    if raw_needed:
        if args.data_path is not None:
            data_path = args.data_path.expanduser().resolve(strict=True)
            if not data_path.is_dir():
                raise NotADirectoryError(data_path)
        elif interactive:
            data_path = prompt_data_path()
        else:
            raise ValueError(
                "Selected range includes Step 1 or Step 2; --data-path is required."
            )

        audit = audit_raw_dataset(
            data_path,
            COMMON_PATTERN,
        )

        if 2 in selected_steps and audit.n_complete_pairs == 0:
            raise ValueError(
                "Step 2 is selected but no complete replicate pairs were found."
            )

    # -------------------------------------------------------------------------
    # Result root.
    # -------------------------------------------------------------------------
    # For raw-data runs, dataset type is needed before deriving the default root.
    if raw_needed:
        if args.dataset_type is not None:
            dataset_type = slug_dataset_type(args.dataset_type)
        elif interactive:
            dataset_type = prompt_dataset_type()
        else:
            raise ValueError(
                "--dataset-type is required in non-interactive raw-data runs."
            )

        results_root = (
            args.results_root.expanduser().resolve()
            if args.results_root is not None
            else default_results_root(
                data_path,  # type: ignore[arg-type]
                dataset_type,
            ).resolve()
        )

        existing_config = maybe_load_existing_config(results_root)

    else:
        if args.results_root is not None:
            results_root = args.results_root.expanduser().resolve()
        elif interactive:
            results_root = prompt_existing_directory(
                "Existing ClonoDynamics results root "
                "(folder containing the Step result directories)"
            )
        else:
            raise ValueError(
                "Starting at Step 3 or later requires --results-root."
            )

        if not results_root.exists() or not results_root.is_dir():
            raise NotADirectoryError(
                f"Existing results root not found: {results_root}"
            )

        existing_config = maybe_load_existing_config(results_root)

        inherited_type = infer_dataset_type_from_config(existing_config)
        if args.dataset_type is not None:
            dataset_type = slug_dataset_type(args.dataset_type)
        elif inherited_type is not None:
            dataset_type = inherited_type
            print(f"[config] Dataset type inherited: {dataset_type}")
        elif interactive:
            dataset_type = prompt_dataset_type()
        else:
            raise ValueError(
                "Dataset type is not available from existing config; "
                "provide --dataset-type."
            )

        data_path = infer_raw_path_from_config(existing_config)

    # -------------------------------------------------------------------------
    # Alpha is asked only if Step 2 is selected.
    # Otherwise retain inherited value for provenance when available.
    # -------------------------------------------------------------------------
    if 2 in selected_steps:
        if args.alpha is not None:
            alpha = float(args.alpha)
            if not (0.0 < alpha < 1.0):
                raise ValueError("--alpha must satisfy 0 < alpha < 1")
        elif interactive:
            alpha = prompt_alpha(0.05)
        else:
            alpha = 0.05
    else:
        if args.alpha is not None:
            alpha = float(args.alpha)
            if not (0.0 < alpha < 1.0):
                raise ValueError("--alpha must satisfy 0 < alpha < 1")
        else:
            alpha = infer_alpha_from_config(existing_config)

    # -------------------------------------------------------------------------
    # Discover only upstream artifacts that must already exist.
    # -------------------------------------------------------------------------
    resolved_upstream = resolve_upstream_inputs(
        selected_steps=selected_steps,
        results_root=results_root,
        step2_per_clone_override=args.step2_per_clone,
        step2_posterior_override=args.step2_posterior_dir,
        step4_trajectories_override=args.step4_trajectories,
        step5_transitions_override=args.step5_transitions,
        interactive=interactive,
    )

    inputs = effective_input_paths(
        selected_steps=selected_steps,
        results_root=results_root,
        resolved_upstream=resolved_upstream,
    )

    # -------------------------------------------------------------------------
    # Max dt: needed only by Step 5.
    # -------------------------------------------------------------------------
    max_dt = resolve_max_dt(
        selected_steps=selected_steps,
        explicit_max_dt=args.max_dt,
        audit=audit,
        inputs=inputs,
        existing_config=existing_config,
        interactive=interactive,
    )

    # -------------------------------------------------------------------------
    # Preflight summary.
    # -------------------------------------------------------------------------
    print_preflight(
        start_step=start_step,
        end_step=end_step,
        dataset_type=dataset_type,
        data_path=data_path if raw_needed else None,
        results_root=results_root,
        alpha=alpha,
        audit=audit,
        resolved_upstream=resolved_upstream,
        max_dt=max_dt,
    )

    # -------------------------------------------------------------------------
    # Fresh/resume applies ONLY to selected Step outputs.
    # -------------------------------------------------------------------------
    mode = choose_selected_run_mode(
        results_root,
        selected_steps,
        forced_restart=bool(args.restart),
        forced_resume=bool(args.resume),
        interactive=interactive,
    )

    if mode == "abort":
        print("Aborted.")
        return 0

    if mode == "fresh":
        clear_selected_outputs(
            results_root,
            selected_steps,
        )

    results_root.mkdir(parents=True, exist_ok=True)
    (results_root / "logs").mkdir(parents=True, exist_ok=True)
    (results_root / "orchestrator_runs").mkdir(
        parents=True,
        exist_ok=True,
    )

    # Recompute effective paths after cleanup. Cleanup touches only selected
    # outputs, not pre-existing upstream dependencies outside the range.
    inputs = effective_input_paths(
        selected_steps=selected_steps,
        results_root=results_root,
        resolved_upstream=resolved_upstream,
    )

    config = build_run_config(
        start_step=start_step,
        end_step=end_step,
        dataset_type=dataset_type,
        data_path=data_path,
        results_root=results_root,
        code_dir=code_dir,
        alpha=alpha,
        audit=audit,
        selected_scripts=selected_scripts,
        resolved_upstream=resolved_upstream,
        effective_inputs=inputs,
        max_dt=max_dt,
        n_jobs=int(args.n_jobs),
        n_jobs_profile=str(args.n_jobs_profile),
        posterior_samples=int(args.posterior_samples),
        seed=int(args.seed),
    )

    config_path = results_root / "00_orchestrator_config.json"
    run_config_path = (
        results_root
        / "orchestrator_runs"
        / (
            f"{timestamp_slug()}_steps_"
            f"{start_step}-{end_step}_config.json"
        )
    )

    payload = json.dumps(
        config,
        indent=2,
        ensure_ascii=False,
    )

    config_path.write_text(payload, encoding="utf-8")
    run_config_path.write_text(payload, encoding="utf-8")

    # -------------------------------------------------------------------------
    # Build commands for selected Steps only.
    # -------------------------------------------------------------------------
    commands: Dict[int, List[str]] = {}

    for step in selected_steps:
        commands[step] = build_step_command(
            step=step,
            code_dir=code_dir,
            results_root=results_root,
            inputs=inputs,
            data_path=data_path if raw_needed else None,
            alpha=alpha,
            file_sep=(
                audit.file_separator
                if audit is not None
                else None
            ),
            max_dt=max_dt,
            n_jobs=int(args.n_jobs),
            n_jobs_profile=str(args.n_jobs_profile),
            posterior_samples=int(args.posterior_samples),
            seed=int(args.seed),
            restart_step5=(
                mode == "fresh" and step == 5
            ),
        )

    outputs = required_outputs(results_root)

    print("\nExecution plan:")
    for step in selected_steps:
        state = (
            "SKIP (complete)"
            if mode == "resume"
            and output_is_complete(step, outputs)
            else "RUN"
        )
        print(f"\nStep {step}: {state}")
        print("  " + shell_join(commands[step]))

    if args.dry_run:
        print(
            "\n[DRY RUN] Dependency resolution and command construction "
            "completed. No Step was executed."
        )
        return 0

    if not args.yes:
        label = (
            f"Step {start_step}"
            if start_step == end_step
            else f"Steps {start_step}-{end_step}"
        )
        answer = input(
            f"\nRun {label} now? [Y/n]: "
        ).strip().lower()

        if answer not in {
            "",
            "y",
            "yes",
            "s",
            "si",
            "sì",
        }:
            print("Aborted.")
            return 0

    # -------------------------------------------------------------------------
    # Execution environment.
    # -------------------------------------------------------------------------
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        str(code_dir)
        if not existing_pythonpath
        else str(code_dir)
        + os.pathsep
        + existing_pythonpath
    )

    manifest_path = (
        results_root
        / "00_orchestrator_step_manifest.csv"
    )
    run_manifest_path = (
        results_root
        / "orchestrator_runs"
        / (
            f"{timestamp_slug()}_steps_"
            f"{start_step}-{end_step}_manifest.csv"
        )
    )

    manifest_rows: List[Dict[str, object]] = []

    # -------------------------------------------------------------------------
    # Execute selected contiguous range.
    # -------------------------------------------------------------------------
    for step in selected_steps:
        script = STEP_SCRIPTS[step]
        log_path = (
            results_root
            / "logs"
            / f"step{step}.log"
        )

        if (
            mode == "resume"
            and output_is_complete(step, outputs)
        ):
            row = {
                "step": step,
                "script": script,
                "status": "skipped_complete",
                "started_at": now_iso(),
                "finished_at": now_iso(),
                "returncode": 0,
                "log_path": str(log_path),
                "command": shell_join(commands[step]),
                "message": (
                    "Required selected-Step outputs already complete."
                ),
            }
            manifest_rows.append(row)
            write_manifest(manifest_path, manifest_rows)
            write_manifest(run_manifest_path, manifest_rows)
            print(
                f"\n[SKIP] Step {step}: "
                "required outputs already complete."
            )
            continue

        step_paths(results_root)[step].mkdir(
            parents=True,
            exist_ok=True,
        )

        # Inputs produced by earlier selected Steps must now exist.
        if step in {3, 4} and 2 in selected_steps:
            inputs["step2_per_clone"] = (
                canonical_upstream_paths(results_root)[
                    "step2_per_clone"
                ]
            )

        if step == 5:
            if 4 in selected_steps:
                inputs["step4_trajectories"] = (
                    canonical_upstream_paths(results_root)[
                        "step4_trajectories"
                    ]
                )
            if 2 in selected_steps:
                inputs["step2_posterior_dir"] = (
                    canonical_upstream_paths(results_root)[
                        "step2_posterior_dir"
                    ]
                )

            # Validate generated/current Step-5 inputs immediately before launch.
            validate_resolved_input(
                inputs["step4_trajectories"],
                label="Step-4 latent_trajectories_long.parquet",
                expect_dir=False,
            )
            validate_resolved_input(
                inputs["step2_posterior_dir"],
                label="Step-2 full posterior directory",
                expect_dir=True,
            )

            # Rebuild Step-5 command using current resolved paths.
            commands[step] = build_step_command(
                step=step,
                code_dir=code_dir,
                results_root=results_root,
                inputs=inputs,
                data_path=data_path if raw_needed else None,
                alpha=alpha,
                file_sep=(
                    audit.file_separator
                    if audit is not None
                    else None
                ),
                max_dt=max_dt,
                n_jobs=int(args.n_jobs),
                n_jobs_profile=str(args.n_jobs_profile),
                posterior_samples=int(args.posterior_samples),
                seed=int(args.seed),
                restart_step5=(
                    mode == "fresh"
                ),
            )

        if step == 6 and 5 in selected_steps:
            inputs["step5_transitions"] = (
                canonical_upstream_paths(results_root)[
                    "step5_transitions"
                ]
            )
            validate_resolved_input(
                inputs["step5_transitions"],
                label="Step-5 latent_transitions.parquet",
                expect_dir=False,
            )
            commands[step] = build_step_command(
                step=step,
                code_dir=code_dir,
                results_root=results_root,
                inputs=inputs,
                data_path=data_path if raw_needed else None,
                alpha=alpha,
                file_sep=(
                    audit.file_separator
                    if audit is not None
                    else None
                ),
                max_dt=max_dt,
                n_jobs=int(args.n_jobs),
                n_jobs_profile=str(args.n_jobs_profile),
                posterior_samples=int(args.posterior_samples),
                seed=int(args.seed),
                restart_step5=False,
            )

        # Validate Step-2 aggregate before Steps 3/4.
        if step in {3, 4}:
            validate_resolved_input(
                inputs["step2_per_clone"],
                label="Step-2 per_clone_latent_subject.parquet",
                expect_dir=False,
            )

        started = now_iso()

        print("\n" + "=" * 76)
        print(f"STEP {step} — {script}")
        print("=" * 76)

        returncode = run_and_tee(
            commands[step],
            log_path=log_path,
            cwd=code_dir,
            env=env,
        )

        finished = now_iso()
        status = (
            "completed"
            if returncode == 0
            else "failed"
        )
        message = ""

        if (
            returncode == 0
            and not output_is_complete(step, outputs)
        ):
            status = "failed_output_check"
            message = (
                "Process exited with code 0 but one or more required "
                "outputs are missing or empty."
            )

        row = {
            "step": step,
            "script": script,
            "status": status,
            "started_at": started,
            "finished_at": finished,
            "returncode": returncode,
            "log_path": str(log_path),
            "command": shell_join(commands[step]),
            "message": message,
        }

        manifest_rows.append(row)
        write_manifest(manifest_path, manifest_rows)
        write_manifest(run_manifest_path, manifest_rows)

        if status != "completed":
            print(
                f"\n[FAIL] Step {step} stopped the selected range."
            )
            print(f"See log:\n  {log_path}")
            if message:
                print(message)
            return (
                returncode
                if returncode != 0
                else 2
            )

        print(
            f"[OK] Step {step} completed and required outputs "
            "were verified."
        )

    # -------------------------------------------------------------------------
    # Done.
    # -------------------------------------------------------------------------
    label = (
        f"Step {start_step}"
        if start_step == end_step
        else f"Steps {start_step}-{end_step}"
    )

    print("\n" + "=" * 76)
    print(f"ClonoDynamics {label} completed successfully.")
    print(f"Results : {results_root}")
    print(f"Config  : {config_path}")
    print(f"Manifest: {manifest_path}")
    print("=" * 76)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
