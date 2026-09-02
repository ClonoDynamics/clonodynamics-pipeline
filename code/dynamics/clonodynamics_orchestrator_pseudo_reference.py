#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClonoDynamics pseudo-reference orchestrator (Steps 7-8).

Managed workflows
-----------------
- Steps 7-8
- Step 7 only
- Step 8 only

Dependencies
------------
Step 7:
    pseudo Step-4 latent_trajectories_long.parquet
    pseudo Step-5 latent_transitions.parquet

Step 8:
    pseudo Step-5 latent_transitions.parquet

Step 8 does not consume Step-7 outputs. Step 7 establishes the conditioning
strategy scientifically; Step 8 independently evaluates displacement
representations with xstar conditioning fixed.

Step 6 is descriptive and is not required.

Default Step-7 production run
-----------------------------
--stage all
--expected-pseudo-subjects 3
--crossfit-dt 1
--n-order-pairs 2000
--seed 123

Default Step-8 production run
-----------------------------
--conditioning-col xstar_latent
--representations latent_frequency=dx_latent,observed_frequency=dx_freq_obs,log_count=dx_count_sum
--relative-frequency-representations latent_frequency,observed_frequency
--closure-col closure_component
--dt-values 1
--classes TT
--min-n 200
--min-common-bins 10
--min-n-subject 50
--min-common-bins-subject 6
--n-subject-bootstrap 2000
--seed 123

No Step-7 edge file is injected into Step 8.

Interactive use
---------------
python3 clonodynamics_orchestrator_pseudo_reference.py
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

VERSION = "v1-pseudo-reference-steps7-8-2026-08-25"

STEP7 = "7-pseudo_longitudinal_conditioning_benchmark.py"
STEP8 = "8-pseudo_longitudinal_representation_assessment.py"

DIR4 = "4-latent_trajectory_construction"
DIR5 = "5-latent_transition_construction"
DIR7 = "7-pseudo_longitudinal_conditioning_benchmark"
DIR8 = "8-pseudo_longitudinal_representation_assessment"

TRAJ = "latent_trajectories_long.parquet"
TRANS = "latent_transitions.parquet"

REQ7 = [
    "00_run_config.json",
    "00_benchmark_summary.csv",
    "00_pipeline_manifest.csv",
    "README_outputs.md",
    "02_replicate_decoupling/forward_drift_bin_edges.csv",
    "03_ordering_robustness/pseudo_reverse_pair_envelope.csv",
]
REQ8 = [
    "00_run_config.json",
    "02_conditioning_bin_edges.csv",
    "03_pooled_representation_curves.csv",
    "04_pooled_representation_assessment.csv",
    "11_representation_assessment.csv",
    "representation_assessment_report.md",
]


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def stamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def shell_join(cmd: Sequence[str]) -> str:
    try:
        return shlex.join([str(x) for x in cmd])
    except AttributeError:
        return " ".join(shlex.quote(str(x)) for x in cmd)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def workflow_from_text(value: str) -> Tuple[int, ...]:
    v = value.strip().lower().replace(" ", "")
    mapping = {
        "1": (7, 8), "7-8": (7, 8), "all": (7, 8),
        "2": (7,), "7": (7,), "step7": (7,),
        "3": (8,), "8": (8,), "step8": (8,),
    }
    if v not in mapping:
        raise ValueError("workflow must be 7-8, 7, or 8")
    return mapping[v]


def workflow_label(steps: Sequence[int]) -> str:
    t = tuple(steps)
    if t == (7, 8):
        return "Steps 7-8"
    return f"Step {t[0]}"


def prompt_workflow() -> Tuple[int, ...]:
    print("\nPseudo-reference workflow:")
    print("  [1] Steps 7-8")
    print("  [2] Step 7 only")
    print("  [3] Step 8 only")
    while True:
        raw = input("Selection [1]: ").strip()
        if not raw:
            return (7, 8)
        try:
            return workflow_from_text(raw)
        except ValueError:
            print("Please enter 1, 2, or 3.")


def looks_like_root(path: Path) -> bool:
    return any(
        (path / name).exists()
        for name in (
            "00_orchestrator_config.json",
            DIR4, DIR5, DIR7, DIR8,
        )
    )


def normalize_root(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)
    if looks_like_root(path):
        return path
    child = path / "pseudo-longitudinal"
    if child.is_dir() and looks_like_root(child):
        print(f"[root] Using pseudo-longitudinal child:\n  {child}")
        return child.resolve()
    return path


def prompt_root() -> Path:
    while True:
        raw = input("\nExisting pseudo-longitudinal results root: ").strip()
        if not raw:
            print("A path is required.")
            continue
        try:
            return normalize_root(Path(raw))
        except Exception as exc:
            print(f"Invalid path: {exc}")


def validate_pseudo_provenance(root: Path) -> Dict[str, object]:
    cfg = root / "00_orchestrator_config.json"
    if not cfg.exists():
        print("[warning] No Steps-1-6 orchestrator config: pseudo provenance cannot be verified.")
        return {"config_found": False, "path": str(cfg), "dataset_type": None}
    with cfg.open("r", encoding="utf-8") as f:
        data = json.load(f)
    dtype = data.get("dataset_type")
    if dtype is not None:
        norm = str(dtype).lower().replace("_", "-")
        if norm not in {"pseudo-longitudinal", "pseudolongitudinal"}:
            raise ValueError(
                f"Result tree is not pseudo-longitudinal: dataset_type={dtype!r}"
            )
    return {"config_found": True, "path": str(cfg), "dataset_type": dtype}


def validate_scripts(code_dir: Path, steps: Sequence[int]) -> Dict[int, Path]:
    out = {}
    if 7 in steps:
        out[7] = code_dir / STEP7
    if 8 in steps:
        out[8] = code_dir / STEP8
    missing = [p for p in out.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing selected scripts:\n" + "\n".join(f"  - {p}" for p in missing)
        )
    return out


def validate_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"{label} missing/empty:\n  {path}")
    return path


def find_candidates(root: Path, basename: str) -> List[Path]:
    return sorted(
        {p.resolve() for p in root.rglob(basename) if p.is_file() and p.stat().st_size > 0},
        key=lambda p: (len(p.parts), str(p)),
    )


def choose_file(label: str, candidates: Sequence[Path]) -> Path:
    print(f"\nRequired upstream input: {label}")
    if candidates:
        for i, p in enumerate(candidates, 1):
            print(f"  [{i}] {p}")
        print("  [0] enter another path")
        while True:
            raw = input("Selection: ").strip()
            try:
                n = int(raw)
            except ValueError:
                print("Enter a number.")
                continue
            if n == 0:
                break
            if 1 <= n <= len(candidates):
                return candidates[n - 1]
            print("Invalid selection.")
    while True:
        raw = input(f"Path for {label}: ").strip()
        try:
            return validate_file(Path(raw), label)
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
    candidates = find_candidates(root, basename)
    if len(candidates) == 1:
        print(f"[input] {label}:\n  {candidates[0]}")
        return candidates[0]
    if interactive:
        return choose_file(label, candidates)
    if not candidates:
        raise FileNotFoundError(f"Could not find {label}; checked {canonical}")
    raise RuntimeError(
        f"Multiple candidates for {label}; use an explicit override:\n"
        + "\n".join(f"  - {p}" for p in candidates)
    )


def resolve_inputs(
    root: Path,
    steps: Sequence[int],
    trajectories_override: Optional[Path],
    transitions_override: Optional[Path],
    interactive: bool,
) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    out["transitions"] = resolve_file(
        root,
        root / DIR5 / TRANS,
        TRANS,
        "Pseudo Step-5 latent_transitions.parquet",
        transitions_override,
        interactive,
    )
    if 7 in steps:
        out["trajectories"] = resolve_file(
            root,
            root / DIR4 / TRAJ,
            TRAJ,
            "Pseudo Step-4 latent_trajectories_long.parquet",
            trajectories_override,
            interactive,
        )
    return out


def prompt_expected_subjects() -> int:
    while True:
        raw = input("\nExpected pseudo biological subjects [3; 0 disables check]: ").strip()
        if not raw:
            return 3
        try:
            n = int(raw)
        except ValueError:
            print("Enter a non-negative integer.")
            continue
        if n >= 0:
            return n
        print("Enter a non-negative integer.")


def outdir(root: Path, step: int) -> Path:
    return root / (DIR7 if step == 7 else DIR8)


def required_outputs(root: Path, step: int) -> List[Path]:
    base = outdir(root, step)
    rels = REQ7 if step == 7 else REQ8
    return [base / x for x in rels]


def complete(root: Path, step: int) -> bool:
    for p in required_outputs(root, step):
        if not p.exists():
            return False
        if p.is_file() and p.stat().st_size <= 0:
            return False
    return True


def selected_outputs_exist(root: Path, steps: Sequence[int]) -> bool:
    return any(outdir(root, s).is_dir() and any(outdir(root, s).iterdir()) for s in steps)


def clear_selected(root: Path, steps: Sequence[int]) -> None:
    for s in steps:
        d = outdir(root, s)
        if d.exists():
            shutil.rmtree(d)
        log = root / "logs" / f"step{s}.log"
        if log.exists():
            log.unlink()


def choose_mode(root: Path, steps: Sequence[int], restart: bool, resume: bool, interactive: bool) -> str:
    if restart and resume:
        raise ValueError("Use only one of --restart or --resume.")
    if not selected_outputs_exist(root, steps):
        return "fresh"
    if restart:
        return "fresh"
    if resume:
        return "resume"
    if not interactive:
        raise RuntimeError("Selected outputs exist; use --restart or --resume.")
    print(f"\nExisting outputs detected for {workflow_label(steps)}:")
    print("  [1] fresh   — delete only selected Step outputs")
    print("  [2] resume  — skip complete selected Steps")
    print("  [3] abort")
    while True:
        raw = input("Selection [2]: ").strip().lower()
        if raw in {"", "2", "resume", "r"}:
            return "resume"
        if raw in {"1", "fresh", "f"}:
            return "fresh"
        if raw in {"3", "abort", "a", "q"}:
            return "abort"
        print("Please enter 1, 2, or 3.")


def cmd7(
    script: Path, inputs: Dict[str, Path], root: Path,
    expected_subjects: int, n_order_pairs: int, seed: int, restart: bool
) -> List[str]:
    cmd = [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--transitions", str(inputs["transitions"]),
        "--trajectories", str(inputs["trajectories"]),
        "--outdir", str(root / DIR7),
        "--stage", "all",
        "--expected-pseudo-subjects", str(expected_subjects),
        "--crossfit-dt", "1",
        "--n-order-pairs", str(n_order_pairs),
        "--seed", str(seed),
    ]
    if restart:
        cmd.append("--restart")
    return cmd


def cmd8(
    script: Path, inputs: Dict[str, Path], root: Path,
    n_boot: int, seed: int
) -> List[str]:
    return [
        str(Path(sys.executable).resolve()), "-u", str(script),
        "--input", str(inputs["transitions"]),
        "--outdir", str(root / DIR8),
        "--conditioning-col", "xstar_latent",
        "--representations",
        "latent_frequency=dx_latent,observed_frequency=dx_freq_obs,log_count=dx_count_sum",
        "--relative-frequency-representations", "latent_frequency,observed_frequency",
        "--closure-col", "closure_component",
        "--dt-values", "1",
        "--classes", "TT",
        "--min-n", "200",
        "--min-common-bins", "10",
        "--min-n-subject", "50",
        "--min-common-bins-subject", "6",
        "--n-subject-bootstrap", str(n_boot),
        "--seed", str(seed),
    ]


def run_tee(cmd: Sequence[str], log: Path, cwd: Path, env: Dict[str, str]) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    print("\n$ " + shell_join(cmd))
    print(f"[log] {log}")
    with log.open("w", encoding="utf-8", buffering=1) as fh:
        p = subprocess.Popen(
            list(cmd), cwd=str(cwd), env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        assert p.stdout is not None
        for line in p.stdout:
            print(line, end="")
            fh.write(line)
        return int(p.wait())


def write_manifest(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = ["step", "script", "status", "started_at", "finished_at",
              "returncode", "log_path", "command", "message"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def print_preflight(
    steps: Sequence[int], root: Path, inputs: Dict[str, Path],
    expected_subjects: Optional[int], n_order_pairs: int, n_boot: int, seed: int
) -> None:
    print("\n" + "=" * 78)
    print(f"ClonoDynamics pseudo-reference — {workflow_label(steps)}")
    print("=" * 78)
    print("Dataset type           : pseudo-longitudinal")
    print(f"Results root           : {root}")
    print("Step 6 required        : no")
    print("Step 8 depends on 7    : no")
    print("\nResolved upstream inputs:")
    if "trajectories" in inputs:
        print(f"  Step-4 trajectories  : {inputs['trajectories']}")
    print(f"  Step-5 transitions   : {inputs['transitions']}")
    if 7 in steps:
        print("\nStep 7:")
        print("  stage                : all")
        print(f"  expected subjects    : {expected_subjects}")
        print("  crossfit dt          : 1")
        print(f"  order/reverse pairs  : {n_order_pairs}")
        print("  bundled engines      : yes")
    if 8 in steps:
        print("\nStep 8:")
        print("  conditioning         : xstar_latent")
        print("  dt / class           : 1 / TT")
        print("  representations      : latent_frequency, observed_frequency, log_count")
        print(f"  subject bootstrap    : {n_boot}")
        print("  Step-7 edges injected: no")
    print(f"\nRandom seed            : {seed}")
    print("=" * 78)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="ClonoDynamics pseudo-reference orchestrator for Steps 7-8.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--workflow", choices=["7-8", "7", "8"], default=None)
    p.add_argument("--results-root", type=Path, default=None)
    p.add_argument("--code-dir", type=Path, default=None)
    p.add_argument("--trajectories", type=Path, default=None)
    p.add_argument("--transitions", type=Path, default=None)
    p.add_argument("--expected-pseudo-subjects", type=int, default=None)
    p.add_argument("--n-order-pairs", type=int, default=2000)
    p.add_argument("--step8-subject-bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--restart", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--yes", action="store_true")
    p.add_argument("--non-interactive", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    print(f"ClonoDynamics pseudo-reference orchestrator | {VERSION}")
    args = parse_args(argv)
    interactive = not args.non_interactive

    if args.workflow:
        steps = workflow_from_text(args.workflow)
    elif interactive:
        steps = prompt_workflow()
    else:
        raise ValueError("--workflow is required in non-interactive mode.")

    if args.results_root:
        root = normalize_root(args.results_root)
    elif interactive:
        root = prompt_root()
    else:
        raise ValueError("--results-root is required in non-interactive mode.")

    provenance = validate_pseudo_provenance(root)

    code_dir = (
        Path(__file__).resolve().parent
        if args.code_dir is None
        else args.code_dir.expanduser().resolve(strict=True)
    )
    scripts = validate_scripts(code_dir, steps)

    inputs = resolve_inputs(
        root, steps, args.trajectories, args.transitions, interactive
    )

    expected_subjects: Optional[int] = None
    if 7 in steps:
        if args.expected_pseudo_subjects is not None:
            expected_subjects = int(args.expected_pseudo_subjects)
            if expected_subjects < 0:
                raise ValueError("--expected-pseudo-subjects must be >= 0")
        elif interactive:
            expected_subjects = prompt_expected_subjects()
        else:
            expected_subjects = 3

    if args.n_order_pairs < 1:
        raise ValueError("--n-order-pairs must be >= 1")
    if args.step8_subject_bootstrap < 1:
        raise ValueError("--step8-subject-bootstrap must be >= 1")

    print_preflight(
        steps, root, inputs, expected_subjects,
        args.n_order_pairs, args.step8_subject_bootstrap, args.seed
    )

    mode = choose_mode(root, steps, args.restart, args.resume, interactive)
    if mode == "abort":
        print("Aborted.")
        return 0
    if mode == "fresh":
        clear_selected(root, steps)

    (root / "logs").mkdir(parents=True, exist_ok=True)
    run_dir = root / "orchestrator_runs_pseudo_reference"
    run_dir.mkdir(parents=True, exist_ok=True)

    commands: Dict[int, List[str]] = {}
    if 7 in steps:
        assert expected_subjects is not None
        commands[7] = cmd7(
            scripts[7], inputs, root, expected_subjects,
            args.n_order_pairs, args.seed, mode == "fresh"
        )
    if 8 in steps:
        commands[8] = cmd8(
            scripts[8], inputs, root, args.step8_subject_bootstrap, args.seed
        )

    cfg = {
        "orchestrator": Path(__file__).name,
        "orchestrator_version": VERSION,
        "generated_at": now(),
        "dataset_type": "pseudo-longitudinal",
        "workflow": workflow_label(steps),
        "selected_steps": list(steps),
        "results_root": str(root),
        "code_dir": str(code_dir),
        "upstream_inputs": {k: str(v) for k, v in inputs.items()},
        "steps_1_6_provenance": provenance,
        "script_sha256": {str(k): sha256(v) for k, v in scripts.items()},
        "step7": {
            "selected": 7 in steps,
            "stage": "all",
            "expected_pseudo_subjects": expected_subjects,
            "crossfit_dt": 1,
            "n_order_pairs": args.n_order_pairs,
            "bundled_engines": True,
        },
        "step8": {
            "selected": 8 in steps,
            "conditioning": "xstar_latent",
            "representations": {
                "latent_frequency": "dx_latent",
                "observed_frequency": "dx_freq_obs",
                "log_count": "dx_count_sum",
            },
            "dt_values": [1],
            "classes": ["TT"],
            "n_subject_bootstrap": args.step8_subject_bootstrap,
            "step7_edges_injected": False,
            "automatic_winner": False,
        },
        "seed": args.seed,
        "dependency_graph": {
            "step7": ["Step-4 trajectories", "Step-5 transitions"],
            "step8": ["Step-5 transitions"],
            "step8_consumes_step7_output": False,
            "step6_required": False,
        },
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
    }

    cfg_text = json.dumps(cfg, indent=2, ensure_ascii=False)
    cfg_path = root / "00_orchestrator_pseudo_reference_config.json"
    cfg_path.write_text(cfg_text, encoding="utf-8")
    (run_dir / f"{stamp()}_{'-'.join(map(str, steps))}_config.json").write_text(
        cfg_text, encoding="utf-8"
    )

    print("\nExecution plan:")
    for s in steps:
        state = "SKIP (complete)" if mode == "resume" and complete(root, s) else "RUN"
        print(f"\nStep {s}: {state}")
        print("  " + shell_join(commands[s]))

    if args.dry_run:
        print("\n[DRY RUN] Inputs resolved and commands built. Nothing executed.")
        return 0

    if not args.yes:
        answer = input(f"\nRun {workflow_label(steps)} now? [Y/n]: ").strip().lower()
        if answer not in {"", "y", "yes", "s", "si", "sì"}:
            print("Aborted.")
            return 0

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    old_pp = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(code_dir) if not old_pp else str(code_dir) + os.pathsep + old_pp

    manifest = root / "00_orchestrator_pseudo_reference_manifest.csv"
    run_manifest = run_dir / f"{stamp()}_{'-'.join(map(str, steps))}_manifest.csv"
    rows: List[Dict[str, object]] = []

    for s in steps:
        script_name = STEP7 if s == 7 else STEP8
        log = root / "logs" / f"step{s}.log"

        if mode == "resume" and complete(root, s):
            row = {
                "step": s, "script": script_name, "status": "skipped_complete",
                "started_at": now(), "finished_at": now(), "returncode": 0,
                "log_path": str(log), "command": shell_join(commands[s]),
                "message": "Required completion markers already present.",
            }
            rows.append(row)
            write_manifest(manifest, rows)
            write_manifest(run_manifest, rows)
            print(f"\n[SKIP] Step {s}: complete.")
            continue

        outdir(root, s).mkdir(parents=True, exist_ok=True)
        started = now()
        print("\n" + "=" * 78)
        print(f"STEP {s} — {script_name}")
        print("=" * 78)

        rc = run_tee(commands[s], log, code_dir, env)
        finished = now()
        status = "completed" if rc == 0 else "failed"
        message = ""

        if rc == 0 and not complete(root, s):
            status = "failed_output_check"
            missing = [
                str(p) for p in required_outputs(root, s)
                if not p.exists() or (p.is_file() and p.stat().st_size <= 0)
            ]
            message = "Missing completion outputs: " + "; ".join(missing)

        row = {
            "step": s, "script": script_name, "status": status,
            "started_at": started, "finished_at": finished, "returncode": rc,
            "log_path": str(log), "command": shell_join(commands[s]),
            "message": message,
        }
        rows.append(row)
        write_manifest(manifest, rows)
        write_manifest(run_manifest, rows)

        if status != "completed":
            print(f"\n[FAIL] Step {s} stopped the workflow.")
            print(f"See log:\n  {log}")
            if message:
                print(message)
            return rc if rc != 0 else 2

        print(f"[OK] Step {s} completed and required outputs were verified.")

    print("\n" + "=" * 78)
    print("ClonoDynamics pseudo-reference workflow completed successfully.")
    print(f"Workflow : {workflow_label(steps)}")
    print(f"Results  : {root}")
    print(f"Config   : {cfg_path}")
    print(f"Manifest : {manifest}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
