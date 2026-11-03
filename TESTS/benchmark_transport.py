#!/usr/bin/env python3
"""Measure paired transport processes using the public mono-energetic CTest fixture."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time

import f90nml

from _ci_utils import StageTimer

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
INPUTS = ("gorilla.inp", "tetra_grid.inp", "gorilla_applets.inp",
          "mono_energetic_transp_coef.inp", "field_divB0.inp", "seed.inp",
          "fluxtubevolume.dat")
EQUILIBRIA = ("g_file_for_test", "convex_wall_for_test.dat", "netcdf_file_for_test.nc")
RTOL, ATOL = 1.e-8, 1.e-12


def command(argv):
    try:
        result = subprocess.run(argv, capture_output=True, text=True, timeout=10)
        return result.stdout.strip() if result.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        return None


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def positive_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def species(value):
    names = {"electron": 1, "deuterium": 2}
    if value in names:
        return names[value]
    if value in ("1", "2"):
        return int(value)
    raise argparse.ArgumentTypeError("choose electron/1 or deuterium/2")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--baseline-binary", type=Path,
                        help="optional older CPU executable; adds baseline-cpu without the new backend key")
    parser.add_argument("--fixture", required=True, type=Path,
                        help="prepared mono_energetic_transport_work from the native CTest")
    parser.add_argument("--output", required=True, type=Path, help="new directory; must not exist")
    parser.add_argument("--particles", nargs="+", type=positive, default=[128, 1024])
    parser.add_argument("--threads", nargs="+", type=positive, default=[1])
    parser.add_argument("--repeats", type=positive, default=3)
    parser.add_argument("--steps", type=positive, help="collisionless measurements (default: 32)")
    parser.add_argument("--time-step", type=positive_float,
                        help="collisionless measurement spacing in seconds (default: 1e-9)")
    parser.add_argument("--species", type=species, default=1)
    parser.add_argument("--backends", nargs="+", choices=["baseline-cpu", "cpu", "cpu-batch", "cuda"],
                        default=["cpu", "cpu-batch", "cuda"])
    parser.add_argument("--collisions", choices=["none", "lorentz"], default="none")
    parser.add_argument("--flight-time-multiplier", type=positive_float, default=1.,
                        help="Lorentz option 2 duration in transport timescales")
    parser.add_argument("--timeout", type=positive_float, default=600., help="seconds per process")
    parser.add_argument("--omp-places", default="cores")
    parser.add_argument("--omp-proc-bind", choices=["close", "spread", "true"], default="close")
    parser.add_argument("--compiler", help="compiler/version/flags supplied by the builder")
    parser.add_argument("--baseline-compiler", help="older executable compiler/version/flags")
    parser.add_argument("--source-revision", action="append", default=[], metavar="NAME=REV",
                        help="actual binary source revision, repeat for core/applets")
    parser.add_argument("--plot", action="store_true", help="write measured summary PNG/PDF")
    args = parser.parse_args()
    if args.collisions == "lorentz" and (args.steps is not None or args.time_step is not None):
        parser.error("option 2 resolves its own steps/time step; use --flight-time-multiplier")
    if args.steps is not None and args.steps < 4:
        parser.error("--steps must be at least 4 for the diffusion fit")
    for entry in args.source_revision:
        if "=" not in entry:
            parser.error("--source-revision requires NAME=REV")
    args.steps = args.steps or 32
    args.time_step = args.time_step or 1.e-9
    if not math.isfinite(args.steps * args.time_step):
        parser.error("total collisionless duration is not finite")
    args.particles = sorted(set(args.particles))
    args.threads = sorted(set(args.threads))
    args.backends = list(dict.fromkeys(args.backends))
    if args.baseline_binary is not None and "baseline-cpu" not in args.backends:
        args.backends.insert(0, "baseline-cpu")
    if "baseline-cpu" in args.backends and args.baseline_binary is None:
        parser.error("baseline-cpu requires --baseline-binary")
    return args


def prepare(case, args, particles, backend):
    case.mkdir()
    for name in INPUTS:
        shutil.copyfile(args.fixture / name, case / name)
    (case / "MHD_EQUILIBRIA").mkdir()
    for name in EQUILIBRIA:
        shutil.copyfile(args.fixture / "MHD_EQUILIBRIA" / name,
                        case / "MHD_EQUILIBRIA" / name)
    core = f90nml.read(case / "gorilla.inp")
    core["gorillanml"].update({
        "eps_phi": 0., "coord_system": 2, "ispecies": args.species,
        "ipusher": 2, "poly_order": 2, "boole_periodic_relocation": True,
        "boole_grid_for_find_tetra": False, "boole_adaptive_time_steps": False,
        "i_precomp": 0, "i_time_tracing_option": 1, "handover_processing_kind": 1,
        "boole_strong_electric_field": False})
    grid = f90nml.read(case / "tetra_grid.inp")
    grid["tetra_grid_nml"].update({
        "grid_kind": 3, "n1": 10, "n2": 10, "n3": 10,
        "boole_n_field_periods": True, "sfc_s_min": .1, "theta_geom_flux": 1,
        "theta0_at_xpoint": True,
        "g_file_filename": "MHD_EQUILIBRIA/g_file_for_test",
        "convex_wall_filename": "MHD_EQUILIBRIA/convex_wall_for_test.dat",
        "netcdf_filename": "MHD_EQUILIBRIA/netcdf_file_for_test.nc"})
    app = f90nml.read(case / "gorilla_applets.inp")
    app["gorilla_applets_nml"].update({
        "i_option": 2 if args.collisions == "lorentz" else 4,
        "filename_fluxtv_precomp": "fluxtubevolume.dat",
        "filename_fluxtv_load": "fluxtubevolume.dat"})
    core["gorillanml"]["boole_helical_pert"] = False
    mono = f90nml.read(case / "mono_energetic_transp_coef.inp")
    mono["transpcoefnml"].update({
        "transport_backend": backend, "i_integrator_type": 1, "n_particles": particles,
        "boole_collisions": args.collisions == "lorentz", "energy_ev": 3000.,
        "boole_random_precalc": True, "seed_option": 2, "random_seed_filename": "seed.inp",
        "idiffcoef_output": 3, "boole_psi_mat": True,
        "boole_write_particle_histories": args.collisions == "lorentz",
        "filename_particle_histories": "particle_histories.dat",
        "filename_transport_metadata": "transport_metadata.dat",
        "filename_transp_diff_coef": "nustar_diffcoef_std.dat",
        "filename_numerical_diff_coef": "numerical_diffcoef.dat",
        "filename_delta_s_squared": "mean_square_displacement.dat",
        "filename_std_dvt_delta_s_squared": "mean_square_uncertainty.dat",
        "nu_star": .1, "v_e": 0., "lambda_start": .9,
        "flight_time_multiplier": args.flight_time_multiplier,
        "nt_steps_numerical_diff": args.steps, "total_mc_time": args.steps * args.time_step})
    if backend == "baseline-cpu":
        del mono["transpcoefnml"]["transport_backend"]
    for filename, nml in (("gorilla.inp", core), ("tetra_grid.inp", grid),
                          ("gorilla_applets.inp", app), ("mono_energetic_transp_coef.inp", mono)):
        nml.write(case / filename, force=True)


def parse_profile(log):
    profile = {}
    for label in ("setup", "upload", "kernel", "download", "replay", "wall"):
        match = re.search(r"Transport batch " + label + r" time \[s\]\s*(" + FLOAT + ")", log)
        if match:
            profile[label + "_seconds"] = float(match[1].replace("D", "E").replace("d", "e"))
    patterns = {
        "attempted": r"Transport batch markers\s+(\d+)",
        "completed": r", completed\s+(\d+)",
        "physically_lost": r", physical losses\s+(\d+)",
        "orbit_advances": r"Transport batch orbit advances\s+(\d+)",
        "cell_crossings": r"Transport batch cell crossings\s+(\d+)",
        "replay_cells": r", replay cells\s+(\d+)",
        "numerically_failed": r"Transport batch numerical failures\s+(\d+)",
        "resolved_steps": r"Number of time steps \(measurements\)\s+(\d+)"}
    counters = {key: int(match[1]) for key, pattern in patterns.items()
                if (match := re.search(pattern, log))}
    return profile, counters


def numeric_rows(path):
    with path.open() as source:
        for line in source:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            row = [float(token.replace("D", "E").replace("d", "e")) for token in line.split()]
            if not all(math.isfinite(value) for value in row):
                raise ValueError(f"non-finite output in {path.name}")
            yield row


def output_files(args):
    common = ["mean_square_displacement.dat", "mean_square_uncertainty.dat", "n_lost_particles.dat"]
    if args.collisions == "lorentz":
        return common + ["nustar_diffcoef_std.dat", "particle_histories.dat", "transport_metadata.dat",
                         "loss_summary.dat", "lost_particle_events.dat"]
    return common + ["numerical_diffcoef.dat"]


def validate(case, args, record):
    errors = record["errors"]
    for filename in output_files(args):
        path = case / filename
        if not path.is_file():
            errors.append(f"missing output: {filename}")
            continue
        try:
            rows = sum(1 for _ in numeric_rows(path))
            if not rows and filename != "lost_particle_events.dat":
                errors.append(f"empty output: {filename}")
        except (ValueError, OSError) as exc:
            errors.append(str(exc))
    counts = record["counters"]
    try:
        loss_rows = list(numeric_rows(case / "n_lost_particles.dat"))
        if len(loss_rows) != 1 or len(loss_rows[0]) != 1 or loss_rows[0][0] != int(loss_rows[0][0]):
            raise ValueError("malformed n_lost_particles.dat")
        lost = int(loss_rows[0][0])
        if not 0 <= lost <= record["particles"]:
            raise ValueError("loss count outside ensemble")
        if record["backend"] in ("cpu", "baseline-cpu"):
            counts.update(attempted=record["particles"], completed=record["particles"] - lost,
                          reported_lost=lost, physically_lost=None, numerically_failed=None)
            record["failure_accounting"] = "unavailable: scalar tracer conflates numerical failure and loss"
        else:
            for key in ("attempted", "completed", "physically_lost", "orbit_advances", "cell_crossings", "replay_cells"):
                if key not in counts:
                    raise ValueError(f"missing batch counter: {key}")
            if (counts["attempted"] != record["particles"] or
                    counts["completed"] + counts["physically_lost"] != counts["attempted"] or
                    counts["physically_lost"] != lost):
                raise ValueError("batch terminal counts do not match ensemble/loss output")
            if counts.get("numerically_failed", 0):
                raise ValueError("batch reported numerical failures")
            if counts.get("numerically_failed") is None:
                counts["numerically_failed"] = 0 if record.get("returncode") == 0 else None
            record["failure_accounting"] = "batch terminal statuses; numerical failure aborts process"
            if len(record["profile"]) != 6:
                errors.append("incomplete batch timing profile")
        if counts.get("resolved_steps", 0) < 4:
            errors.append("fewer than four resolved measurements")
        if args.collisions == "none" and counts.get("resolved_steps") != args.steps:
            errors.append("resolved measurement count differs from requested count")
        if record["backend"] not in ("cpu", "baseline-cpu") and lost == 0:
            if counts["orbit_advances"] != record["particles"] * counts.get("resolved_steps", 0):
                raise ValueError("batch did not advance every no-loss marker at every measurement")
        cadence = args.time_step
        if args.collisions == "lorentz":
            rows = list(numeric_rows(case / "transport_metadata.dat"))
            if len(rows) != 1 or len(rows[0]) != 11:
                raise ValueError("malformed transport metadata")
            metadata = rows[0]
            if (metadata[7] != counts.get("resolved_steps") or
                    metadata[8] != record["particles"] or metadata[9] != 2 or metadata[6] <= 0):
                raise ValueError("resolved transport settings differ from the requested ensemble/seed")
            record["resolved_time_step_seconds"] = metadata[6]
            record["resolved_duration_seconds"] = metadata[6] * metadata[7]
            cadence = metadata[6]
            last_step = {}
            for row in numeric_rows(case / "particle_histories.dat"):
                if len(row) != 5 or row[1] != int(row[1]) or row[2] != int(row[2]):
                    raise ValueError("malformed particle history")
                marker, step = int(row[1]), int(row[2])
                if not 1 <= marker <= record["particles"] or step != last_step.get(marker, -1) + 1:
                    raise ValueError("particle history has missing or repeated samples")
                if abs(row[3] - step * metadata[6]) > 1.e-12:
                    raise ValueError("particle history time differs from the resolved cadence")
                last_step[marker] = step
            expected = {marker: int(metadata[7]) for marker in range(1, record["particles"] + 1)}
            events = list(numeric_rows(case / "lost_particle_events.dat"))
            if len(events) != lost:
                raise ValueError("loss event count differs from loss summary")
            seen = set()
            for row in events:
                if len(row) != 8 or row[1] != int(row[1]) or row[3] != int(row[3]):
                    raise ValueError("malformed loss event")
                marker, step = int(row[1]), int(row[3])
                if marker not in expected or marker in seen or not 1 <= step <= metadata[7]:
                    raise ValueError("loss event marker/step outside ensemble")
                seen.add(marker)
                expected[marker] = step - 1
            if last_step != expected:
                raise ValueError("particle histories do not cover every completed/lost marker")
        for filename in ("mean_square_displacement.dat", "mean_square_uncertainty.dat"):
            rows = numeric_rows(case / filename)
            for step in range(1, counts.get("resolved_steps", 0) + 1):
                row = next(rows, None)
                if row is None or len(row) != 2 or abs(row[0] - step * cadence) > 1.e-12:
                    raise ValueError(f"{filename}: missing measurement or inconsistent physical time")
            if next(rows, None) != [0., 0.] or next(rows, None) is not None:
                raise ValueError(f"{filename}: malformed end of measurement curve")
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    if any(not math.isfinite(value) or value < 0 for value in record["profile"].values()):
        errors.append("invalid profile duration")


def compare(case, reference, args):
    errors, maxima = [], {}
    from itertools import zip_longest
    for filename in output_files(args):
        maximum = 0.
        try:
            for index, (actual, expected) in enumerate(zip_longest(
                    numeric_rows(case / filename), numeric_rows(reference / filename)), 1):
                if actual is None or expected is None or len(actual) != len(expected):
                    raise ValueError(f"{filename}: row shape differs at {index}")
                for column, (a, b) in enumerate(zip(actual, expected)):
                    maximum = max(maximum, abs(a - b))
                    exact = (filename == "n_lost_particles.dat" or
                             filename == "loss_summary.dat" and column in (1, 2) or
                             filename == "lost_particle_events.dat" and column in (1, 2, 3))
                    if abs(a - b) > (0. if exact else ATOL + RTOL * abs(b)):
                        raise ValueError(f"{filename}: mismatch at row {index}: {a} versus {b}")
            maxima[filename] = maximum
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    return errors, maxima


def run_case(args, particles, threads, backend, replicate):
    phase = "cold_process" if replicate == 0 else "repeat_process"
    case = args.output / f"n{particles}_{backend}_t{threads}_{replicate:02d}"
    record = {"case": case.name, "backend": backend, "particles": particles,
              "threads": threads, "replicate": replicate, "phase": phase, "errors": [],
              "failure_accounting": "unavailable until successful terminal/output validation"}
    binary = args.baseline_binary if backend == "baseline-cpu" else args.binary
    record["binary"] = str(binary)
    timer = StageTimer("transport_benchmark", case / "timings.json")
    with timer.time("prepare_inputs"):
        prepare(case, args, particles, backend)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", OMP_PLACES=args.omp_places,
               OMP_PROC_BIND=args.omp_proc_bind, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    record["environment"] = {name: env.get(name) for name in (
        "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PLACES", "OMP_PROC_BIND", "OMP_THREAD_LIMIT",
        "OMP_WAIT_POLICY", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "CUDA_VISIBLE_DEVICES")}
    start = time.perf_counter()
    log = ""
    with (case / "process.log").open("w") as capture:
        try:
            process = subprocess.run([str(binary)], cwd=case, env=env, stdout=capture,
                                     stderr=subprocess.STDOUT, timeout=args.timeout)
            record["returncode"] = process.returncode
            if process.returncode:
                record["errors"].append(f"process exited {process.returncode}")
        except subprocess.TimeoutExpired:
            record["errors"].append(f"process timed out after {args.timeout} seconds")
            record["returncode"] = None
        except OSError as exc:
            record["errors"].append(str(exc))
            record["returncode"] = None
    record["wall_seconds"] = time.perf_counter() - start
    timer.stages["transport_process"] = record["wall_seconds"]
    timer.write()
    log = (case / "process.log").read_text(errors="replace")
    record["profile"], record["counters"] = parse_profile(log)
    record["counters"].setdefault("numerically_failed", None)
    if "Transport numerical failure" in log or "ERROR STOP" in log:
        record["errors"].append("transport reported numerical failure or error stop")
    validate(case, args, record)
    for key in ("attempted", "completed", "physically_lost", "numerically_failed"):
        record["counters"].setdefault(key, None)
    record["valid"] = not record["errors"]
    record["speedup_eligible"] = record["valid"] and backend not in ("cpu", "baseline-cpu")
    (case / "record.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def summarize(records):
    groups = {}
    for record in records:
        if record["phase"] == "repeat_process":
            groups.setdefault((record["particles"], record["backend"], record["threads"]), []).append(record)
    summary = []
    for (particles, backend, threads), rows in groups.items():
        valid = all(row["valid"] for row in records if
                    (row["particles"], row["backend"], row["threads"]) == (particles, backend, threads))
        values = [row["wall_seconds"] for row in rows]
        summary.append({"particles": particles, "backend": backend, "threads": threads,
                        "valid": valid, "repeats": len(rows),
                        "speedup_eligible": valid and all(row["speedup_eligible"] for row in rows),
                        "median_wall_seconds": statistics.median(values) if valid else None,
                        "min_wall_seconds": min(values) if valid else None,
                        "max_wall_seconds": max(values) if valid else None,
                        "speedup_vs_best_measured_cpu": None,
                        "speedup_vs_best_measured_cpu_batch": None})
    for particles in sorted({row["particles"] for row in summary}):
        for name, backends in (("cpu", ("baseline-cpu", "cpu", "cpu-batch")), ("cpu_batch", ("cpu-batch",))):
            baselines = [row for row in summary if row["particles"] == particles and
                         row["backend"] in backends and row["speedup_eligible"]]
            if not baselines:
                continue
            best = min(baselines, key=lambda row: row["median_wall_seconds"])
            for row in summary:
                if row["particles"] == particles and row["speedup_eligible"]:
                    row[name + "_baseline_threads"] = best["threads"]
                    row[name + "_baseline_backend"] = best["backend"]
                    row["speedup_vs_best_measured_" + name] = best["median_wall_seconds"] / row["median_wall_seconds"]
        for row in summary:
            if row["particles"] != particles or not row["speedup_eligible"]:
                continue
            for backend in ("baseline-cpu", "cpu", "cpu-batch"):
                peer = next((item for item in summary if item["particles"] == particles and
                             item["threads"] == row["threads"] and item["backend"] == backend and
                             item["speedup_eligible"]), None)
                row["speedup_vs_same_threads_" + backend.replace("-", "_")] = (
                    peer["median_wall_seconds"] / row["median_wall_seconds"] if peer else None)
    return summary


def plot_summary(args, summary, metadata):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("--plot requires matplotlib; raw JSON/CSV have been saved") from exc
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {"baseline-cpu": "#000000", "cpu": "#0072B2", "cpu-batch": "#009E73", "cuda": "#D55E00"}
    linestyles = {"baseline-cpu": ":", "cpu": "-", "cpu-batch": "--", "cuda": "-."}
    markers = ("o", "s", "^", "D", "v", "P", "X", "*")
    keys = sorted({(row["backend"], row["threads"]) for row in summary})
    for backend, threads in keys:
        rows = sorted((row for row in summary if row["valid"] and
                       (row["backend"], row["threads"]) == (backend, threads)), key=lambda row: row["particles"])
        color, style = colors[backend], linestyles[backend]
        marker = markers[args.threads.index(threads) % len(markers)]
        label = f"{backend}, {threads} CPU thread(s)"
        if rows:
            axes[0].errorbar([r["particles"] for r in rows], [r["median_wall_seconds"] for r in rows],
                             yerr=[[r["median_wall_seconds"] - r["min_wall_seconds"] for r in rows],
                                   [r["max_wall_seconds"] - r["median_wall_seconds"] for r in rows]],
                             color=color, marker=marker, linestyle=style, label=label, capsize=3)
        paired = [r for r in rows if r["speedup_vs_best_measured_cpu"] is not None]
        if paired:
            axes[1].plot([r["particles"] for r in paired],
                         [r["speedup_vs_best_measured_cpu"] for r in paired],
                         color=color, marker=marker, linestyle=style, label=label)
    for ax in axes:
        ax.set_xlabel("Markers (count)")
        ax.grid(alpha=.25)
        if ax.lines:
            ax.legend(fontsize=8)
    axes[0].set_ylabel("Whole-process wall time (s)")
    axes[0].set_ylim(bottom=0.)
    axes[0].set_title("Measured median; range across repeats")
    axes[1].set_ylabel("Speedup (dimensionless)")
    axes[1].set_title("Relative to fastest validated CPU configuration")
    if not any(row["valid"] for row in summary):
        axes[0].text(.5, .5, "No valid measurements", ha="center", transform=axes[0].transAxes)
    if not any(row["speedup_vs_best_measured_cpu"] is not None for row in summary):
        axes[1].text(.5, .5, "No valid paired measurements", ha="center", transform=axes[1].transAxes)
    axes[1].axhline(1., color=".5", linewidth=.8, linestyle=":")
    inventory = metadata.get("gpu")
    gpu = ", ".join(sorted({line.split(",")[0] for line in inventory.splitlines()})) if inventory else "GPU unavailable or not queried"
    fig.suptitle(f"{metadata['host']} | GPU inventory: {gpu}", fontsize=10)
    fig.savefig(args.output / "transport_benchmark.png", dpi=180)
    fig.savefig(args.output / "transport_benchmark.pdf")
    plt.close(fig)


def main():
    args = arguments()
    args.binary, args.fixture, args.output = (path.resolve() for path in (args.binary, args.fixture, args.output))
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        sys.exit("--binary must be an executable file")
    if args.baseline_binary is not None:
        args.baseline_binary = args.baseline_binary.resolve()
        if not args.baseline_binary.is_file() or not os.access(args.baseline_binary, os.X_OK):
            sys.exit("--baseline-binary must be an executable file")
    for path in [args.fixture / name for name in INPUTS] + [args.fixture / "MHD_EQUILIBRIA" / name for name in EQUILIBRIA]:
        if not path.is_file():
            sys.exit(f"missing prepared fixture file: {path}")
    try:
        args.output.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        sys.exit("--output already exists; choose a new directory")
    metadata = {
        "host": platform.node(), "platform": platform.platform(), "cpu": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "allowed_cpus": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "cpu_detail": command(["lscpu"]) or command(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "gpu": command(["nvidia-smi", "--query-gpu=name,uuid,driver_version,memory.total", "--format=csv,noheader"]),
        "nvcc_on_path": command(["nvcc", "--version"]), "binary": str(args.binary),
        "binary_sha256": sha256(args.binary), "compiler_supplied": args.compiler,
        "baseline_binary": str(args.baseline_binary) if args.baseline_binary else None,
        "baseline_binary_sha256": sha256(args.baseline_binary) if args.baseline_binary else None,
        "baseline_compiler_supplied": args.baseline_compiler,
        "source_revisions_supplied": dict(entry.split("=", 1) for entry in args.source_revision),
        "harness_revision": command(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"]),
        "harness_sha256": sha256(Path(__file__).resolve()),
        "fixture": str(args.fixture), "fixture_sha256": {
            str(path.relative_to(args.fixture)): sha256(path)
            for path in [args.fixture / name for name in INPUTS] +
            [args.fixture / "MHD_EQUILIBRIA" / name for name in EQUILIBRIA]},
        "settings": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "comparison": {"rtol": RTOL, "atol": ATOL,
                       "loss_fields": "counts and marker/reason/step identifiers exact; times and positions use numeric tolerances"},
        "wall_scope": "fresh executable: field/mesh setup, transport, transfers, recovery, diagnostics, shutdown",
        "repeat_scope": "one cold-process sample then independent processes; later processes may benefit from OS caches; no resident device warmup",
        "phase_scope": "batch profile only; CUDA device kernel events differ from host whole-process wall",
        "baseline_scope": "fastest validated baseline-cpu/cpu/cpu-batch configuration actually measured on this host; strict cpu-batch baseline also reported; scalar admitted only with zero reported losses and paired agreement with explicitly accounted no-loss cpu-batch",
        "backend_scope": "cpu/cpu-batch/cuda use --binary; baseline-cpu uses the separately identified older --baseline-binary with transport_backend omitted",
        "cold_scope": "first sample per configuration; OS/device cache state is not controlled",
        "scientific_scope": "fixed public fixture; numerical equivalence, not a physical convergence study"}
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    records = []
    for particles in args.particles:
        for replicate in range(args.repeats + 1):
            for backend in args.backends:
                for threads in args.threads:
                    record = run_case(args, particles, threads, backend, replicate)
                    records.append(record)
                    print(f"{record['case']}: {record['wall_seconds']:.6g} s; valid={record['valid']}", flush=True)
                    (args.output / "records.json").write_text(json.dumps(records, indent=2) + "\n")
    for particles in args.particles:
        candidates = [row for row in records if row["particles"] == particles and row["valid"]]
        reference = next((row for row in candidates if row["backend"] == "cpu-batch"), None)
        if reference is None:
            reference = next(iter(candidates), None)
        for row in candidates:
            row["reference_case"] = reference["case"]
            errors, maxima = compare(args.output / row["case"], args.output / reference["case"], args)
            row["max_absolute_output_errors"] = maxima
            row["errors"].extend(errors)
            row["valid"] = not row["errors"]
            accounted_peer = reference["backend"] == "cpu-batch" and reference["valid"]
            scalar_no_loss = (row["counters"].get("reported_lost") == 0 and
                              reference["counters"].get("physically_lost") == 0)
            scalar = row["backend"] in ("cpu", "baseline-cpu")
            row["speedup_eligible"] = row["valid"] and accounted_peer and (not scalar or scalar_no_loss)
            if scalar and row["speedup_eligible"]:
                row["completion_evidence"] = "inferred: no reported losses and output agrees with explicitly accounted no-loss cpu-batch; scalar numerical counter unavailable"
            (args.output / row["case"] / "record.json").write_text(json.dumps(row, indent=2) + "\n")
    summary = summarize(records)
    (args.output / "records.json").write_text(json.dumps(records, indent=2) + "\n")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output / "records.csv").open("w", newline="") as out:
        columns = ["case", "backend", "particles", "threads", "replicate", "phase", "wall_seconds", "valid", "speedup_eligible", "errors", "counters", "profile"]
        writer = csv.DictWriter(out, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in records:
            writer.writerow({key: json.dumps(value) if isinstance(value, (dict, list)) else value for key, value in row.items()})
    if args.plot:
        plot_summary(args, summary, metadata)
    invalid = sum(not row["valid"] for row in records)
    print(f"Results: {args.output}; invalid processes or pairs: {invalid}")
    return 1 if invalid else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as error:
        sys.exit(str(error))
