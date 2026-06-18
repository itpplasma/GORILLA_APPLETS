"""Shared utilities for GORILLA_APPLETS CI tests.

Provides:
  - numeric file comparison against a committed reference, with tolerances;
  - timing helper that writes a JSON record consumable by the workflow.
"""
import json
import sys
import time
from pathlib import Path


def _tokens(line: str):
    return line.split()


def _parse_floats(tokens):
    try:
        return [float(t) for t in tokens]
    except ValueError:
        return None


def compare_numeric_file(
    observed: Path,
    reference: Path,
    rtol: float = 1.0e-10,
    atol: float = 0.0,
    sort_lines: bool = False,
) -> None:
    """Compare two whitespace-separated numeric files.

    Each line is tokenised; tokens are compared as floats with the given
    relative/absolute tolerance. Non-numeric tokens must match exactly.
    If `sort_lines` is True, both files are sorted line-wise before comparison
    (used for outputs whose row order depends on OpenMP scheduling).
    Exits the process with a non-zero status on mismatch.
    """
    if not reference.exists():
        sys.exit(f"FAIL: reference file missing: {reference}")
    if not observed.exists():
        sys.exit(f"FAIL: observed file missing: {observed}")

    obs_lines = [ln for ln in observed.read_text().splitlines() if ln.strip()]
    ref_lines = [ln for ln in reference.read_text().splitlines() if ln.strip()]

    if sort_lines:
        obs_lines = sorted(obs_lines)
        ref_lines = sorted(ref_lines)

    if len(obs_lines) != len(ref_lines):
        sys.exit(
            f"FAIL: line count mismatch in {observed.name}: "
            f"observed {len(obs_lines)} vs reference {len(ref_lines)}"
        )

    for i, (obs, ref) in enumerate(zip(obs_lines, ref_lines), start=1):
        obs_tok = _tokens(obs)
        ref_tok = _tokens(ref)
        if len(obs_tok) != len(ref_tok):
            sys.exit(
                f"FAIL: token count mismatch in {observed.name} line {i}: "
                f"observed {len(obs_tok)} vs reference {len(ref_tok)}"
            )
        obs_f = _parse_floats(obs_tok)
        ref_f = _parse_floats(ref_tok)
        if obs_f is None or ref_f is None:
            if obs_tok != ref_tok:
                sys.exit(
                    f"FAIL: non-numeric mismatch in {observed.name} line {i}: "
                    f"observed {obs_tok} vs reference {ref_tok}"
                )
            continue
        for j, (a, b) in enumerate(zip(obs_f, ref_f), start=1):
            diff = abs(a - b)
            tol = atol + rtol * abs(b)
            if diff > tol:
                sys.exit(
                    f"FAIL: value mismatch in {observed.name} "
                    f"line {i} token {j}: observed {a!r} vs reference {b!r} "
                    f"(|diff|={diff:.3e}, tol={tol:.3e})"
                )

    print(f"PASS: {observed.name} matches reference within tolerance.")


class StageTimer:
    """Accumulates wall-clock timings for named stages.

    If `output_path` is given, the JSON file is rewritten after each stage
    completes. This guarantees timings are persisted even if a later step
    (e.g. a reference-comparison check) raises SystemExit.
    """

    def __init__(self, test_name: str, output_path: "Path | None" = None):
        self.test_name = test_name
        self.stages = {}
        self.output_path = output_path

    def time(self, label: str):
        return _StageContext(self, label)

    def write(self, path: "Path | None" = None) -> None:
        target = path or self.output_path
        if target is None:
            return
        target.write_text(
            json.dumps(
                {"test": self.test_name, "stages": self.stages},
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )


class _StageContext:
    def __init__(self, timer: StageTimer, label: str):
        self.timer = timer
        self.label = label

    def __enter__(self):
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        elapsed = time.perf_counter() - self.t0
        self.timer.stages[self.label] = elapsed
        print(f"TIMING: {self.timer.test_name}/{self.label} = {elapsed:.3f} s")
        self.timer.write()
