#!/usr/bin/env python3
"""Compare per-stage CI test timings between a baseline and a candidate.

Usage:
    compare_timings.py <baseline_dir> <candidate_dir> [--warn PCT] [--fail PCT]
                       [--summary PATH]

Both directories must contain one or more `<test>_timings.json` files of the
shape written by `TESTS/_ci_utils.py:StageTimer.write`, optionally with a
`.run<N>.json` suffix indicating one of several repeated runs. The best
(minimum) wall-clock time per stage is used for each side. The script prints
a Markdown table, writes it to --summary if given, and exits non-zero when
any stage exceeds the --fail percentage slowdown.
"""
import argparse
import json
import re
import sys
from pathlib import Path


RUN_SUFFIX_RE = re.compile(r"\.run\d+$")


def load_best(directory: Path) -> dict[str, dict[str, float]]:
    """Return {test_name: {stage: best_seconds}} over all runs in `directory`."""
    best: dict[str, dict[str, float]] = {}
    for path in sorted(directory.glob("*.json")):
        data = json.loads(path.read_text())
        test = data["test"]
        stages = data["stages"]
        slot = best.setdefault(test, {})
        for stage, secs in stages.items():
            if stage not in slot or secs < slot[stage]:
                slot[stage] = float(secs)
    return best


def render(
    baseline,
    candidate,
    warn_pct: float,
    fail_pct: float,
    min_duration: float,
) -> tuple[str, bool]:
    """Return (markdown_table, any_failure).

    Stages whose baseline runtime is below `min_duration` seconds are reported
    but never marked as regression (timing noise on short stages exceeds the
    percentage thresholds on hosted CI runners).
    """
    lines = [
        "| Test | Stage | Baseline (s) | Candidate (s) | Δ (s) | Δ (%) | Status |",
        "|------|-------|-------------:|--------------:|------:|------:|:------:|",
    ]
    any_fail = False
    tests = sorted(set(baseline) | set(candidate))
    for test in tests:
        stages = sorted(set(baseline.get(test, {})) | set(candidate.get(test, {})))
        for stage in stages:
            b = baseline.get(test, {}).get(stage)
            c = candidate.get(test, {}).get(stage)
            if b is None or c is None:
                lines.append(
                    f"| {test} | {stage} | "
                    f"{'-' if b is None else f'{b:.3f}'} | "
                    f"{'-' if c is None else f'{c:.3f}'} | "
                    f"- | - | :grey_question: missing |"
                )
                continue
            diff = c - b
            pct = 100.0 * diff / b if b > 0 else 0.0
            if b < min_duration:
                status = ":information_source: too short to gate"
            elif pct >= fail_pct:
                status = ":x: regression"
                any_fail = True
            elif pct >= warn_pct:
                status = ":warning: slow"
            elif pct <= -warn_pct:
                status = ":rocket: faster"
            else:
                status = ":white_check_mark: ok"
            lines.append(
                f"| {test} | {stage} | {b:.3f} | {c:.3f} | "
                f"{diff:+.3f} | {pct:+.1f}% | {status} |"
            )
    return "\n".join(lines) + "\n", any_fail


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline_dir", type=Path)
    ap.add_argument("candidate_dir", type=Path)
    ap.add_argument("--warn", type=float, default=20.0, help="warn threshold (%%)")
    ap.add_argument("--fail", type=float, default=50.0, help="fail threshold (%%)")
    ap.add_argument(
        "--min-duration",
        type=float,
        default=0.5,
        help="stages with baseline runtime below this (seconds) are reported "
        "but never gated",
    )
    ap.add_argument("--summary", type=Path, help="write markdown report to this path")
    args = ap.parse_args()

    baseline = load_best(args.baseline_dir)
    candidate = load_best(args.candidate_dir)

    table, any_fail = render(
        baseline, candidate, args.warn, args.fail, args.min_duration
    )
    header = (
        f"## Performance comparison\n\n"
        f"Baseline (best of N): `{args.baseline_dir}`\n"
        f"Candidate (best of N): `{args.candidate_dir}`\n"
        f"Thresholds: warn ≥ {args.warn:.0f}%, fail ≥ {args.fail:.0f}%, "
        f"gate skipped for stages < {args.min_duration:.2f} s\n\n"
    )
    report = header + table
    print(report)
    if args.summary:
        args.summary.write_text(report)
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
