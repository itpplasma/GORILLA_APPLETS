"""Regenerate the Homma-Hatayama derivation outputs and compare with the repository.

DERIVATIONS/homma_hatayama.wl writes, from the kernel proved in homma_hatayama_proofs.wl,
SRC/homma_series_coefficients.inc and TESTS/REFERENCE/homma/*.csv. This test
reruns it into a temporary directory (about 6 s) and requires the
outputs to be byte-identical to the committed files.

Environment: APPLETS_ROOT (repository), WOLFRAM_MATH (the math executable).
"""

import filecmp
import os
import subprocess
import sys
import tempfile
from pathlib import Path

GENERATED = [
    "SRC/homma_series_coefficients.inc",
    "TESTS/REFERENCE/homma/B01_oracle.csv",
    "TESTS/REFERENCE/homma/B02_oracle.csv",
    "TESTS/REFERENCE/homma/B03_oracle.csv",
    "TESTS/REFERENCE/homma/B07_targets.csv",
]


def main():
    root = Path(os.environ["APPLETS_ROOT"])
    math = os.environ.get("WOLFRAM_MATH", "math")
    script = root / "DERIVATIONS" / "homma_hatayama.wl"

    with tempfile.TemporaryDirectory() as out:
        run = subprocess.run([math, "-script", str(script), out],
                             capture_output=True, text=True)
        print(run.stdout)
        if "all checks passed" not in run.stdout:
            print(run.stderr)
            print(f"FAIL: derivation script did not pass its checks (return code {run.returncode})")
            return 1
        if run.returncode != 0:
            # The Wolfram kernel occasionally fails while shutting down, after
            # the script has written everything; the file comparison decides.
            print(f"warning: math returned {run.returncode} after all checks passed")

        stale = [rel for rel in GENERATED
                 if not filecmp.cmp(Path(out) / rel, root / rel, shallow=False)]
        if stale:
            print("FAIL: generated files differ from the repository:", *stale, sep="\n  ")
            print("Regenerate with: math -script DERIVATIONS/homma_hatayama.wl")
            return 1

    print("PASS: Homma derivation outputs are current")
    return 0


if __name__ == "__main__":
    sys.exit(main())
