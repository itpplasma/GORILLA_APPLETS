#!/bin/bash
# Sourced by the native Slurm driver; both sides share the same allocation.
set -euo pipefail
export RUNS=3
export GITHUB_STEP_SUMMARY=$root/perf-summary.md
cd "$root"
git clone -q --no-hardlinks --dissociate "$src" candidate
git clone -q --no-hardlinks --dissociate "$src" baseline
git -C baseline checkout "${SCI_BASE_SHA:-origin/main}"
if [ -s candidate/GORILLA_REF ]; then
    read -r core_ref < candidate/GORILLA_REF
    [[ "$core_ref" =~ ^[0-9a-f]{40}$ ]]
else
    core_ref=${SCI_HEAD_REF:-main}
    git ls-remote --exit-code --heads https://github.com/itpplasma/GORILLA.git "$core_ref" >/dev/null || core_ref=main
fi
git clone https://github.com/itpplasma/GORILLA.git GORILLA
git -C GORILLA checkout "$core_ref"
python3 -m venv "$root/venv"
. "$root/venv/bin/activate"
python -m pip install --upgrade pip
pip install f90nml
mkdir -p "$src/ci-out"
trap 'if [ -d "$root/perf-results" ]; then cp -a "$root/perf-results" "$src/ci-out/"; fi' EXIT
echo 'Sync TESTS/ from candidate into baseline'
(cd "$root/."
rm -rf baseline/TESTS
cp -r candidate/TESTS baseline/TESTS
)
echo 'Build candidate'
(cd "$root/candidate"
make
)
echo 'Build baseline'
(cd "$root/baseline"
make configure
cmake --build BUILD --target gorilla_applets_main.x
)
echo 'Benchmark candidate (best of $RUNS)'
(cd "$root/candidate"
mkdir -p ../perf-results/candidate
for i in $(seq 1 "$RUNS"); do
  echo "::group::candidate run $i"
  # Run tests; ignore failures so we still capture timings for
  # regressions that affect physics but not run-to-completion.
  ctest --test-dir BUILD --output-on-failure || true
  for f in BUILD/TESTS/*/timings.json; do
    [ -f "$f" ] || continue
    name=$(basename "$(dirname "$f")" _work)
    cp "$f" "../perf-results/candidate/${name}.run${i}.json"
  done
  echo "::endgroup::"
done
)
echo 'Benchmark baseline (best of $RUNS)'
(cd "$root/baseline"
mkdir -p ../perf-results/baseline
for i in $(seq 1 "$RUNS"); do
  echo "::group::baseline run $i"
  ctest --test-dir BUILD --output-on-failure || true
  for f in BUILD/TESTS/*/timings.json; do
    [ -f "$f" ] || continue
    name=$(basename "$(dirname "$f")" _work)
    cp "$f" "../perf-results/baseline/${name}.run${i}.json"
  done
  echo "::endgroup::"
done
)
echo 'Compare timings'
(cd "$root/."
python3 candidate/TESTS/compare_timings.py \
  perf-results/baseline perf-results/candidate \
  --warn 20 --fail 50 \
  --summary perf-results/report.md
cat perf-results/report.md >> "$GITHUB_STEP_SUMMARY"
)
