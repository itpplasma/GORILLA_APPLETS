#!/usr/bin/env bash
# Mono-energetic radial transport coefficient: full two-pass run.
#
# Pass 1 (i_option=1) precomputes the flux-tube volume that seeds the
# Monte Carlo start positions. Pass 2 (i_option=3) scans the normalized
# collisionality nu* and writes the transport coefficient. Set MODE=single
# to compute a single nu* (i_option=2) instead of the scan.
#
# All work happens in a fresh run_output/ subdirectory so the tracked
# input files in this example stay untouched.
set -euo pipefail
cd "$(dirname "$0")"

: "${OMP_NUM_THREADS:=$(nproc)}"
export OMP_NUM_THREADS
MODE="${MODE:-scan}"   # scan -> i_option=3, single -> i_option=2

exe=./gorilla_applets_main.x
if [ ! -x "$exe" ]; then
    exe=../../BUILD/gorilla_applets_main.x
fi
if [ ! -x "$exe" ]; then
    echo "error: gorilla_applets_main.x not found. Build first: (cd ../.. && make)" >&2
    exit 1
fi
exe="$(cd "$(dirname "$exe")" && pwd)/$(basename "$exe")"

run=run_output
rm -rf "$run"
mkdir "$run"
for f in gorilla.inp tetra_grid.inp mono_energetic_transp_coef.inp seed.inp \
         field_divB0.inp preload_for_SYNCH.inp MHD_EQUILIBRIA; do
    ln -s "../$f" "$run/$f"
done
cp gorilla_applets.inp "$run/gorilla_applets.inp"
cd "$run"

# Portable in-place edit: GNU `sed -i -E` and BSD `sed -i '' -E` are
# incompatible, so write to a temp file and move it back instead.
set_option() {
    local tmp
    tmp="$(mktemp)"
    sed -E "s/^( *i_option *=).*/\1 $1 ,/" gorilla_applets.inp > "$tmp"
    mv "$tmp" gorilla_applets.inp
}

echo "== pass 1: flux-tube volume (i_option=1) =="
set_option 1
"$exe"

if [ "$MODE" = single ]; then
    echo "== pass 2: single nu* transport coefficient (i_option=2) =="
    set_option 2
else
    echo "== pass 2: nu* scan (i_option=3) =="
    set_option 3
fi
"$exe"

echo
echo "done. results in $(pwd):"
echo "  nustar_diffcoef_std.dat   transport coefficient D11 vs nu*"
echo "  psi2.dat, std_psi2.dat    (Delta s)^2 trends used for the fit"
