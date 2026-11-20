#!/usr/bin/env bash
# Regenerate the delta-f kernels into a temporary directory and fail on any
# difference from the committed SRC/generated output. Needs fpm and SymEngine.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
committed="${here}/../../SRC/generated"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
(cd "${here}" && fpm run gen_deltaf_kernels -- "${tmp}")
status=0
for f in "${tmp}"/*.f90; do
    name="$(basename "${f}")"
    if ! diff -u "${committed}/${name}" "${f}"; then
        status=1
    fi
done
for f in "${committed}"/*.f90; do
    [[ -e "${tmp}/$(basename "${f}")" ]] || { echo "stale: ${f}"; status=1; }
done
exit "${status}"
