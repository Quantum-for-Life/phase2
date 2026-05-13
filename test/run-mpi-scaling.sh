#!/usr/bin/env bash
# run-mpi-scaling.sh
#
# Verify that the coeff_matrix Trotter path produces identical
# values when run with 1, 2, 4, 8 MPI ranks.  Catches ownership-
# filter bugs in state_prep_coeff_expand.
#
# Expects a pre-built test/t-circ_trott_coeff binary.

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
BIN="$HERE/t-circ_trott_coeff"
if [[ ! -x "$BIN" ]]; then
	echo "missing $BIN; run 'make build-test'" >&2
	exit 1
fi

MPIRUN=${MPIRUN:-mpirun}
OUT_DIR=$(mktemp -d)
trap 'rm -rf "$OUT_DIR"' EXIT

PASS=0
FAIL=0
for R in 1 2 4 8; do
	OUT="$OUT_DIR/r$R.log"
	if "$MPIRUN" -n "$R" "$BIN" > "$OUT" 2>&1; then
		echo "R=$R OK"
		PASS=$((PASS + 1))
	else
		echo "R=$R FAIL"
		cat "$OUT"
		FAIL=$((FAIL + 1))
	fi
done

if (( FAIL > 0 )); then
	exit 1
fi
echo "scaling OK: $PASS configurations"
