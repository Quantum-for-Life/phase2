#!/usr/bin/env bash
# check-docs.sh
#
# Verify that every new public symbol declared in
# include/phase2/*.h is documented in either
# doc/simul-h5-specs.md or doc/state-prep.md.
#
# This is a coverage gate, not a content check: a name appearing
# in a doc text is enough.  Failures are intended to catch
# silently added APIs.

set -euo pipefail

HERE=$(cd "$(dirname "$0")/.." && pwd)
INC="$HERE/include/phase2"
DOCS=("$HERE/doc/simul-h5-specs.md" "$HERE/doc/state-prep.md")

# Symbols introduced by this change set; widen the list when
# adding new public APIs.
SYMBOLS=(
	data_state_prep_kind
	data_coeff_matrix_getnums
	data_coeff_matrix_read
	data_coeff_matrix_csf_count
	data_coeff_matrix_csf_read
	state_prep_coeff_expand
	state_prep_coeff_expand_all
	state_prep_coeff_inner
	det_small
	combinations_init
	combinations_next
	circ_coeff_init
	circ_coeff_free
	STPREP_MULTIDET
	STPREP_COEFF_MATRIX
)

FAIL=0
for sym in "${SYMBOLS[@]}"; do
	if ! grep -q -F "$sym" "${DOCS[@]}"; then
		echo "MISSING: $sym not mentioned in any doc" >&2
		FAIL=$((FAIL + 1))
	fi
done

if (( FAIL > 0 )); then
	exit 1
fi
echo "all new symbols are documented"
