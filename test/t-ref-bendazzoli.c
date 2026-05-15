#include "c23_compat.h"
#include <complex.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"

#include "test.h"

/*
 * Reference test for state_prep_coeff_expand_all on a small
 * multi-component CSF fixture (n_qubits = 8, 2 closed-shell-per-
 * component blocks).
 *
 * Fixtures (committed under test/data/bendazzoli/n4_oss_k0/):
 *
 *   n4_oss_k0.h5            -- /state_prep/coeff_matrix/ with
 *                              two csf/<k>/ blocks; the only
 *                              input phase2 reads.
 *   n4_oss_k0_expected.txt  -- pre-computed dense JW-basis
 *                              amplitudes (idx<TAB>re<TAB>im
 *                              per line; 17-digit Float64
 *                              precision; ascending idx).
 *                              Comment lines start with '#'.
 *
 * Both files are produced offline by an independent reference
 * implementation and committed verbatim; this test treats them
 * as opaque inputs.  Asserts the expansion matches the reference
 * to TOL absolute amplitude difference and that the recovered
 * trial state has the same support (no extra non-zero indices
 * above the sparsity floor).
 */

#define WD_SEED UINT64_C(0x8f2a73bd1c4e5066)

#define TOL             1.0e-12
#define SPARSITY_FLOOR  1.0e-12

#define MAX_EXPECTED    1024

struct ref_amp {
	uint64_t idx;
	double   re;
	double   im;
};

static int read_expected(const char *path,
	struct ref_amp *out, size_t cap, size_t *n_out)
{
	FILE *f = fopen(path, "r");
	if (!f) {
		fprintf(stderr, "open %s: %m\n", path);
		return -1;
	}
	size_t n = 0;
	char line[256];
	while (fgets(line, sizeof line, f)) {
		/* skip blank lines and comments */
		const char *p = line;
		while (*p == ' ' || *p == '\t')
			p++;
		if (*p == '\0' || *p == '\n' || *p == '#')
			continue;
		if (n == cap) {
			fprintf(stderr,
				"too many ref amps (cap=%zu)\n", cap);
			fclose(f);
			return -1;
		}
		uint64_t idx;
		double re, im;
		if (sscanf(p, "%" SCNu64 " %lg %lg",
			&idx, &re, &im) != 3) {
			fprintf(stderr,
				"malformed ref line: %s", line);
			fclose(f);
			return -1;
		}
		out[n].idx = idx;
		out[n].re  = re;
		out[n].im  = im;
		n++;
	}
	fclose(f);
	*n_out = n;
	return 0;
}

static void t_n4_oss_k0_expansion(void)
{
	data_id fid =
		data_open(PH2_TESTDIR "/data/bendazzoli/n4_oss_k0/n4_oss_k0.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open n4_oss_k0.h5");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.nqb, 8u);
	TEST_EQ(cm.n_sites,  4u);
	TEST_EQ(cm.n_alpha,  2u);
	TEST_EQ(cm.n_beta,   2u);
	TEST_EQ(cm.n_components, 2u);

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, cm.nqb), 0);
	qreg_zero(&reg);
	TEST_EQ(state_prep_coeff_expand_all(&reg, &cm), 0);

	struct ref_amp *ref = malloc(sizeof(*ref) * MAX_EXPECTED);
	TEST_ASSERT(ref != NULL, "malloc ref");
	size_t n_ref = 0;
	TEST_EQ(read_expected(
		PH2_TESTDIR "/data/bendazzoli/n4_oss_k0/n4_oss_k0_expected.txt",
		ref, MAX_EXPECTED, &n_ref), 0);
	TEST_ASSERT(n_ref > 0, "empty reference");

	/* Mark which indices the reference covers so the
	 * "no extra amplitudes" check below can detect bins
	 * that phase2 populated but the reference did not. */
	const uint64_t namp = UINT64_C(1) << cm.nqb;
	int *covered = calloc(namp, sizeof *covered);
	TEST_ASSERT(covered != NULL, "calloc covered");

	double worst_diff  = 0.0;
	uint64_t worst_idx = 0;
	for (size_t i = 0; i < n_ref; i++) {
		TEST_ASSERT(ref[i].idx < namp,
			"ref idx %" PRIu64 " >= namp %" PRIu64,
			ref[i].idx, namp);
		_Complex double z;
		qreg_getamp(&reg, ref[i].idx, &z);
		const double dre = creal(z) - ref[i].re;
		const double dim = cimag(z) - ref[i].im;
		const double diff = sqrt(dre * dre + dim * dim);
		if (diff > worst_diff) {
			worst_diff = diff;
			worst_idx  = ref[i].idx;
		}
		covered[ref[i].idx] = 1;
	}
	TEST_ASSERT(worst_diff <= TOL,
		"worst amp diff %.3e at idx %" PRIu64
		" > %.0e", worst_diff, worst_idx, TOL);

	/* Conversely: phase2 must not produce any non-zero
	 * amplitude at an index the reference did not list. */
	uint64_t extras = 0;
	double worst_extra = 0.0;
	for (uint64_t i = 0; i < namp; i++) {
		if (covered[i])
			continue;
		_Complex double z;
		qreg_getamp(&reg, i, &z);
		const double m = sqrt(
			creal(z) * creal(z) + cimag(z) * cimag(z));
		if (m > SPARSITY_FLOOR) {
			extras++;
			if (m > worst_extra)
				worst_extra = m;
		}
	}
	TEST_ASSERT(extras == 0,
		"phase2 produced %" PRIu64 " extra amplitudes "
		"(worst |amp| = %.3e) at indices the reference "
		"does not cover", extras, worst_extra);

	free(covered);
	free(ref);
	qreg_free(&reg);
	data_coeff_matrix_free(&cm);
	data_close(fid);

	printf("t-ref-bendazzoli: %zu reference amplitudes, "
	       "worst diff %.3e\n", n_ref, worst_diff);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	t_n4_oss_k0_expansion();
	world_free();
}
