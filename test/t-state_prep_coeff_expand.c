#include "c23_compat.h"
#include <complex.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "combinations.h"
#include "det_small.h"
#include "phase2/data.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xae3c8829b51d04ee)

/*
 * Reference expansion: naive double loop, builds an in-memory
 * map of (idx -> amplitude).  Only used in tests; not on the
 * hot path.
 */
struct ref_amp {
	uint64_t idx;
	double cf;
};

static int cmp_amp(const void *a, const void *b)
{
	const struct ref_amp *x = a;
	const struct ref_amp *y = b;
	if (x->idx < y->idx)
		return -1;
	if (x->idx > y->idx)
		return 1;
	return 0;
}

static size_t expand_ref(uint32_t n_sites, uint32_t na, uint32_t nb,
	const double *Ca, const double *Cb, int tapered, struct ref_amp *out,
	size_t cap)
{
	struct combo it_a;
	combinations_init(&it_a, n_sites, na);
	uint32_t ta[COMBINATIONS_MAX_K];
	uint32_t tb[COMBINATIONS_MAX_K];

	size_t n = 0;

	const double *Cbb = Cb ? Cb : Ca;

	while (combinations_next(&it_a, ta) == 0) {
		double sub_a[DET_SMALL_MAX_N * DET_SMALL_MAX_N];
		double da;
		if (na == 0) {
			da = 1.0;
		} else {
			for (uint32_t r = 0; r < na; r++)
				for (uint32_t c = 0; c < na; c++)
					sub_a[r * na + c] =
						Ca[ta[r] * na + c];
			da = det_small(sub_a, na);
		}
		struct combo it_b;
		combinations_init(&it_b, n_sites, nb);
		while (combinations_next(&it_b, tb) == 0) {
			double sub_b[DET_SMALL_MAX_N * DET_SMALL_MAX_N];
			double db;
			if (nb == 0) {
				db = 1.0;
			} else {
				for (uint32_t r = 0; r < nb; r++)
					for (uint32_t c = 0; c < nb; c++)
						sub_b[r * nb + c] =
							Cbb[tb[r] * nb + c];
				db = det_small(sub_b, nb);
			}
			const double cf = da * db;
			if (fabs(cf) < 1e-12)
				continue;
			uint64_t idx = 0;
			for (uint32_t i = 0; i < na; i++)
				idx |= UINT64_C(1) << ta[i];
			for (uint32_t i = 0; i < nb; i++)
				idx |= UINT64_C(1) << (n_sites + tb[i]);
			if (tapered) {
				uint64_t lo =
					(idx >> 1)
					& ((UINT64_C(1) << (n_sites - 1)) - 1);
				uint64_t hi = idx >> (n_sites + 1);
				idx = lo | (hi << (n_sites - 1));
			}
			if (n >= cap)
				TEST_FAIL("ref cap exceeded");
			out[n].idx = idx;
			out[n].cf = cf;
			n++;
		}
	}
	return n;
}

/*
 * Open a fixture, expand into a fresh qreg, compare every
 * amplitude to the reference expansion.
 */
static void t_fixture(const char *path)
{
	data_id fid = data_open(path);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open %s", path);

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.n_components, (size_t)0);

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, cm.nqb), 0);
	qreg_zero(&reg);

	TEST_EQ(state_prep_coeff_expand(&reg, cm.n_sites, cm.n_alpha, cm.n_beta,
			cm.C_alpha, cm.C_beta, 1.0, cm.tapered, 0),
		0);

	struct ref_amp *refs =
		malloc(sizeof(struct ref_amp) * (size_t)1 << 18);
	const size_t n_ref = expand_ref(cm.n_sites, cm.n_alpha, cm.n_beta,
		cm.C_alpha, cm.closed_shell ? NULL : cm.C_beta, cm.tapered,
		refs, (size_t)1 << 18);
	qsort(refs, n_ref, sizeof refs[0], cmp_amp);

	double total_diff = 0.0;
	const uint64_t namp = UINT64_C(1) << cm.nqb;
	for (uint64_t idx = 0; idx < namp; idx++) {
		_Complex double z;
		qreg_getamp(&reg, idx, &z);
		double expect = 0.0;
		for (size_t i = 0; i < n_ref; i++)
			if (refs[i].idx == idx) {
				expect = refs[i].cf;
				break;
			}
		const double d = fabs(creal(z) - expect)
				 + fabs(cimag(z));
		total_diff += d;
		TEST_ASSERT(d < 1e-12,
			"%s idx=%lu got=%f+%fi expect=%f", path,
			(unsigned long)idx, creal(z), cimag(z), expect);
	}
	(void)total_diff;

	free(refs);
	data_coeff_matrix_free(&cm);
	qreg_free(&reg);
	data_close(fid);
}

static void t_identity(void)
{
	const uint32_t n_sites = 4, na = 2, nb = 2;
	double Ca[4 * 2] = { 1, 0, 0, 1, 0, 0, 0, 0 };
	struct qreg reg;
	TEST_EQ(qreg_init(&reg, 2 * n_sites), 0);
	qreg_zero(&reg);

	TEST_EQ(state_prep_coeff_expand(&reg, n_sites, na, nb, Ca, NULL, 1.0,
			0, 0),
		0);

	const uint64_t hf = (1u << 0) | (1u << 1) | (1u << (n_sites + 0))
			    | (1u << (n_sites + 1));
	const uint64_t namp = UINT64_C(1) << (2 * n_sites);
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(&reg, i, &z);
		if (i == hf)
			TEST_ASSERT(fabs(creal(z) - 1.0) < 1e-14,
				"HF amp got %f+%fi", creal(z), cimag(z));
		else
			TEST_ASSERT(cabs(z) < 1e-14,
				"i=%lu got %f+%fi", (unsigned long)i,
				creal(z), cimag(z));
	}
	qreg_free(&reg);
}

static void t_boundary_zero(void)
{
	/* n_alpha = 0, n_beta = 0 -> the vacuum state. */
	const uint32_t n_sites = 4;
	struct qreg reg;
	TEST_EQ(qreg_init(&reg, 2 * n_sites), 0);
	qreg_zero(&reg);
	TEST_EQ(state_prep_coeff_expand(
			&reg, n_sites, 0, 0, NULL, NULL, 1.0, 0, 0),
		0);
	_Complex double z;
	qreg_getamp(&reg, 0, &z);
	TEST_ASSERT(fabs(creal(z) - 1.0) < 1e-14, "vacuum amp[0]=%f", creal(z));
	qreg_free(&reg);
}

static void t_boundary_full(void)
{
	/* n_alpha = n_sites = 4 -> the |1111> alpha block. */
	const uint32_t n_sites = 4;
	double Ca[4 * 4] = {
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
	};
	struct qreg reg;
	TEST_EQ(qreg_init(&reg, 2 * n_sites), 0);
	qreg_zero(&reg);
	TEST_EQ(state_prep_coeff_expand(
			&reg, n_sites, n_sites, 0, Ca, NULL, 1.0, 0, 0),
		0);
	const uint64_t target = 0xfu;
	_Complex double z;
	qreg_getamp(&reg, target, &z);
	TEST_ASSERT(fabs(creal(z) - 1.0) < 1e-14, "full amp got %f+%fi",
		creal(z), cimag(z));
	qreg_free(&reg);
}

/*
 * Dump every non-zero amplitude of the expanded register to
 * `out` as `idx re im` lines, one per line, sorted by idx
 * ascending.  Only rank 0 writes; other ranks remain silent
 * (qreg_getamp broadcasts so rank 0 sees all owned slots).
 *
 * Used by `t-ref-coeff_matrix.py` for the strict bit-for-bit
 * cross-validation against the Python reference at
 * `test/ref/coeff_matrix_reference.py`.
 */
static int dump_expand(const char *fixture, FILE *out)
{
	data_id fid = data_open(fixture);
	if (fid == DATA_INVALID_FID) {
		fprintf(stderr, "dump_expand: cannot open %s\n", fixture);
		return -1;
	}

	/* `data_coeff_matrix_load` reads both top-level C matrices
	 * and any csf/<k>/ subgroups via the same path
	 * `circ_prepst` uses; we dispatch through
	 * `state_prep_coeff_expand_all` so the dump covers
	 * single-block AND CSF fixtures with one code path. */
	struct data_coeff_matrix cm;
	if (data_coeff_matrix_load(fid, &cm) < 0) {
		fprintf(stderr, "dump_expand: data_coeff_matrix_load failed\n");
		data_close(fid);
		return -1;
	}

	struct qreg reg;
	if (qreg_init(&reg, cm.nqb) < 0) {
		data_coeff_matrix_free(&cm);
		data_close(fid);
		return -1;
	}
	qreg_zero(&reg);

	if (state_prep_coeff_expand_all(&reg, &cm) < 0) {
		qreg_free(&reg);
		data_coeff_matrix_free(&cm);
		data_close(fid);
		return -1;
	}

	const uint64_t namp = UINT64_C(1) << cm.nqb;
	for (uint64_t idx = 0; idx < namp; idx++) {
		_Complex double z;
		qreg_getamp(&reg, idx, &z);
		if (reg.wd.rank != 0)
			continue;
		const double re = creal(z);
		const double im = cimag(z);
		if (fabs(re) + fabs(im) == 0.0)
			continue;
		fprintf(out, "%" PRIu64 " %.17g %.17g\n",
			idx, re, im);
	}
	if (reg.wd.rank == 0)
		fflush(out);

	qreg_free(&reg);
	data_coeff_matrix_free(&cm);
	data_close(fid);
	return 0;
}

int main(int argc, char **argv)
{
	world_init(&argc, &argv, WD_SEED);

	/*
	 * --dump <out-path> <fixture-h5>: emit the expanded
	 * (idx, re, im) triples for the C side, then exit.  Used
	 * by the strict cross-validation harness.
	 */
	if (argc == 4 && strcmp(argv[1], "--dump") == 0) {
		const char *out_path = argv[2];
		const char *fixture = argv[3];
		FILE *out = NULL;
		int rc = 0;
		if (strcmp(out_path, "-") == 0) {
			out = stdout;
		} else {
			out = fopen(out_path, "w");
			if (!out) {
				fprintf(stderr, "cannot open %s for write\n",
					out_path);
				world_free();
				return 2;
			}
		}
		rc = dump_expand(fixture, out);
		if (out != stdout)
			fclose(out);
		world_free();
		return rc == 0 ? 0 : 1;
	}

	t_identity();
	t_boundary_zero();
	t_boundary_full();

	t_fixture(PH2_TESTDIR "/data/N4_closed.h5");
	t_fixture(PH2_TESTDIR "/data/N4_open.h5");
	t_fixture(PH2_TESTDIR "/data/N8_untapered.h5");
	t_fixture(PH2_TESTDIR "/data/N8_tapered.h5");

	world_free();
}
