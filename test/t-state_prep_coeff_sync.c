/*
 * Backend-coverage regression test for the host/device
 * sync invariant in state_prep_coeff_expand_all and
 * state_prep_coeff_inner.  Targets the two sites where
 * state_prep_coeff bypasses qreg_setamp / qreg_getamp by
 * touching reg->amp directly for performance:
 *
 *   1.  state_prep_coeff_expand_all must leave the
 *       canonical state populated on every backend.
 *       Test: load a known fixture, expand, read back via
 *       the backend-aware qreg_getamp, compare to the
 *       reference amplitudes.  On a CUDA backend without
 *       the sync this returns the all-zero device buffer
 *       and the test fails on the first non-zero ref.
 *
 *   2.  state_prep_coeff_inner must see the canonical
 *       state, including any post-state-prep device
 *       evolution.  Test: state-prep, apply one Pauli
 *       rotation (qreg_paulirot runs on-device on CUDA),
 *       call state_prep_coeff_inner against the same C
 *       matrices, expect a value distinct from the
 *       trivial <psi|psi>=1.  Without the device->host
 *       sync at the top of inner, reg->amp is the
 *       pre-rotation host shadow and the assertion fails.
 *
 * Both checks pass trivially on the CPU backend (where
 * reg->amp IS the canonical state) and pass on the CUDA
 * backend only with the bulk sync in place.
 */

#include "c23_compat.h"
#include <complex.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7c2bd7af09c4e51e)
#define TOL     1.0e-12

#define MAX_EXPECTED 1024

struct ref_amp {
	uint64_t idx;
	double   re;
	double   im;
};

static int read_expected(const char *path, struct ref_amp *out,
	size_t cap, size_t *n_out)
{
	FILE *f = fopen(path, "r");
	if (!f)
		return -1;
	size_t n = 0;
	char line[256];
	while (fgets(line, sizeof line, f)) {
		const char *p = line;
		while (*p == ' ' || *p == '\t')
			p++;
		if (*p == '\0' || *p == '\n' || *p == '#')
			continue;
		if (n == cap) {
			fclose(f);
			return -1;
		}
		uint64_t idx;
		double re, im;
		if (sscanf(p, "%" SCNu64 " %lg %lg",
			    &idx, &re, &im) != 3) {
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

/*
 * After state_prep_coeff_expand_all the backend-aware
 * qreg_getamp must return the chemist's trial-state
 * amplitudes -- whether the canonical state lives in
 * reg->amp (CPU) or cu->damp (CUDA).
 */
static void t_expand_all_visible_via_qreg_getamp(void)
{
	data_id fid = data_open(
		PH2_TESTDIR "/data/bendazzoli/n4_oss_k0/n4_oss_k0.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open n4_oss_k0.h5");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, cm.nqb), 0);
	struct state_prep_coeff_scratch sc;
	TEST_EQ(state_prep_coeff_scratch_init(&sc, cm.n_sites,
		cm.n_alpha, cm.n_beta), 0);
	TEST_EQ(state_prep_coeff_expand_all(&reg, &sc, &cm), 0);

	struct ref_amp *ref = malloc(sizeof(*ref) * MAX_EXPECTED);
	TEST_ASSERT(ref != NULL, "malloc ref");
	size_t n_ref = 0;
	TEST_EQ(read_expected(
		PH2_TESTDIR "/data/bendazzoli/n4_oss_k0/n4_oss_k0_expected.txt",
		ref, MAX_EXPECTED, &n_ref), 0);
	TEST_ASSERT(n_ref > 0, "empty reference");

	double worst = 0.0;
	for (size_t i = 0; i < n_ref; i++) {
		_Complex double z;
		qreg_getamp(&reg, ref[i].idx, &z);
		const double dre = creal(z) - ref[i].re;
		const double dim = cimag(z) - ref[i].im;
		const double d = sqrt(dre * dre + dim * dim);
		if (d > worst)
			worst = d;
	}
	TEST_ASSERT(worst <= TOL,
		"expand_all amp mismatch via qreg_getamp: %.3e > %.0e",
		worst, TOL);

	free(ref);
	state_prep_coeff_scratch_free(&sc);
	qreg_free(&reg);
	data_coeff_matrix_free(&cm);
	data_close(fid);

	printf("t-state_prep_coeff_sync: expand_all via "
	       "qreg_getamp -- worst diff %.3e\n", worst);
}

/*
 * After state-prep, an on-device kernel (qreg_paulirot)
 * evolves the state.  state_prep_coeff_inner must see the
 * evolved state, not the pre-rotation host shadow.
 *
 * Rotation: exp(i theta Z_0) on the n4_oss_k0 trial
 * state.  Z_0 is diagonal in the JW computational basis,
 * so the rotation imprints a per-amplitude phase of
 *   exp(+- i theta) depending on whether bit 0 of the
 * basis state is 0 or 1.  The inner product
 *   <psi | exp(i theta Z_0) | psi> = sum_i |c_i|^2 * e^(+- i theta)
 * has imaginary part that is non-zero for generic theta
 * (only zero if the trial state is supported entirely on
 * even-bit-0 or entirely on odd-bit-0 basis states, which
 * the OSS K=0 trial is not).
 *
 * Pre-fix CUDA: inner reads the stale pre-rotation
 * reg->amp; returns the unrotated <psi|psi> = 1 (real,
 * zero imaginary part).
 * Post-fix CUDA: inner syncs device->host, sees the
 * rotated state, returns a complex number with
 * non-trivial imaginary part.
 */
static void t_inner_after_paulirot(void)
{
	data_id fid = data_open(
		PH2_TESTDIR "/data/bendazzoli/n4_oss_k0/n4_oss_k0.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open n4_oss_k0.h5");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, cm.nqb), 0);
	struct state_prep_coeff_scratch sc;
	TEST_EQ(state_prep_coeff_scratch_init(&sc, cm.n_sites,
		cm.n_alpha, cm.n_beta), 0);
	TEST_EQ(state_prep_coeff_expand_all(&reg, &sc, &cm), 0);

	/* Sanity: pre-rotation <psi|psi> = 1.  This passes
	 * even without the device->host sync because the
	 * host shadow is in sync with the device immediately
	 * after expand_all (no intervening kernel). */
	_Complex double inner_pre = 0.0;
	TEST_EQ(state_prep_coeff_inner(&reg, &sc,
		cm.n_components ? cm.blocks[0].C_alpha : cm.C_alpha,
		NULL, 1.0, cm.tapered, &inner_pre), 0);

	/* exp(i theta Z_0): single-term paulirot with code
	 * representing Z on qubit 0, identity on the rest.
	 * The on-CUDA-backend execution runs on the device
	 * buffer; reg->amp becomes stale. */
	struct paulis code_hi = paulis_new();
	struct paulis code_lo = paulis_new();
	paulis_set(&code_lo, PAULI_Z, 0);
	const double theta = 0.37;
	qreg_paulirot(&reg, code_hi, &code_lo, &theta, 1);

	_Complex double inner_post = 0.0;
	TEST_EQ(state_prep_coeff_inner(&reg, &sc,
		cm.n_components ? cm.blocks[0].C_alpha : cm.C_alpha,
		NULL, 1.0, cm.tapered, &inner_post), 0);

	/* The rotated inner product must differ from the
	 * unrotated one.  The pre-fix CUDA bug returns the
	 * unrotated value (inner reads stale reg->amp); the
	 * delta would be exactly zero and the assertion
	 * fires. */
	const double dre = creal(inner_post) - creal(inner_pre);
	const double dim = cimag(inner_post) - cimag(inner_pre);
	const double delta = sqrt(dre * dre + dim * dim);
	TEST_ASSERT(delta > 1.0e-3,
		"inner unchanged by qreg_paulirot: pre=(%.6f,%.6f) "
		"post=(%.6f,%.6f) delta=%.3e -- the device->host "
		"sync at the top of state_prep_coeff_inner is "
		"missing",
		creal(inner_pre), cimag(inner_pre),
		creal(inner_post), cimag(inner_post), delta);

	state_prep_coeff_scratch_free(&sc);
	qreg_free(&reg);
	data_coeff_matrix_free(&cm);
	data_close(fid);

	printf("t-state_prep_coeff_sync: inner after paulirot "
	       "-- delta %.3e (theta=%.2f)\n", delta, theta);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	t_expand_all_visible_via_qreg_getamp();
	t_inner_after_paulirot();
	world_free();
}
