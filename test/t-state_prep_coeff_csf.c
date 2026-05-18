#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xb2c3a98d4e110097)

static void t_two_component(void)
{
	data_id fid = data_open(PH2_TESTDIR "/data/N4_csf.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open csf");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.n_components, 2u);
	TEST_NEAR(cm.blocks[0].cf, 0.6, 1e-15);
	TEST_NEAR(cm.blocks[1].cf, 0.8, 1e-15);

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, cm.nqb), 0);
	qreg_zero(&reg);
	struct state_prep_coeff_scratch sc;
	TEST_EQ(state_prep_coeff_scratch_init(&sc, cm.n_sites,
		cm.n_alpha, cm.n_beta), 0);
	TEST_EQ(state_prep_coeff_expand_all(&reg, &sc, &cm), 0);
	state_prep_coeff_scratch_free(&sc);

	/* Sanity: the superposed state must have at least one
	 * non-zero amplitude. */
	double sumsq = 0.0;
	const uint64_t namp = UINT64_C(1) << cm.nqb;
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(&reg, i, &z);
		sumsq += creal(z) * creal(z) + cimag(z) * cimag(z);
	}
	TEST_ASSERT(sumsq > 0.01, "CSF state is zero (sumsq=%f)", sumsq);

	qreg_free(&reg);
	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_single_block_csf_equiv(void)
{
	/* state_prep_coeff_expand_all on n_components=0 must
	 * produce the same state as state_prep_coeff_expand
	 * with accumulate=0. */
	data_id fid = data_open(PH2_TESTDIR "/data/N4_closed.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open closed");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.n_components, 0u);

	struct qreg r0, r1;
	TEST_EQ(qreg_init(&r0, cm.nqb), 0);
	TEST_EQ(qreg_init(&r1, cm.nqb), 0);
	qreg_zero(&r0);
	qreg_zero(&r1);

	struct state_prep_coeff_scratch sc;
	TEST_EQ(state_prep_coeff_scratch_init(&sc, cm.n_sites,
		cm.n_alpha, cm.n_beta), 0);
	TEST_EQ(state_prep_coeff_expand_all(&r0, &sc, &cm), 0);
	TEST_EQ(state_prep_coeff_expand(&r1, &sc, cm.C_alpha, NULL,
		1.0, cm.tapered, 0), 0);
	state_prep_coeff_scratch_free(&sc);

	const uint64_t namp = UINT64_C(1) << cm.nqb;
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double a, b;
		qreg_getamp(&r0, i, &a);
		qreg_getamp(&r1, i, &b);
		TEST_ASSERT(cabs(a - b) < 1e-15,
			"wrapper transparency: i=%lu a=%f b=%f",
			(unsigned long)i, creal(a), creal(b));
	}

	qreg_free(&r1);
	qreg_free(&r0);
	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_csf_empty_rejected(void)
{
	/*
	 * A coeff_matrix fixture whose csf/ subgroup exists but
	 * advertises n_components=0 is malformed.
	 * data_coeff_matrix_load must refuse to populate the
	 * struct.
	 */
	data_id fid = data_open(PH2_TESTDIR "/data/N4_csf_empty.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open csf-empty");

	struct data_coeff_matrix cm;
	const int rc = data_coeff_matrix_load(fid, &cm);
	TEST_ASSERT(rc < 0,
		"data_coeff_matrix_load must reject empty csf/, got %d", rc);

	data_close(fid);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_two_component();
	t_single_block_csf_equiv();
	t_csf_empty_rejected();

	world_free();
}
