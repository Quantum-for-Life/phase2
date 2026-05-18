#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "ph2run/data.h"
#include "phase2.h"

#include "circ_cache.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7766aa11d8dfbc4c)

/* Build a fresh circ from the N4_closed fixture.  Caller
 * owns ct and the open data_id; must call circ_free +
 * data_close. */
static data_id open_n4_closed(struct circ *ct)
{
	data_id fid = data_open(PH2_TESTDIR "/data/N4_closed.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open closed");

	struct circ_hamil hm;
	TEST_EQ(data_hamil_load(fid, &hm), 0);
	enum stprep_kind k;
	TEST_EQ(data_state_prep_kind(fid, &k), 0);
	TEST_EQ((int)k, (int)STPREP_COEFF_MATRIX);
	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);

	memset(ct, 0, sizeof *ct);
	TEST_EQ(circ_init(ct, hm, k, &cm, 1), 0);
	return fid;
}

static double sumsq_state(struct qreg *reg, uint32_t nqb)
{
	double s = 0.0;
	const uint64_t namp = UINT64_C(1) << nqb;
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(reg, i, &z);
		s += creal(z) * creal(z) + cimag(z) * cimag(z);
	}
	return s;
}

static void t_norm_and_structure(void)
{
	struct circ ct;
	data_id fid = open_n4_closed(&ct);

	TEST_EQ(circ_prepst(&ct), 0);

	const uint64_t namp = UINT64_C(1) << ct.cm.nqb;
	double sumsq = 0.0;
	int nonzero = 0;
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(&ct.reg, i, &z);
		const double m2 = creal(z) * creal(z) + cimag(z) * cimag(z);
		sumsq += m2;
		if (m2 > 1e-24)
			nonzero++;
	}
	TEST_ASSERT(nonzero > 0, "no non-zero amplitudes");
	TEST_ASSERT(sumsq > 1e-6,
		"state effectively zero (sumsq=%f)", sumsq);

	circ_free(&ct);
	data_close(fid);
}

/*
 * Two back-to-back circ_prepst calls against the same
 * cached scratch on struct circ must produce bit-
 * identical amplitudes.
 */
static void t_prepst_idempotent(void)
{
	struct circ ct;
	data_id fid = open_n4_closed(&ct);

	const uint64_t namp = UINT64_C(1) << ct.cm.nqb;
	_Complex double *snap = malloc(sizeof *snap * namp);
	TEST_ASSERT(snap != NULL, "alloc snap");

	TEST_EQ(circ_prepst(&ct), 0);
	for (uint64_t i = 0; i < namp; i++)
		qreg_getamp(&ct.reg, i, &snap[i]);

	TEST_EQ(circ_prepst(&ct), 0);
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(&ct.reg, i, &z);
		TEST_ASSERT(z == snap[i],
			"prepst idempotent: i=%lu first=%g+%gi"
			" second=%g+%gi", (unsigned long)i,
			creal(snap[i]), cimag(snap[i]),
			creal(z), cimag(z));
	}

	free(snap);
	circ_free(&ct);
	data_close(fid);
}

/*
 * Evolution between prepst calls must not contaminate
 * the second prepst.  Gates the expand_all-zeros-itself
 * contract from commit 4c4be70.
 */
static void t_prepst_after_evolution(void)
{
	struct circ ct;
	data_id fid = open_n4_closed(&ct);

	const uint64_t namp = UINT64_C(1) << ct.cm.nqb;
	_Complex double *snap = malloc(sizeof *snap * namp);
	TEST_ASSERT(snap != NULL, "alloc snap");

	TEST_EQ(circ_prepst(&ct), 0);
	for (uint64_t i = 0; i < namp; i++)
		qreg_getamp(&ct.reg, i, &snap[i]);

	/* Apply one paulirot to dirty the register.  Use the
	 * first Hamiltonian term, split into hi / lo. */
	TEST_ASSERT(ct.hm.len > 0, "Hamiltonian has terms");
	struct paulis lo, hi;
	paulis_split(ct.hm.terms[0].op, ct.reg.qb_lo, ct.reg.qb_hi,
		&lo, &hi);
	const double phi = 0.3;
	qreg_paulirot(&ct.reg, hi, &lo, &phi, 1);

	TEST_EQ(circ_prepst(&ct), 0);
	for (uint64_t i = 0; i < namp; i++) {
		_Complex double z;
		qreg_getamp(&ct.reg, i, &z);
		TEST_ASSERT(z == snap[i],
			"prepst after evolution: i=%lu before=%g+%gi"
			" after=%g+%gi", (unsigned long)i,
			creal(snap[i]), cimag(snap[i]),
			creal(z), cimag(z));
	}

	free(snap);
	circ_free(&ct);
	data_close(fid);
}

/*
 * After prepst, the trial-state inner product matches
 * the squared norm of the register (since the trial
 * state IS the register at this point).
 */
static void t_prepst_then_measure(void)
{
	struct circ ct;
	data_id fid = open_n4_closed(&ct);

	TEST_EQ(circ_prepst(&ct), 0);
	const double sumsq = sumsq_state(&ct.reg, ct.cm.nqb);

	const _Complex double inner = circ_measure(&ct);
	TEST_NEAR(creal(inner), sumsq, 1e-12);
	TEST_NEAR(cimag(inner), 0.0, 1e-12);

	circ_free(&ct);
	data_close(fid);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	t_norm_and_structure();
	t_prepst_idempotent();
	t_prepst_after_evolution();
	t_prepst_then_measure();
	world_free();
}
