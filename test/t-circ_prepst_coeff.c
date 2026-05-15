#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "log.h"
#include "ph2run/data.h"
#include "phase2.h"

#include "circ_cache.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7766aa11d8dfbc4c)

static void t_norm_and_structure(void)
{
	data_id fid = data_open(PH2_TESTDIR "/data/N4_closed.h5");
	TEST_ASSERT(fid != DATA_INVALID_FID, "open closed");

	struct circ_hamil hm;
	TEST_EQ(circ_hamil_load(fid, &hm), 0);
	enum stprep_kind k;
	TEST_EQ(data_state_prep_kind(fid, &k), 0);
	TEST_EQ((int)k, (int)STPREP_COEFF_MATRIX);
	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);

	struct circ ct = { 0 };
	TEST_EQ(circ_init(&ct, hm, k, &cm, 1), 0);

	TEST_EQ(circ_prepst(&ct), 0);

	double sumsq = 0.0;
	int nonzero = 0;
	const uint64_t namp = UINT64_C(1) << ct.cm.nqb;
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

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	t_norm_and_structure();
	world_free();
}
