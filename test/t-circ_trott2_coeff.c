#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ/trott2.h"
#include "log.h"
#include "ph2run/data.h"
#include "phase2.h"

#include "t-data.h"
#include "test.h"

#define WD_SEED UINT64_C(0xd1bc44e2918887aa)
#define MARGIN (1.0e-13)
#define STEPS (2)

static _Complex double run_trott2(const char *src)
{
	const char *tmp = "/tmp/t-circ_trott2_coeff.tmp.h5";
	test_fixture_copy(src, tmp);

	struct trott2 tt;
	struct trott2_data td = { .delta = 0.35, .steps = STEPS };

	data_id fid = data_open(tmp);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open %s", tmp);

	struct circ_hamil hm;
	TEST_EQ(data_hamil_load(fid, &hm), 0);
	enum stprep_kind k;
	TEST_EQ(data_state_prep_kind(fid, &k), 0);
	struct circ_muldet md = { 0 };
	struct data_coeff_matrix cm = { 0 };
	const void *sp = nullptr;
	switch (k) {
	case STPREP_MULTIDET:
		TEST_EQ(data_muldet_load(fid, &md), 0);
		sp = &md;
		break;
	case STPREP_COEFF_MATRIX:
		TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
		sp = &cm;
		break;
	}
	TEST_EQ(trott2_init(&tt, &td, hm, k, sp, nullptr), 0);
	TEST_EQ(trott2_simul(&tt), 0);
	const _Complex double v = tt.ct.vals.z[STEPS - 1];
	trott2_free(&tt);
	data_close(fid);

	struct world_info wd;
	world_info(&wd);
	if (wd.rank == 0)
		remove(tmp);
	return v;
}

static void t_n4_equivalence(void)
{
	const _Complex double v_cm = run_trott2(
		PH2_TESTDIR "/data/N4_closed.h5");
	const _Complex double v_md = run_trott2(
		PH2_TESTDIR "/data/N4_multidet.h5");

	TEST_ASSERT(cabs(v_cm - v_md) < MARGIN,
		"trott2 eq: cm=%f+%fi md=%f+%fi diff=%g", creal(v_cm),
		cimag(v_cm), creal(v_md), cimag(v_md),
		cabs(v_cm - v_md));
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	t_n4_equivalence();
	world_free();
}
