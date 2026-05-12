#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ/trott2.h"
#include "log.h"
#include "phase2.h"

#include "test.h"

#define WD_SEED UINT64_C(0xd1bc44e2918887aa)
#define MARGIN (1.0e-13)
#define STEPS (2)

static _Complex double run_trott2(const char *src)
{
	struct trott2 tt;
	struct trott2_data td = { .delta = 0.35, .steps = STEPS };

	data_id fid = data_open(src);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open %s", src);
	TEST_EQ(trott2_init(&tt, &td, fid), 0);
	TEST_EQ(trott2_simul(&tt), 0);
	const _Complex double v = tt.ct.vals.z[STEPS - 1];
	trott2_free(&tt);
	data_close(fid);
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
