#include "c23_compat.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x6f3119a7c8771d22)

#define F_N4_CLOSED PH2_TESTDIR "/data/N4_closed.h5"
#define F_N4_OPEN PH2_TESTDIR "/data/N4_open.h5"
#define F_N4_CSF PH2_TESTDIR "/data/N4_csf.h5"
#define F_N8_TAP PH2_TESTDIR "/data/N8_tapered.h5"
#define F_N4_BOTH PH2_TESTDIR "/data/N4_both.h5"
#define F_N4_CSF_EMPTY PH2_TESTDIR "/data/N4_csf_empty.h5"
#define F_MD PH2_TESTDIR "/data/case-d9f603dc.h5_solved"

static void t_dispatch(void)
{
	data_id fid;
	enum stprep_kind k;

	fid = data_open(F_N4_CLOSED);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open closed");
	TEST_EQ(data_state_prep_kind(fid, &k), 0);
	TEST_EQ((int)k, (int)STPREP_COEFF_MATRIX);
	data_close(fid);

	fid = data_open(F_MD);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open multidet ref");
	TEST_EQ(data_state_prep_kind(fid, &k), 0);
	TEST_EQ((int)k, (int)STPREP_MULTIDET);
	data_close(fid);

	fid = data_open(F_N4_BOTH);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open both");
	const int rc = data_state_prep_kind(fid, &k);
	TEST_ASSERT(rc < 0, "both-present must error, got %d", rc);
	data_close(fid);
}

static void t_closed(void)
{
	data_id fid = data_open(F_N4_CLOSED);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open closed");

	uint32_t nqb, ns, na, nb;
	int cs, tap;
	TEST_EQ(data_coeff_matrix_getnums(fid, &nqb, &ns, &na, &nb, &cs, &tap),
		0);
	TEST_EQ(nqb, 8u);
	TEST_EQ(ns, 4u);
	TEST_EQ(na, 2u);
	TEST_EQ(nb, 2u);
	TEST_EQ(cs, 1);
	TEST_EQ(tap, 0);

	double Ca[4 * 2];
	TEST_EQ(data_coeff_matrix_read(fid, Ca, NULL), 0);
	double sumsq = 0.0;
	for (int i = 0; i < 8; i++)
		sumsq += Ca[i] * Ca[i];
	TEST_ASSERT(sumsq > 0.1, "C_alpha looks empty (sumsq=%f)", sumsq);

	size_t n_comp;
	TEST_EQ(data_coeff_matrix_csf_count(fid, &n_comp), 0);
	TEST_EQ(n_comp, 0u);

	data_close(fid);
}

static void t_open(void)
{
	data_id fid = data_open(F_N4_OPEN);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open open-shell");

	uint32_t nqb, ns, na, nb;
	int cs, tap;
	TEST_EQ(data_coeff_matrix_getnums(fid, &nqb, &ns, &na, &nb, &cs, &tap),
		0);
	TEST_EQ(cs, 0);
	TEST_EQ(na, 2u);
	TEST_EQ(nb, 1u);

	double Ca[4 * 2];
	double Cb[4 * 1];
	TEST_EQ(data_coeff_matrix_read(fid, Ca, Cb), 0);

	/* Passing the wrong shape of NULL is a programmer error. */
	TEST_ASSERT(data_coeff_matrix_read(fid, Ca, NULL) < 0,
		"open-shell must require C_beta buffer");

	data_close(fid);
}

static void t_tapered(void)
{
	data_id fid = data_open(F_N8_TAP);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open tapered");

	uint32_t nqb, ns, na, nb;
	int cs, tap;
	TEST_EQ(data_coeff_matrix_getnums(fid, &nqb, &ns, &na, &nb, &cs, &tap),
		0);
	TEST_EQ(nqb, 14u);
	TEST_EQ(ns, 8u);
	TEST_EQ(tap, 1);

	data_close(fid);
}

static void t_csf(void)
{
	data_id fid = data_open(F_N4_CSF);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open csf");

	size_t n_comp;
	TEST_EQ(data_coeff_matrix_csf_count(fid, &n_comp), 0);
	TEST_EQ(n_comp, 2u);

	double cf0, cf1;
	double Ca0[4 * 2], Ca1[4 * 2];
	TEST_EQ(data_coeff_matrix_csf_read(fid, 0, &cf0, Ca0, NULL), 0);
	TEST_EQ(data_coeff_matrix_csf_read(fid, 1, &cf1, Ca1, NULL), 0);
	TEST_ASSERT(fabs(cf0 - 0.6) < 1e-15, "csf[0].coef=%f", cf0);
	TEST_ASSERT(fabs(cf1 - 0.8) < 1e-15, "csf[1].coef=%f", cf1);

	data_close(fid);
}

static void t_csf_empty(void)
{
	/*
	 * An explicit `csf/` subgroup with n_components=0 must be
	 * rejected so a misbuilt pak cannot run on a silently
	 * empty state.
	 */
	data_id fid = data_open(F_N4_CSF_EMPTY);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open csf-empty");

	size_t n_comp = 999;
	const int rc = data_coeff_matrix_csf_count(fid, &n_comp);
	TEST_ASSERT(rc < 0, "csf n_components=0 must error, got %d", rc);

	data_close(fid);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_dispatch();
	t_closed();
	t_open();
	t_tapered();
	t_csf();
	t_csf_empty();

	world_free();
}
