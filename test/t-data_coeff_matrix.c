#include "c23_compat.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/circ.h"
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

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.nqb, 8u);
	TEST_EQ(cm.n_sites, 4u);
	TEST_EQ(cm.n_alpha, 2u);
	TEST_EQ(cm.n_beta, 2u);
	TEST_EQ(cm.closed_shell, 1);
	TEST_EQ(cm.tapered, 0);
	TEST_EQ(cm.n_components, (size_t)0);
	TEST_ASSERT(cm.C_alpha != NULL, "C_alpha must be allocated");
	TEST_ASSERT(cm.C_beta == NULL, "closed shell: C_beta must be NULL");
	TEST_ASSERT(cm.blocks == NULL, "single block: blocks must be NULL");
	double sumsq = 0.0;
	for (size_t i = 0; i < (size_t)cm.n_sites * cm.n_alpha; i++)
		sumsq += cm.C_alpha[i] * cm.C_alpha[i];
	TEST_ASSERT(sumsq > 0.1, "C_alpha looks empty (sumsq=%f)", sumsq);

	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_open(void)
{
	data_id fid = data_open(F_N4_OPEN);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open open-shell");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.closed_shell, 0);
	TEST_EQ(cm.n_alpha, 2u);
	TEST_EQ(cm.n_beta, 1u);
	TEST_EQ(cm.n_components, (size_t)0);
	TEST_ASSERT(cm.C_alpha != NULL, "C_alpha must be allocated");
	TEST_ASSERT(cm.C_beta != NULL, "open shell: C_beta must be allocated");

	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_tapered(void)
{
	data_id fid = data_open(F_N8_TAP);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open tapered");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.nqb, 14u);
	TEST_EQ(cm.n_sites, 8u);
	TEST_EQ(cm.tapered, 1);

	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_csf(void)
{
	data_id fid = data_open(F_N4_CSF);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open csf");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.n_components, (size_t)2);
	TEST_ASSERT(cm.blocks != NULL, "csf: blocks must be allocated");
	TEST_ASSERT(cm.C_alpha == NULL,
		"csf: top-level C_alpha must be NULL (use blocks[k])");
	TEST_ASSERT(cm.C_beta == NULL,
		"csf: top-level C_beta must be NULL (use blocks[k])");
	TEST_ASSERT(fabs(cm.blocks[0].cf - 0.6) < 1e-15, "csf[0].cf=%f",
		cm.blocks[0].cf);
	TEST_ASSERT(fabs(cm.blocks[1].cf - 0.8) < 1e-15, "csf[1].cf=%f",
		cm.blocks[1].cf);
	TEST_ASSERT(cm.blocks[0].C_alpha != NULL, "csf[0].C_alpha allocated");
	TEST_ASSERT(cm.blocks[1].C_alpha != NULL, "csf[1].C_alpha allocated");

	data_coeff_matrix_free(&cm);
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

	struct data_coeff_matrix cm;
	const int rc = data_coeff_matrix_load(fid, &cm);
	TEST_ASSERT(rc < 0, "csf n_components=0 must error, got %d", rc);
	TEST_ASSERT(cm.C_alpha == NULL && cm.blocks == NULL,
		"failed load must leave the struct empty");

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
