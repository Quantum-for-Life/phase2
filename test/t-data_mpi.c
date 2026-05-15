/*
 * t-data_mpi -- MPI-specific behaviour of the data layer:
 *
 *   1. data_open returns a real fid on rank 0 and
 *      DATA_FOLLOWER_FID on other ranks.
 *   2. The collective Hamiltonian read via data_hamil_load
 *      produces byte-equal packed terms on every rank
 *      (Bcast worked).  Same property is locked down for
 *      data_muldet_load and data_coeff_matrix_load: a
 *      second bcast of rank-0's locally-loaded values is
 *      compared element-by-element against each follower's
 *      load.
 *   3. data_circ_writer_init pre-allocates the values dataset
 *      with NaN padding and caches the open handle;
 *      data_circ_write_step fills rows one at a time --
 *      after writing rows 0..k-1, the remaining rows
 *      0..N-1 are still NaN.  Rank 0 is the only writer.
 *
 * Run via `make check-mpi MPIRANKS>=2`.  Single-rank
 * invocation also passes (the rank-0 vs follower scenarios
 * fold into rank-0 paths).
 */
#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <hdf5.h>

#include "mpi.h"

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/paulis.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x21847b9ef0d3c2a4)

#define N_TERMS (5)
#define N_QB (3)
#define N_STEPS (8)
#define N_WRITTEN (4)

#define N_DETS (2)
#define MD_N_QB (3)

#define CM_N_QB (4)
#define CM_N_SITES (2)
#define CM_N_ALPHA (1)
#define CM_N_BETA (1)

static char *FILENAME = "/tmp/t-data_mpi.h5";
static char *FILENAME_MD = "/tmp/t-data_mpi-md.h5";
static char *FILENAME_CM = "/tmp/t-data_mpi-cm.h5";

/* Rank 0 builds a tiny pauli_hamil group; the test then
 * reads it through data_hamil_load on every rank and
 * confirms every rank ends up with the same packed terms. */
static int build_hamil_fixture(void)
{
	const hid_t f = H5Fcreate(
		FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (f == H5I_INVALID_HID)
		return -1;
	const hid_t g = H5Gcreate(f, "pauli_hamil", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);

	const hid_t asp = H5Screate(H5S_SCALAR);
	const hid_t aid = H5Acreate2(g, "normalization", H5T_IEEE_F64LE, asp,
		H5P_DEFAULT, H5P_DEFAULT);
	const double anv = 1.0;
	H5Awrite(aid, H5T_NATIVE_DOUBLE, &anv);
	H5Aclose(aid);
	H5Sclose(asp);

	const hsize_t cd[1] = { N_TERMS };
	const hid_t csp = H5Screate_simple(1, cd, NULL);
	const hid_t cds = H5Dcreate2(g, "coeffs", H5T_IEEE_F64LE, csp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const double cv[N_TERMS] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
	H5Dwrite(cds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cv);
	H5Dclose(cds);
	H5Sclose(csp);

	const hsize_t pd[2] = { N_TERMS, N_QB };
	const hid_t psp = H5Screate_simple(2, pd, NULL);
	const hid_t pds = H5Dcreate2(g, "paulis", H5T_NATIVE_UCHAR, psp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const unsigned char pv[N_TERMS * N_QB] = {
		1, 0, 0,
		0, 2, 0,
		0, 0, 3,
		1, 2, 0,
		2, 0, 3,
	};
	H5Dwrite(pds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, pv);
	H5Dclose(pds);
	H5Sclose(psp);

	H5Gclose(g);
	H5Fclose(f);
	return 0;
}

/* Rank 0 builds a tiny /state_prep/multidet group; the test
 * then reads it through data_muldet_load on every rank and
 * confirms every rank ends up with the same packed dets. */
static int build_muldet_fixture(void)
{
	const hid_t f = H5Fcreate(
		FILENAME_MD, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (f == H5I_INVALID_HID)
		return -1;
	const hid_t sp = H5Gcreate(f, "state_prep", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);
	const hid_t g = H5Gcreate(sp, "multidet", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);

	const hsize_t cd[2] = { N_DETS, 2 };
	const hid_t csp = H5Screate_simple(2, cd, NULL);
	const hid_t cds = H5Dcreate2(g, "coeffs", H5T_IEEE_F64LE, csp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const double cv[N_DETS * 2] = {
		0.6, 0.1,
		0.8, -0.2,
	};
	H5Dwrite(cds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cv);
	H5Dclose(cds);
	H5Sclose(csp);

	const hsize_t dd[2] = { N_DETS, MD_N_QB };
	const hid_t dsp = H5Screate_simple(2, dd, NULL);
	const hid_t dds = H5Dcreate2(g, "dets", H5T_NATIVE_UCHAR, dsp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const unsigned char dv[N_DETS * MD_N_QB] = {
		1, 0, 1,
		0, 1, 1,
	};
	H5Dwrite(dds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, dv);
	H5Dclose(dds);
	H5Sclose(dsp);

	H5Gclose(g);
	H5Gclose(sp);
	H5Fclose(f);
	return 0;
}

/* Helper: write a scalar uint32 attribute to grp. */
static void write_u32_attr(hid_t grp, const char *name, uint32_t v)
{
	const hid_t sp = H5Screate(H5S_SCALAR);
	const hid_t aid = H5Acreate2(grp, name, H5T_STD_U32LE, sp,
		H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, H5T_NATIVE_UINT32, &v);
	H5Aclose(aid);
	H5Sclose(sp);
}

/* Helper: write a scalar uint8 attribute to grp. */
static void write_u8_attr(hid_t grp, const char *name, uint8_t v)
{
	const hid_t sp = H5Screate(H5S_SCALAR);
	const hid_t aid = H5Acreate2(grp, name, H5T_STD_U8LE, sp,
		H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, H5T_NATIVE_UINT8, &v);
	H5Aclose(aid);
	H5Sclose(sp);
}

/* Rank 0 builds a tiny single-block /state_prep/coeff_matrix
 * group; the test then reads it through data_coeff_matrix_load
 * on every rank and confirms every rank ends up with the same
 * scalar fields and C_alpha values. */
static int build_coeffmat_fixture(void)
{
	const hid_t f = H5Fcreate(
		FILENAME_CM, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (f == H5I_INVALID_HID)
		return -1;
	const hid_t sp = H5Gcreate(f, "state_prep", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);
	const hid_t g = H5Gcreate(sp, "coeff_matrix", H5P_DEFAULT,
		H5P_DEFAULT, H5P_DEFAULT);

	write_u32_attr(g, "n_qubits", CM_N_QB);
	write_u32_attr(g, "n_sites", CM_N_SITES);
	write_u32_attr(g, "n_alpha", CM_N_ALPHA);
	write_u32_attr(g, "n_beta", CM_N_BETA);
	write_u8_attr(g, "closed_shell", 1);
	write_u8_attr(g, "tapered", 0);

	const hsize_t cd[2] = { CM_N_SITES, CM_N_ALPHA };
	const hid_t csp = H5Screate_simple(2, cd, NULL);
	const hid_t cds = H5Dcreate2(g, "C_alpha", H5T_IEEE_F64LE, csp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const double cv[CM_N_SITES * CM_N_ALPHA] = {
		0.6, 0.8,
	};
	H5Dwrite(cds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cv);
	H5Dclose(cds);
	H5Sclose(csp);

	H5Gclose(g);
	H5Gclose(sp);
	H5Fclose(f);
	return 0;
}

static void t_open_follower_sentinel(int rank)
{
	const data_id fid = data_open(FILENAME);
	if (rank == 0) {
		TEST_ASSERT(fid > 0,
			"rank 0: expected positive fid, got %lld",
			(long long)fid);
	} else {
		TEST_EQ(fid, DATA_FOLLOWER_FID);
	}
	data_close(fid);
}

static void t_bcast_buffers_match(void)
{
	const data_id fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open for bcast test");

	struct circ_hamil hm;
	TEST_EQ(data_hamil_load(fid, &hm), 0);
	TEST_EQ(hm.qb, (uint32_t)N_QB);
	TEST_EQ(hm.len, (size_t)N_TERMS);

	/* Send rank-0's packed values to all ranks via a second
	 * bcast, compare against what each rank loaded. */
	double ref_cfs[N_TERMS];
	struct paulis ref_ops[N_TERMS];
	for (size_t i = 0; i < N_TERMS; i++) {
		ref_cfs[i] = hm.terms[i].cf;
		ref_ops[i] = hm.terms[i].op;
	}
	MPI_Bcast(ref_cfs, N_TERMS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ref_ops, N_TERMS * (int)sizeof(struct paulis),
		MPI_BYTE, 0, MPI_COMM_WORLD);
	for (size_t i = 0; i < N_TERMS; i++) {
		TEST_ASSERT(ref_cfs[i] == hm.terms[i].cf,
			"term[%zu].cf differs between rank 0 and this rank",
			i);
		TEST_ASSERT(paulis_eq(ref_ops[i], hm.terms[i].op),
			"term[%zu].op differs between rank 0 and this rank",
			i);
	}

	circ_hamil_free(&hm);
	data_close(fid);
}

static void t_bcast_muldet_match(void)
{
	const data_id fid = data_open(FILENAME_MD);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open for muldet bcast test");

	struct circ_muldet md;
	TEST_EQ(data_muldet_load(fid, &md), 0);
	TEST_EQ(md.len, (size_t)N_DETS);

	/* Send rank-0's packed dets to all ranks via a second
	 * bcast, compare against what each rank loaded. */
	uint64_t ref_idx[N_DETS];
	_Complex double ref_cf[N_DETS];
	for (size_t i = 0; i < N_DETS; i++) {
		ref_idx[i] = md.dets[i].idx;
		ref_cf[i] = md.dets[i].cf;
	}
	MPI_Bcast(ref_idx, N_DETS, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(ref_cf, N_DETS * (int)sizeof(_Complex double),
		MPI_BYTE, 0, MPI_COMM_WORLD);
	for (size_t i = 0; i < N_DETS; i++) {
		TEST_ASSERT(ref_idx[i] == md.dets[i].idx,
			"det[%zu].idx differs between rank 0 and this rank",
			i);
		TEST_ASSERT(ref_cf[i] == md.dets[i].cf,
			"det[%zu].cf differs between rank 0 and this rank",
			i);
	}

	circ_muldet_free(&md);
	data_close(fid);
}

static void t_bcast_coeffmat_match(void)
{
	const data_id fid = data_open(FILENAME_CM);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open for coeffmat bcast test");

	struct data_coeff_matrix cm;
	TEST_EQ(data_coeff_matrix_load(fid, &cm), 0);
	TEST_EQ(cm.nqb, (uint32_t)CM_N_QB);
	TEST_EQ(cm.n_sites, (uint32_t)CM_N_SITES);
	TEST_EQ(cm.n_alpha, (uint32_t)CM_N_ALPHA);
	TEST_EQ(cm.n_beta, (uint32_t)CM_N_BETA);
	TEST_EQ(cm.closed_shell, 1);
	TEST_EQ(cm.tapered, 0);
	TEST_EQ(cm.n_components, (size_t)0);
	TEST_ASSERT(cm.C_alpha != NULL, "C_alpha must be allocated");

	/* Send rank-0's scalars + C_alpha to all ranks via a
	 * second bcast, compare against what each rank loaded. */
	const size_t sz_a = (size_t)cm.n_sites * cm.n_alpha;
	uint32_t ref_scl[4] = {
		cm.nqb, cm.n_sites, cm.n_alpha, cm.n_beta,
	};
	int ref_flags[2] = { cm.closed_shell, cm.tapered };
	double ref_ca[CM_N_SITES * CM_N_ALPHA];
	for (size_t i = 0; i < sz_a; i++)
		ref_ca[i] = cm.C_alpha[i];
	MPI_Bcast(ref_scl, 4, MPI_UINT32_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(ref_flags, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ref_ca, (int)sz_a, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	TEST_ASSERT(ref_scl[0] == cm.nqb, "nqb differs");
	TEST_ASSERT(ref_scl[1] == cm.n_sites, "n_sites differs");
	TEST_ASSERT(ref_scl[2] == cm.n_alpha, "n_alpha differs");
	TEST_ASSERT(ref_scl[3] == cm.n_beta, "n_beta differs");
	TEST_ASSERT(ref_flags[0] == cm.closed_shell, "closed_shell differs");
	TEST_ASSERT(ref_flags[1] == cm.tapered, "tapered differs");
	for (size_t i = 0; i < sz_a; i++)
		TEST_ASSERT(ref_ca[i] == cm.C_alpha[i],
			"C_alpha[%zu] differs between rank 0 and this rank",
			i);

	data_coeff_matrix_free(&cm);
	data_close(fid);
}

static void t_nan_padding_and_partial_write(int rank)
{
	const data_id fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID, "open for write test");

	struct data_circ_writer wr;
	TEST_EQ(data_circ_writer_init(fid, "circ_trott", N_STEPS, &wr), 0);
	for (size_t i = 0; i < N_WRITTEN; i++) {
		const _Complex double z = CMPLX((double)(i + 1), -(double)i);
		TEST_EQ(data_circ_write_step(&wr, i, z), 0);
	}
	data_circ_writer_close(&wr);
	data_close(fid);

	/* Rank 0 re-reads the dataset and asserts the written
	 * rows match while the trailing rows are NaN. */
	if (rank == 0) {
		const hid_t f = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
		const hid_t g = H5Gopen(f, "circ_trott", H5P_DEFAULT);
		const hid_t d = H5Dopen2(g, "values", H5P_DEFAULT);
		double buf[N_STEPS * 2];
		H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			buf);
		for (size_t i = 0; i < N_WRITTEN; i++) {
			TEST_ASSERT(buf[2 * i] == (double)(i + 1),
				"row %zu real: got %g, expected %g",
				i, buf[2 * i], (double)(i + 1));
			TEST_ASSERT(buf[2 * i + 1] == -(double)i,
				"row %zu imag: got %g, expected %g",
				i, buf[2 * i + 1], -(double)i);
		}
		for (size_t i = N_WRITTEN; i < N_STEPS; i++) {
			TEST_ASSERT(isnan(buf[2 * i]),
				"row %zu real: expected NaN, got %g",
				i, buf[2 * i]);
			TEST_ASSERT(isnan(buf[2 * i + 1]),
				"row %zu imag: expected NaN, got %g",
				i, buf[2 * i + 1]);
		}
		H5Dclose(d);
		H5Gclose(g);
		H5Fclose(f);
	}
}

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	if (wd.rank == 0) {
		TEST_EQ(build_hamil_fixture(), 0);
		TEST_EQ(build_muldet_fixture(), 0);
		TEST_EQ(build_coeffmat_fixture(), 0);
	}
	/* All ranks must wait until the files exist. */
	MPI_Barrier(MPI_COMM_WORLD);

	t_open_follower_sentinel(wd.rank);
	t_bcast_buffers_match();
	t_bcast_muldet_match();
	t_bcast_coeffmat_match();
	t_nan_padding_and_partial_write(wd.rank);

	if (wd.rank == 0) {
		remove(FILENAME);
		remove(FILENAME_MD);
		remove(FILENAME_CM);
	}
	world_free();
	return 0;
}
