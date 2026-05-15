/*
 * t-data_mpi -- MPI-specific behaviour of the data layer:
 *
 *   1. data_open returns a real fid on rank 0 and
 *      DATA_FOLLOWER_FID on other ranks.
 *   2. The collective Hamiltonian read via circ_hamil_load
 *      produces byte-equal packed terms on every rank
 *      (Bcast worked).
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
#include "phase2/data.h"
#include "phase2/paulis.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x21847b9ef0d3c2a4)

#define N_TERMS (5)
#define N_QB (3)
#define N_STEPS (8)
#define N_WRITTEN (4)

static char *FILENAME = "/tmp/t-data_mpi.h5";

/* Rank 0 builds a tiny pauli_hamil group; the test then
 * reads it through circ_hamil_load on every rank and
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
	TEST_EQ(circ_hamil_load(fid, &hm), 0);
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

	if (wd.rank == 0)
		TEST_EQ(build_hamil_fixture(), 0);
	/* All ranks must wait until the file exists. */
	MPI_Barrier(MPI_COMM_WORLD);

	t_open_follower_sentinel(wd.rank);
	t_bcast_buffers_match();
	t_nan_padding_and_partial_write(wd.rank);

	if (wd.rank == 0)
		remove(FILENAME);
	world_free();
	return 0;
}
