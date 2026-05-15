/*
 * t-data_hamil_validate -- exercise the error paths of
 * data_hamil_load against malformed /pauli_hamil fixtures
 * built in-process.  Hamil counterpart to
 * t-data_dets_validate.
 *
 * Each case rebuilds a tiny pauli_hamil group on rank 0
 * with one specific defect, opens the file via data_open,
 * and asserts that data_hamil_load returns -1 and leaves
 * the output struct without dangling buffers.
 */
#include "c23_compat.h"
#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "mpi.h"

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xab1d9e3f01c47a92)

#define N_TERMS (3)
#define N_QB (2)

static char *FILENAME = "/tmp/t-data_hamil_validate.h5";

static void write_norm(hid_t g)
{
	const hid_t asp = H5Screate(H5S_SCALAR);
	const hid_t aid = H5Acreate2(g, "normalization", H5T_IEEE_F64LE, asp,
		H5P_DEFAULT, H5P_DEFAULT);
	const double v = 1.0;
	H5Awrite(aid, H5T_NATIVE_DOUBLE, &v);
	H5Aclose(aid);
	H5Sclose(asp);
}

static void write_coeffs(hid_t g, hsize_t n)
{
	const hsize_t dims[1] = { n };
	const hid_t sp = H5Screate_simple(1, dims, NULL);
	const hid_t ds = H5Dcreate2(g, "coeffs", H5T_IEEE_F64LE, sp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	double *buf = malloc(sizeof *buf * n);
	for (hsize_t i = 0; i < n; i++)
		buf[i] = (double)(i + 1) * 0.1;
	H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	free(buf);
	H5Dclose(ds);
	H5Sclose(sp);
}

static void write_paulis(hid_t g, hsize_t nterms, hsize_t nqb,
	int bad_byte, int paulis_1d)
{
	const hsize_t dims2[2] = { nterms, nqb };
	const hsize_t dims1[1] = { nterms * nqb };
	const hid_t sp = paulis_1d
		? H5Screate_simple(1, dims1, NULL)
		: H5Screate_simple(2, dims2, NULL);
	const hid_t ds = H5Dcreate2(g, "paulis", H5T_NATIVE_UCHAR, sp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	unsigned char *buf = malloc(sizeof *buf * nterms * nqb);
	for (hsize_t i = 0; i < nterms * nqb; i++)
		buf[i] = (unsigned char)((i % 3) + 1);  /* X, Y, Z */
	if (bad_byte)
		buf[1] = 4;  /* paulis[0][1] = 4 */
	H5Dwrite(ds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	free(buf);
	H5Dclose(ds);
	H5Sclose(sp);
}

/*
 * Build /pauli_hamil with all pieces, then selectively
 * omit one (`omit_*`), inject an out-of-range Pauli byte
 * (`bad_pauli_byte`), or write paulis as a 1-D dataset
 * (`paulis_1d`) -- the dim-2 reader rejects anything else.
 */
static int build_fixture(int omit_norm, int omit_coeffs, int omit_paulis,
	int bad_pauli_byte, int paulis_1d)
{
	const hid_t f = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT,
		H5P_DEFAULT);
	if (f == H5I_INVALID_HID)
		return -1;
	const hid_t g = H5Gcreate(f, "pauli_hamil", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);

	if (!omit_norm)
		write_norm(g);
	if (!omit_coeffs)
		write_coeffs(g, (hsize_t)N_TERMS);
	if (!omit_paulis)
		write_paulis(g, (hsize_t)N_TERMS, (hsize_t)N_QB,
			bad_pauli_byte, paulis_1d);

	H5Gclose(g);
	H5Fclose(f);
	return 0;
}

static void run_case(const char *what, int omit_norm, int omit_coeffs,
	int omit_paulis, int bad_pauli_byte, int paulis_1d)
{
	struct world_info wd;
	world_info(&wd);

	if (wd.rank == 0) {
		remove(FILENAME);
		TEST_EQ(build_fixture(omit_norm, omit_coeffs, omit_paulis,
			bad_pauli_byte, paulis_1d), 0);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	const data_id fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID, "%s: data_open", what);

	struct circ_hamil hm = { 0 };
	const int rt = data_hamil_load(fid, &hm);
	TEST_ASSERT(rt < 0, "%s: expected rt<0, got %d", what, rt);
	TEST_ASSERT(hm.terms == NULL && hm.len == 0,
		"%s: failed load left dangling buffer", what);

	circ_hamil_free(&hm);
	data_close(fid);
}

static void t_missing_norm(void)
{
	run_case("missing normalization", 1, 0, 0, 0, 0);
}

static void t_missing_coeffs(void)
{
	run_case("missing coeffs", 0, 1, 0, 0, 0);
}

static void t_missing_paulis(void)
{
	run_case("missing paulis", 0, 0, 1, 0, 0);
}

static void t_bad_pauli_byte(void)
{
	run_case("bad pauli byte", 0, 0, 0, 1, 0);
}

static void t_paulis_wrong_rank(void)
{
	run_case("paulis wrong rank", 0, 0, 0, 0, 1);
}

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	t_missing_norm();
	t_missing_coeffs();
	t_missing_paulis();
	t_bad_pauli_byte();
	t_paulis_wrong_rank();

	if (wd.rank == 0)
		remove(FILENAME);
	world_free();
	return 0;
}
