/*
 * t-data_dets_validate -- verify that circ_muldet_load rejects
 * determinant bytes outside {0, 1} with a -1 return and a
 * log_error line (instead of silently corrupting the basis-
 * state index).
 *
 * Single-rank: builds a tiny multidet group in-process with
 * a deliberately malformed `dets` value, then runs the loader
 * and asserts the rejection plus a zeroed output struct.
 */
#include "c23_compat.h"
#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xa9c45d18f37c2090)

static char *FILENAME = "/tmp/t-data_dets_validate.h5";

static int build_bad_multidet(void)
{
	const hid_t f = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT,
		H5P_DEFAULT);
	if (f == H5I_INVALID_HID)
		return -1;
	const hid_t sp = H5Gcreate(f, "state_prep", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);
	const hid_t md = H5Gcreate(sp, "multidet", H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT);

	/* coeffs: shape (1, 2), real + imag. */
	const hsize_t cd[2] = { 1, 2 };
	const hid_t csp = H5Screate_simple(2, cd, NULL);
	const hid_t cds = H5Dcreate2(md, "coeffs", H5T_IEEE_F64LE, csp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const double cv[2] = { 1.0, 0.0 };
	H5Dwrite(cds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cv);
	H5Dclose(cds);
	H5Sclose(csp);

	/* dets: shape (1, 3), one row with a bad byte (2). */
	const hsize_t dd[2] = { 1, 3 };
	const hid_t dsp = H5Screate_simple(2, dd, NULL);
	const hid_t dds = H5Dcreate2(md, "dets", H5T_NATIVE_UCHAR, dsp,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const unsigned char dv[3] = { 1, 2, 0 };  /* the 2 is illegal */
	H5Dwrite(dds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, dv);
	H5Dclose(dds);
	H5Sclose(dsp);

	H5Gclose(md);
	H5Gclose(sp);
	H5Fclose(f);
	return 0;
}

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	if (wd.rank == 0)
		TEST_EQ(build_bad_multidet(), 0);

	const data_id fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID, "data_open");

	struct circ_muldet md = { 0 };
	const int rt = circ_muldet_load(fid, &md);
	TEST_ASSERT(rt < 0,
		"circ_muldet_load must reject dets[1]=2 (got rt=%d)", rt);
	TEST_ASSERT(md.dets == NULL && md.len == 0,
		"failed load must leave the struct without dangling buffers");

	data_close(fid);
	if (wd.rank == 0)
		remove(FILENAME);
	world_free();
	return 0;
}
