/*
 * Test routine: data_trott_write() from data.h by creating a fresh H5 file
 * Check if the values are written correctly.  Remove temporary file.
 */
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x3a1e27cce387af66)

#define SIZE (5)
#define MARGIN (1.0e-6)

_Complex double tst_vals[SIZE] = {
	CMPLX(0.0, 000.01),
	CMPLX(1.8, 000.11),
	CMPLX(2.7, 001.11),
	CMPLX(3.6, 011.11),
	CMPLX(4.5, 111.11),
};

static double delta = 0.3224;

static char *FILENAME = "/tmp/G1w1Clar2ZLovBir2cGYUbCxgIaV4";

static int prepare_test_file(hid_t file_id)
{
	int rt = -1;

	/* Create main group */
	const hid_t grp_id = H5Gcreate(
		file_id, DATA_CIRCTROTT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID)
		goto ex_create;

	// create a scalar (singleton) attribute
	if (data_attr_write(
		    file_id, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA, delta) < 0)
		goto ex_attr;

	rt = 0;

ex_attr:
	H5Gclose(grp_id);
ex_create:
	return rt;
}

static void TEST_MAIN(void)
{
	struct world wd;
	world_init((void *)0, (void *)0, WD_SEED);
	world_info(&wd);

	if (wd.rank == 0)
		remove(FILENAME);

	/*
	 * Set up file access property list with parallel I/O access
	 */
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

	hid_t file_id =
		H5Fcreate(FILENAME, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("create H5 file");
		goto ex_create;
	}

	if (prepare_test_file(file_id) < 0) {
		TEST_FAIL("prepare test file");
		goto ex_prepare;
	};
	H5Fclose(file_id);

	data_id fid = data_open(FILENAME);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("data: reopen file");
		goto ex_dat2_open;
	}

	data_res_write(
		fid, DATA_CIRCTROTT, DATA_CIRCTROTT_VALUES, tst_vals, SIZE);
	data_close(fid);

	file_id = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("open H5 file");
		goto ex_create;
	}

	double d;
	if (data_attr_read(file_id, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA, &d) <
		0)
		TEST_FAIL("read delta");
	if (fabs(d - delta) > MARGIN)
		TEST_FAIL("wrong value of delta: %f", d);

	hid_t grp_id = H5Gopen(file_id, DATA_CIRCTROTT, H5P_DEFAULT);
	hid_t dset = H5Dopen2(grp_id, DATA_CIRCTROTT_VALUES, H5P_DEFAULT);
	_Complex double val_read[SIZE];
	if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    val_read) < 0) {
		goto ex_readh5;
	}

	for (size_t i = 0; i < SIZE; i++) {
		if (cabs(tst_vals[i] - val_read[i]) > MARGIN) {
			TEST_FAIL("wrong value: %f+%fi (expected: %f+%fi)",
				creal(val_read[i]), cimag(val_read[i]),
				creal(tst_vals[i]), cimag(tst_vals[i]));
			goto ex_readh5;
		}
	}

ex_readh5:
	H5Dclose(dset);
	H5Gclose(grp_id);
	H5Fclose(file_id);
	fid = data_open(FILENAME);
	data_close(fid);
ex_dat2_open:
ex_prepare:
	/* Delete temporary file */
	if (wd.rank == 0 && remove(FILENAME) != 0)
		TEST_FAIL("remove temp file");
ex_create:
	world_free();
}
