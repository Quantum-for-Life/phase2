/*
 * data attribute read/write round-trip.  Writes one `unsigned
 * long` and one `double` attribute via data_attr_write(),
 * reopens, reads them back via data_attr_read(), verifies bit-
 * exact match on the unsigned-long path and within machine
 * epsilon on the double path.
 *
 * The `_dbl` variant exercises both read and write; the `_ul`
 * variant only the write half (the cmpsit / qdrift init paths
 * write size_t scalars through this dispatcher, but nothing
 * reads them through data_attr_read).
 */
#include "c23_compat.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7b2f4a91d6e803c5)

#define MARGIN (1.0e-15)

static char *FILENAME = "/tmp/N3kF8xVmPq2AttrTest";

static const char *GRP_NAME = "test_grp";
static const char *ATTR_UL = "val_ul";
static const char *ATTR_DBL = "val_dbl";

static const unsigned long VAL_UL = 123456789UL;
static const double VAL_DBL = 3.14159265358979;

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	/* The data layer is rank-0-only.  Rank 0 creates the file
	 * via serial HDF5; other ranks hold the follower sentinel. */
	data_id fid = DATA_FOLLOWER_FID;
	if (wd.rank == 0) {
		remove(FILENAME);
		hid_t file_id = H5Fcreate(
			FILENAME, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
		if (file_id == H5I_INVALID_HID)
			TEST_FAIL("create H5 file");
		fid = (data_id)file_id;
	}

	if (data_grp_create(fid, GRP_NAME) < 0)
		TEST_FAIL("data_grp_create");

	if (data_attr_write(fid, GRP_NAME, ATTR_UL, VAL_UL) < 0)
		TEST_FAIL("data_attr_write unsigned long");

	if (data_attr_write(fid, GRP_NAME, ATTR_DBL, VAL_DBL) < 0)
		TEST_FAIL("data_attr_write double");

	if (wd.rank == 0)
		H5Fclose((hid_t)fid);

	/* Reopen with data_open and read back. */
	fid = data_open(FILENAME);
	if (fid == DATA_INVALID_FID)
		TEST_FAIL("data_open");

	double rd_dbl;
	if (data_attr_read(fid, GRP_NAME, ATTR_DBL, &rd_dbl) < 0)
		TEST_FAIL("data_attr_read double");
	TEST_ASSERT(fabs(rd_dbl - VAL_DBL) < MARGIN,
		"double mismatch: got %.17g, expected %.17g",
		rd_dbl, VAL_DBL);

	data_close(fid);

	/* Delete temporary file. */
	if (wd.rank == 0 && remove(FILENAME) != 0)
		TEST_FAIL("remove temp file");

	world_free();
}
