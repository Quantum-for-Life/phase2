/*
 * t-data_idempotence -- cover the idempotent paths in the
 * data layer:
 *   - data_grp_create called twice on the same path succeeds.
 *   - data_close on an already-closed fid is a no-op (no SEGV).
 *   - data_open after data_close reopens the file cleanly.
 *
 * Single-rank by design; the MPI-side scenarios live in
 * t-data_mpi.
 */
#include "c23_compat.h"
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "ph2run/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x6c8b3eed7f29ad11)

static char *FILENAME = "/tmp/t-data_idempotence.h5";
static const char *GRP_NAME = "test_grp";

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	/* Single-rank only; followers participate via the collective
	 * data_* calls.  Rank 0 sets up an empty file. */
	if (wd.rank == 0) {
		remove(FILENAME);
		const hid_t f = H5Fcreate(FILENAME, H5F_ACC_EXCL,
			H5P_DEFAULT, H5P_DEFAULT);
		TEST_ASSERT(f != H5I_INVALID_HID, "H5Fcreate");
		H5Fclose(f);
	}

	/* grp_create twice -- second call succeeds. */
	data_id fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID,
		"data_open first: got %lld", (long long)fid);
	TEST_EQ(data_grp_create(fid, GRP_NAME), 0);
	TEST_EQ(data_grp_create(fid, GRP_NAME), 0);  /* idempotent */
	data_close(fid);

	/* data_close on already-closed sentinel (DATA_INVALID_FID)
	 * is a no-op.  Should not crash. */
	data_close(DATA_INVALID_FID);
	data_close(DATA_FOLLOWER_FID);

	/* data_open after data_close reopens cleanly. */
	fid = data_open(FILENAME);
	TEST_ASSERT(fid != DATA_INVALID_FID,
		"data_open second: got %lld", (long long)fid);
	data_close(fid);

	if (wd.rank == 0)
		remove(FILENAME);
	world_free();
	return 0;
}
