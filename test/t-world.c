#define LOG_SUBSYS "test"

#include "c23_compat.h"
#include <stdint.h>

#include "log.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x84a06b714640f7dc)

int main(void)
{
	struct world_info wd;

	/* Before init, status must be WORLD_UNDEF. */
	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_UNDEF, "wrong status");

	/* world_free on a never-initialised world: should not
	 * crash, MPI_Finalize is skipped, status stays
	 * WORLD_UNDEF. */
	TEST_ASSERT(world_free() == WORLD_UNDEF,
		"free without init must return WORLD_UNDEF");
	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_UNDEF, "status drifted on bare free");

	/* Zero seed must be rejected. */
	TEST_ASSERT(world_init(nullptr, nullptr, 0) == WORLD_ERR,
		"zero seed should fail");
	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_ERR, "wrong status after zero seed");

	/* Valid initialisation. */
	TEST_ASSERT(world_init(nullptr, nullptr, WD_SEED) == WORLD_READY,
		"error initializing the world");

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_READY, "wrong status");
	TEST_ASSERT(wd.seed == WD_SEED, "wrong seed");

	const int first_rank = wd.rank;
	const int first_size = wd.size;

	/* Re-init without an intervening free: MPI is already
	 * up, so MPI_Init is skipped; world_init repopulates
	 * the snapshot and the second call also returns
	 * WORLD_READY with the same rank / size. */
	TEST_ASSERT(world_init(nullptr, nullptr, WD_SEED) == WORLD_READY,
		"re-init must succeed");
	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_READY, "re-init dropped status");
	TEST_ASSERT(wd.rank == first_rank, "re-init changed rank");
	TEST_ASSERT(wd.size == first_size, "re-init changed size");

	if (wd.rank == 0)
		log_info("This is rank no. %d", wd.rank);

	world_free();

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_DONE, "wrong status");
}
