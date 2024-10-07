#include <stdint.h>

#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x84a06b714640f7dc)

static void TEST_MAIN(void)
{
	struct world wd;

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_UNDEF, "wrong status");

	TEST_ASSERT(world_init((void *)0, (void *)0, WD_SEED) == WORLD_READY,
		"error initializing the world");

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_READY, "wrong status");
	TEST_ASSERT(wd.seed == WD_SEED, "wrong seed");

	if (wd.rank == 0)
		log_info("This is rank no. %d", wd.rank);

	world_destroy();

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_DONE, "wrong status");
}
