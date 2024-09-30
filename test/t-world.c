#include "phase2/world.h"

#include "test.h"

void TEST_MAIN(void)
{
	struct world wd;

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_UNDEF, "wrong status");

	TEST_ASSERT(world_init((void *)0, (void *)0) == WORLD_READY,
		"error initializing the world");

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_READY, "wrong status");

	if (wd.rank == 0)
		log_info("This is rank no. %d", wd.rank);

	world_fin();

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_DONE, "wrong status");
}
