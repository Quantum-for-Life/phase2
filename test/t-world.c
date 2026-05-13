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

	if (wd.rank == 0)
		log_info("This is rank no. %d", wd.rank);

	world_free();

	world_info(&wd);
	TEST_ASSERT(wd.stat == WORLD_DONE, "wrong status");
}
