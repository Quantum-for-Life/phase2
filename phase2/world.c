#include "c23_compat.h"
#include <stdlib.h>

#include "mpi.h"

#include "log.h"
#include "xoshiro256ss.h"

#include "phase2/world.h"

static struct world_info WORLD = {
	.stat = WORLD_UNDEF,
	.size = 0,
	.rank = 0,
	.seed = UINT64_C(0x77dd8e60521fb661),
};

static struct xoshiro256ss rng;

int world_backend_init(const struct world_info *wd);
void world_backend_destroy(const struct world_info *wd);

int world_init(int *argc, char ***argv, uint64_t seed)
{
	int init, sz, rk;

	if (seed == 0)
		goto err;
	if (log_init() < 0)
		goto err;

	MPI_Initialized(&init);
	if (!init && MPI_Init(argc, argv) != MPI_SUCCESS)
		goto err;
	if (MPI_Comm_size(MPI_COMM_WORLD, &sz) != MPI_SUCCESS)
		goto err;
	if (sz == 0 || (sz & (sz - 1)) != 0) {
		log_error("Number of MPI processes (%u) must"
			  " be a power of two.",
			sz);
		goto err;
	}
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rk) != MPI_SUCCESS)
		goto err;

	WORLD.size = sz;
	WORLD.rank = rk;

	WORLD.seed = seed;
	xoshiro256ss_init(&rng, WORLD.seed);
	/* Split the state for parrallel distributed computation. */
	for (int i = 0; i < WORLD.rank; i++)
		xoshiro256ss_longjump(&rng);

	if (world_backend_init(&WORLD) < 0)
		goto err;

	WORLD.stat = WORLD_READY;
	log_info("*** Init ***");
	log_info("World size: %d", WORLD.size);
	log_info("Backend: %s", WORLD_BACKEND);

	return WORLD.stat;

err:
	return WORLD.stat = WORLD_ERR;
}

int world_free(void)
{
	if (WORLD.stat == WORLD_READY) {
		if (MPI_Finalize() == MPI_SUCCESS)
			WORLD.stat = WORLD_DONE;
		else
			WORLD.stat = WORLD_ERR;
	}

	world_backend_destroy(&WORLD);

	return WORLD.stat;
}

int world_info(struct world_info *wd)
{
	if (wd) {
		wd->stat = WORLD.stat;
		wd->size = WORLD.size;
		wd->rank = WORLD.rank;
		wd->seed = WORLD.seed;
	}

	return WORLD.stat;
}

#if PHASE2_BACKEND == 0 /* qreg */
inline int world_backend_init(const struct world_info *wd)
{
	(void)wd;

	return 0;
}

inline void world_backend_destroy(const struct world_info *wd)
{
	(void)wd;
}
#endif /* PHASE2_BACKEND == 0 */
