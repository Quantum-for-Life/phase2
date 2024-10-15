#include <stdlib.h>

#include "mpi.h"

#include "phase2/world.h"

int log_init(void);

#if PHASE2_BACKEND == 1 /* QuEST */
int world_QuEST_init(struct world *wd);
int world_QuEST_destroy(struct world *wd);
#else
static int world_QuEST_init(struct world *wd)
{
	(void)wd;
	return 0;
}
static int world_QuEST_destroy(struct world *wd)
{
	(void)wd;
	return 0;
}
#endif /* PHASE2_BACKEND == 1 */

#if PHASE2_BACKEND == 2 /* cuQuantum */
int world_cuQuantum_init(struct world *wd);
int world_cuQuantum_destroy(struct world *wd);
#else
static int world_cuQuantum_init(struct world *wd)
{
	(void)wd;
	return 0;
}
static int world_cuQuantum_destroy(struct world *wd)
{
	(void)wd;
	return 0;
}
#endif /* PHASE2_BACKEND == 1 */

static struct world WORLD = {
	.stat = WORLD_UNDEF,
	.size = 0,
	.rank = 0,
	.seed = UINT64_C(1111),
	.data = (void *)0,
};

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
	if (sz == 0)
		goto err;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rk) != MPI_SUCCESS)
		goto err;

	WORLD.size = sz;
	WORLD.rank = rk;

	WORLD.seed = seed;
	xoshiro256ss_init(&WORLD.rng, WORLD.seed);
	/* Split the state for parrallel distributed computation. */
	for (int i = 0; i < WORLD.rank; i++)
		xoshiro256ss_longjump(&WORLD.rng);

	if (world_QuEST_init(&WORLD) < 0)
		goto err;
	if (world_cuQuantum_init(&WORLD) < 0)
		goto err;

	return WORLD.stat = WORLD_READY;

err:
	return WORLD.stat = WORLD_ERR;
}

int world_destroy(void)
{
	if (WORLD.stat == WORLD_READY) {
		if (MPI_Finalize() == MPI_SUCCESS)
			WORLD.stat = WORLD_DONE;
		else
			WORLD.stat = WORLD_ERR;
	}

	if (world_QuEST_destroy(&WORLD) < 0)
		WORLD.stat = WORLD_ERR;
	if (world_cuQuantum_destroy(&WORLD) < 0)
		WORLD.stat = WORLD_ERR;

	return WORLD.stat;
}

int world_info(struct world *wd)
{
	wd->size = WORLD.size;
	wd->rank = WORLD.rank;
	wd->seed = WORLD.seed;
	wd->rng = WORLD.rng;
	wd->data = WORLD.data;

	return wd->stat = WORLD.stat;
}
