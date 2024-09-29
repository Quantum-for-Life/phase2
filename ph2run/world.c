#include "mpi.h"

#include "world.h"

int world_init(struct world *wd)
{
	int init, nrk, rk;

	MPI_Initialized(&init);
	if (!init && MPI_Init(NULL, NULL) != MPI_SUCCESS)
		return -1;

	if (MPI_Comm_size(MPI_COMM_WORLD, &nrk) != MPI_SUCCESS)
		return -1;
	if (nrk == 0)
		return -1;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rk) != MPI_SUCCESS)
		return -1;

	wd->num_ranks = nrk;
	wd->rank = rk;

	return 0;
}

