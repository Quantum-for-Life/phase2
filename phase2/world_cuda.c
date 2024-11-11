#include "c23_compat.h"
#include <stdlib.h>

#include <cuda_runtime_api.h>

#include "mpi.h"

#include "phase2/world.h"

#include "world_cuda.h"

int world_backend_init(struct world *wd)
{
	struct world_cuda *cu = malloc(sizeof *cu);
	if (!cu)
		return -1;

	/* Determine the local MPI rank (within the node). */
	int loc_rank, loc_size;
	MPI_Comm loc_comm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, wd->rank,
		MPI_INFO_NULL, &loc_comm);
	MPI_Comm_size(loc_comm, &loc_size);
	MPI_Comm_rank(loc_comm, &loc_rank);

	/* Assume one GPU per process. */
	cudaSetDevice(loc_rank % loc_size);

	cu->loc_rank = loc_rank;
	cu->loc_size = loc_size;
	wd->data = cu;

	return 0;
}

void world_backend_destroy(struct world *wd)
{
	struct world_cuda *cu = wd->data;
	free(cu);
}
