#include <stdlib.h>

#include <cuda_runtime_api.h>

#include "mpi.h"

#include "phase2/world.h"
#include "world_cuda.h"

int world_cuQuantum_init(struct world *wd)
{
	struct world_cuQuantum *cu = malloc(sizeof *cu);
	if (!cu)
		return -1;

	/* Determine the local MPI rank (within the node). */
	int local_rank, local_size;
	MPI_Comm local_comm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,
			wd->rank, MPI_INFO_NULL, &local_comm);
	MPI_Comm_size(local_comm, &local_size);
	MPI_Comm_rank(local_comm, &local_rank);

	/* Assume one GPU per process. */
	cudaSetDevice(local_rank % local_size);

 	cu->local_rank = local_rank;
	cu->local_size = local_size;
	wd->data = cu;

	return 0;
}

int world_cuQuantum_destroy(struct world *wd)
{
	struct world_cuQuantum *cu = wd->data;
	if (!cu)
		return -1;

	free(cu);

	return 0;
}
