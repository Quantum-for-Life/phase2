#include <mpi.h>

#include "error.h"

#include "ev.h"

int ev_init(struct ev *ev)
{
	int initialized, num_ranks, rank;

	MPI_Initialized(&initialized);
	if (!initialized && MPI_Init(NULL, NULL) != MPI_SUCCESS)
		return -EMPI;

	if (MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) != MPI_SUCCESS)
		return -EMPI;
	if (num_ranks == 0)
		return -ESIZE;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
		return -EMPI;

	ev->num_ranks = num_ranks;
	ev->rank      = rank;

	return OK;
}

int ev_destroy(struct ev *ev)
{
	(void)ev;

	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized && MPI_Finalize() != MPI_SUCCESS)
		return -EMPI;

	ev->num_ranks = 0;

	return OK;
}
