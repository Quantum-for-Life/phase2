#include <stdlib.h>

#include "custatevec.h"

#include "phase2/world.h"
#include "world_cuQuantum.h"

int world_cuQuantum_init(struct world *wd)
{
	struct world_cuQuantum *cu = malloc(sizeof *cu);
	if (!cu)
		return -1;

	if (custatevecCreate(&cu->handle) != CUSTATEVEC_STATUS_SUCCESS)
		goto err;
	wd->data = cu;

	return 0;
	
err:
	free(cu);
	return -1;
}

int world_cuQuantum_destroy(struct world *wd)
{
	struct world_cuQuantum *cu = wd->data;
	if (!cu)
		return -1;

	if (custatevecDestroy(cu->handle) != CUSTATEVEC_STATUS_SUCCESS)
		return -1;

	free(cu);

	return 0;
}
