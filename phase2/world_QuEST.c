#include <stdlib.h>

#include "QuEST.h"

#include "phase2/world.h"
#include "world_QuEST.h"

int world_QuEST_init(struct world *wd) {

	struct world_QuEST *qe = malloc(sizeof *qe);
	if (!qe)
		return -1;

	qe->env = createQuESTEnv();
	seedQuEST(&qe->env, &wd->seed, 1);

	wd->data = qe;

	return 0;
}

int world_QuEST_destroy(struct world *wd)
{
	struct world_QuEST *qe = wd->data;
	if (!qe)
		return -1;

	destroyQuESTEnv(qe->env);
	free(qe);

	return 0;
}
