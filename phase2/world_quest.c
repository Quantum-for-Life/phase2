#include "c23_compat.h"
#include <stdlib.h>

#include "QuEST.h"

#include "phase2/world.h"

#include "world_impl.h"
#include "world_quest.h"

int world_backend_init(struct world *wd)
{
	struct world_quest *w = malloc(sizeof *w);
	if (w == nullptr)
		return -1;

	w->env = createQuESTEnv();
	seedQuEST(&w->env, &wd->seed, 1);

	wd->data = w;

	return 0;
}

void world_backend_destroy(struct world *wd)
{
	struct world_quest *w = wd->data;

	destroyQuESTEnv(w->env);
	free(w);
}
