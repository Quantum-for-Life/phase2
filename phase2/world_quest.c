#include "c23_compat.h"
#include <stdlib.h>

#include "QuEST.h"

#include "phase2/world.h"

#include "world.h"
#include "world_quest.h"

int world_backend_init(struct world *wd)
{
	struct world_quest *w = malloc(sizeof *w);
	if (!w)
		return -1;

	w->env = createQuESTEnv();
	seedQuEST(&w->env, &wd->seed, 1);

	wd->data = w;

	return 0;
}

void world_backend_destroy(struct world *wd)
{
	struct world_quest *w = wd->data;
	if (w) {
		destroyQuESTEnv(w->env);
		free(w);
	}
}
