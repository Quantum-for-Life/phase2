#include "c23_compat.h"
#include <stdlib.h>

#include "QuEST.h"

#include "phase2/world.h"
#include "world_quest.h"

int world_quest_init(struct world *wd)
{
	struct world_quest *w = malloc(sizeof *w);
	if (w == nullptr)
		return -1;

	w->env = createQuESTEnv();
	seedQuEST(&w->env, &wd->seed, 1);

	wd->data = w;

	return 0;
}

int world_quest_destroy(struct world *wd)
{
	struct world_quest *w = wd->data;
	if (w == nullptr)
		return -1;

	destroyQuESTEnv(w->env);
	free(w);

	return 0;
}
