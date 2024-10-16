#ifndef WORLD_QUEST_H
#define WORLD_QUEST_H

#include "c23_compat.h"

#include "QuEST.h"

#ifdef __cplusplus
extern "C" {
#endif

struct world_quest {
	QuESTEnv env;
};

int world_quest_init(struct world *wd);
int world_quest_destroy(struct world *wd);

#ifdef __cplusplus
}
#endif

#endif /* WORLD_QUEST_H */
