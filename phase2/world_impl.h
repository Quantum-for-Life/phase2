#ifndef WORLD_IMPL_H
#define WORLD_IMPL_H

#include "phase2/world.h"

#ifdef __cplusplus
extern "C" {
#endif

int world_backend_init(struct world *wd);

void world_backend_destroy(struct world *wd);

#ifdef __cplusplus
}
#endif

#endif /* WORLD_IMPL_H */
