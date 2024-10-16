#ifndef WORLD_CUDA_H
#define WORLD_CUDA_H

#include "c23_compat.h"

#ifdef __cplusplus
extern "C" {
#endif

struct world_cuda {
	int loc_rank;
	int loc_size;
};

int world_cuda_init(struct world *wd);
int world_cuda_destroy(struct world *wd);

#ifdef __cplusplus
}
#endif

#endif /* WORLD_CUDA_H */
