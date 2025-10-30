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

#ifdef __cplusplus
}
#endif

#endif /* WORLD_CUDA_H */
