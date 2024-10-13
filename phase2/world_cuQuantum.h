#ifndef WORLD_CUQUANTUM_H
#define WORLD_CUQUANTUM_H

#ifdef __cplusplus
extern "C" {
#endif

#include "custatevec.h"

struct world_cuQuantum {
	int local_size, local_rank;
	custatevecHandle_t handle;
};

#ifdef __cplusplus
}
#endif

#endif /* WORLD_CUQUANTUM_H */
