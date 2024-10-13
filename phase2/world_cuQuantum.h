#include "custatevec.h"

struct world_cuQuantum {
	int local_size, local_rank;
	custatevecHandle_t handle;
};
