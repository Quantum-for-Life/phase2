#ifndef PHASE2_LINEN_H
#define PHASE2_LINEN_H

#include "circ.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB 4
#define LINEN_DEFAULT_NUM_SYS_QB 8
#define LINEN_DEFAULT_NUM_ANC_QB 4

struct circuit linen_circuit_factory(void *data);

#endif //PHASE2_LINEN_H
