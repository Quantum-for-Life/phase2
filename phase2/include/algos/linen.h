#ifndef PHASE2_LINEN_H
#define PHASE2_LINEN_H

#include "circ.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB (4)
#define LINEN_DEFAULT_NUM_SYS_QB (8)
#define LINEN_DEFAULT_NUM_ANC_QB (4)

struct linen_circuit_data {
	int val_prepst;
	int val_effect;
	int val_measure;
};

struct linen_circ_data {
	int pass_prepst;
	int pass_effect;
	int pass_measure;
};

void linen_circuit_init(struct circuit *ct, struct linen_circuit_data *ct_dat);

int linen_simulate(void);

#endif //PHASE2_LINEN_H
