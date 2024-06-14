#ifndef PHASE2_CIRC_TROTT_H
#define PHASE2_CIRC_TROTT_H

#include "circ.h"
#include "data.h"

struct circ_trott_data {
	struct circ_hamil    hamil;
	struct circ_multidet multidet;

	double time_factor;

	double *trott_steps[2];
	size_t	num_trott_steps;
};

int  circ_trott_data_init(struct circ_trott_data *cd, size_t num_steps);
void circ_trott_data_destroy(struct circ_trott_data *cd);
int  circ_trott_data_from_file(struct circ_trott_data *cd, data_id fid);
int  circ_trott_simulate(const struct circ_trott_data *cd);

#endif // PHASE2_CIRC_TROTT_H
