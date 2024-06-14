#ifndef PHASE2_CIRC_QDRIFT_H
#define PHASE2_CIRC_QDRIFT_H

#include "circ.h"
#include "data.h"

struct circ_qdrift_data {
	struct circ_hamil    hamil;
	struct circ_multidet multidet;

	double step_size;
	size_t depth;

	double *samples[2];
	size_t	num_samples;
};

int  circ_qdrift_data_init(struct circ_qdrift_data *cd, data_id fid);
void circ_qdrift_data_destroy(struct circ_qdrift_data *cd);
int  circ_data_from_file(struct circ_qdrift_data *cd, data_id fid);
int  circ_qdrift_simulate(const struct circ_qdrift_data *cd);

#endif // PHASE2_CIRC_QDRIFT_H
