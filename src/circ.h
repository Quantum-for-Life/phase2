#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdint.h>

#include "data.h"

struct circ_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	struct paulis *paulis;
};

struct circ_multidet {
	struct {
		uint64_t idx;
		double	 coeff[2];
	}     *dets;
	size_t num_dets;
};

struct circ_data {
	struct circ_hamil    hamil;
	struct circ_multidet multidet;

	double time_factor;

	double *samples[2];
	size_t	num_samples;

	double step_size;
	size_t depth;
};

int
circ_data_init(struct circ_data *cd, data_id fid);

void
circ_data_destroy(struct circ_data *cd);

int
circ_simulate(const struct circ_data *cd);

#endif // PHASE2_CIRC_H
