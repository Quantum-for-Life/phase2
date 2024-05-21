#ifndef PHASE2_CIRC_QDRIFT_H
#define PHASE2_CIRC_QDRIFT_H

#include <stdint.h>

#include "data.h"

struct circ_qdrift_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	struct paulis *paulis;
};

struct circ_qdrift_multidet {
	struct {
		uint64_t idx;
		double	 coeff[2];
	}     *dets;
	size_t num_dets;
};

struct circ_qdrift_data {
	struct circ_qdrift_hamil hamil;
	struct circ_qdrift_multidet multidet;

	double *samples[2];
	size_t	num_samples;
	double	step_size;
	size_t	depth;
};

int
circ_qdrift_data_init(struct circ_qdrift_data *cd, data_id fid);

void
circ_qdrift_data_destroy(struct circ_qdrift_data *cd);

int
circ_qdrift_simulate(const struct circ_qdrift_data *cd);

#endif // PHASE2_CIRC_QDRIFT_H
