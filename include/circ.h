#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "data.h"
#include "qreg.h"

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

void circ_hamil_init(struct circ_hamil *h);
void circ_hamil_destroy(struct circ_hamil *h);
int  circ_hamil_from_file(struct circ_hamil *h, data_id fid);

void circ_multidet_init(struct circ_multidet *md);
void circ_multidet_destroy(struct circ_multidet *md);
int  circ_multidet_from_file(struct circ_multidet *md, data_id fid);

/* Circuit: trott */
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

/* Circuit: qdrift */
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
int  circ_qdirft_data_from_file(struct circ_qdrift_data *cd, data_id fid);
int  circ_qdrift_simulate(const struct circ_qdrift_data *cd);

#endif // CIRC_H
