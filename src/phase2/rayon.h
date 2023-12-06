#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#include <stdint.h>

#include "data.h"

#define RAYON_NAME "rayon"
#define RAYON_NUM_MEA_QB (1)
#define RAYON_NUM_ANC_QB (0)

typedef unsigned long pauli_pak_t;

const int a = sizeof(unsigned);

struct rayon_data_hamil {
	size_t num_qubits;
	size_t num_terms;
	double *coeffs;
	pauli_pak_t *pak;
};

struct rayon_data_multidet {
	size_t num_dets;
	struct {
		long long index;
		double coeff_re;
		double coeff_im;
	} *dets;
};

struct rayon_data_times {
	size_t num_steps;
	struct {
		double t;
		double val_re;
		double val_im;
	} *steps;
};

struct rayon_data {
	struct rayon_data_hamil hamil;
	struct rayon_data_multidet multidet;
	struct rayon_data_times times;
};

void rayon_data_init(struct rayon_data *rd);

void rayon_data_destroy(struct rayon_data *rd);

int rayon_data_from_data(struct rayon_data *rd, const struct data *dat);

void rayon_data_write_times(struct data_time_series *dat,
			    struct rayon_data_times *rt);

int rayon_simulate(struct circ_env *env, struct rayon_data *rd);

#endif //PHASE2_RAYON_H
