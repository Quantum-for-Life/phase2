#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#include "data.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
//#define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0

struct rayon_data_hamil {
	size_t num_qubits;
	size_t num_terms;
	double *coeffs;
	int *paulis;
};

struct rayon_data_multidet {
	size_t num_dets;
	struct {
		long long det;
		double coeff_real;
		double coeff_imag;
	} *dets;
};

struct rayon_data {
	struct rayon_data_hamil hamil;
	struct rayon_data_multidet multidet;
};

void rayon_data_init(struct rayon_data *ct_dat);

void rayon_data_destroy(struct rayon_data *ct_dat);

int rayon_data_from_data(struct rayon_data *ct_dat, const struct data *dat);

int rayon_simulate(struct circ_env env, const struct rayon_data *ct_dat,
		   const struct data_time_series *dat_ts);

#endif //PHASE2_RAYON_H
