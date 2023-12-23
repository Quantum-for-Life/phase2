#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#include "circ.h"
#include "data.h"

#define RAYON_NAME "rayon"
#define RAYON_NUM_MEA_QB (1)
#define RAYON_NUM_ANC_QB (0)

struct rayon_data_multidet {
	size_t num_dets;
	struct {
		long long	index;
		_Complex double coeff;
	} *dets;
};

struct rayon_data_times {
	size_t num_steps;
	struct {
		double		t;
		_Complex double val;
	} *steps;
};

struct rayon_data {
	struct circ_hamil	   hamil;
	struct rayon_data_multidet multidet;
	struct rayon_data_times	   times;
};

void
rayon_data_init(struct rayon_data *rd);

void
rayon_data_destroy(struct rayon_data *rd);

int
rayon_data_from_data(
	struct rayon_data *rd, const struct data *dat, data_id fid);

void
rayon_data_write_times(
	struct data_time_series *dat, const struct rayon_data_times *rt);

int
rayon_simulate(const struct rayon_data *rd);

#endif // PHASE2_RAYON_H
