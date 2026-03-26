#ifndef CMPSIT_H
#define CMPSIT_H

#include <stddef.h>

#include "phase2.h"
#include "xoshiro256ss.h"

struct cmpsit_data {
	uint64_t seed;
	size_t length;
	size_t depth;
	size_t steps;
	double angle_det;
	double angle_rand;
	size_t samples;
};

/* Sampled circuit. */
struct cmpsit_ranct {
	struct circ_hamil hm_det, hm_ran, hm_smpl;
	struct prob_cdf cdf;
	double lambda_r;
	size_t depth;
	double angle_rand;
};

struct cmpsit {
	struct circ ct;
	struct cmpsit_data dt;
	struct cmpsit_ranct ranct;
	struct xoshiro256ss rng;
};

int cmpsit_init(struct cmpsit *cp, const struct cmpsit_data *dt, data_id fid);

void cmpsit_free(struct cmpsit *cp);

int cmpsit_simul(struct cmpsit *cp);

int cmpsit_write_res(struct cmpsit *cp, data_id fid);

#endif // CMPSIT_H
