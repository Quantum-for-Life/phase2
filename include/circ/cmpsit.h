#ifndef CMPSIT_H
#define CMPSIT_H

#include <stddef.h>

#include "phase2.h"
#include "xoshiro256ss.h"

#define CMPSIT_TRUNC_DIST (99UL)

struct cmpsit_data {
	size_t depth;
	size_t length;
	size_t samples;
	double step_size;
	size_t steps;
};

/* Sampled circuit. */
struct cmpsit_ranct {
	struct circ_hamil hm_ran;
	struct prob_cdf cdf;
};

/* Results */
struct cmpsit_smpl {
	_Complex double *z;
	size_t len;
};

struct cmpsit {
	struct circ ct;
	struct cmpsit_data dt;
	struct cmpsit_ranct ranct;
	struct cmpsit_smpl smpl;
	struct xoshiro256ss rng;
};

int cmpsit_init(struct cmpsit *cp, const struct cmpsit_data *dt, data_id fid);

void cmpsit_free(struct cmpsit *cp);

int cmpsit_simul(struct cmpsit *cp);

int cmpsit_write_res(struct cmpsit *cp, data_id fid);

#endif // CMPSIT_H
