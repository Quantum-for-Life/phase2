#ifndef QDRIFT_H
#define QDRIFT_H

#include <stddef.h>

#include "phase2/circ.h"
#include "phase2/data.h"
#include "phase2/prob.h"
#include "xoshiro256ss.h"

struct qdrift_data {
	size_t depth;
	size_t samples;
	double step_size;
};

struct qdrift_ranct {
	struct circ_hamil hm_ran; /* randomized */
	struct prob_cdf cdf;
};

struct qdrift_smpl {
	_Complex double *z;
	size_t len;
};

struct qdrift {
	struct circ ct;
	struct qdrift_data dt;
	struct qdrift_ranct ranct;
	struct qdrift_smpl smpl;
	struct xoshiro256ss rng;
};

int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt, data_id fid);

void qdrift_free(struct qdrift *qd);

int qdrift_write_res(struct qdrift *qd, data_id fid);

#endif // QDRIFT_H
