#ifndef QDRIFT_H
#define QDRIFT_H

#include <stddef.h>

#include "phase2/circ.h"
#include "phase2/data.h"
#include "xoshiro256ss.h"

struct qdrift_data {
	size_t depth;
	double step_size;
	size_t samples;
};

struct qdrift_rct {
	struct circ_hamil rhm;
};

struct qdrift_samples {
	_Complex double *z;
	size_t len;
};

struct qdrift {
	struct circ ct;
	struct qdrift_data dt;
	struct qdrift_rct rct;
	struct qdrift_samples smp;
	struct xoshiro256ss rng;
};

int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt, data_id fid);

void qdrift_destroy(struct qdrift *qd);

int qdrift_write_res(struct qdrift *qd, data_id fid);

#endif // QDRIFT_H
