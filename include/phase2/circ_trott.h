#ifndef CIRC_TROTT_H
#define CIRC_TROTT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "phase2/qreg.h"

struct circ_trott_data {
	double delta;
	size_t nsteps;
};

struct circ_trott {
	struct qreg reg;
	struct circ circ;
	struct circ_cache cache;

	double delta;

	struct {
		_Complex double *steps;
		size_t nsteps;
	} res;
};

int circ_trott_init(
	struct circ_trott *tt, struct circ_trott_data *data, data_id fid);

void circ_trott_destroy(struct circ_trott *tt);

#ifdef __cplusplus
}
#endif

#endif // CIRC_TROTT_H
