#ifndef CIRC_TROTT_H
#define CIRC_TROTT_H

#include "phase2/circ.h"
#include "phase2/data.h"

struct trott_data {
	double delta;
	size_t nsteps;
};

struct trott_steps {
	_Complex double *z;
	size_t len;
};

struct trott {
	struct circ circ;
	double delta;
	struct trott_steps steps;
};

int trott_init(struct trott *tt, struct trott_data *data, data_id fid);

void trott_destroy(struct trott *tt);

#endif // CIRC_TROTT_H
