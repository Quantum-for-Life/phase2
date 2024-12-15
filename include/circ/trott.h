#ifndef TROTT_H
#define TROTT_H

#include "phase2/circ.h"
#include "phase2/data.h"

struct trott_data {
	double delta;
	size_t steps;
};

struct trott_steps {
	_Complex double *z;
	size_t len;
};

struct trott {
	struct circ ct;
	struct trott_data dt;
};

int trott_init(struct trott *tt, const struct trott_data *dt, data_id fid);

void trott_free(struct trott *tt);

int trott_simul(struct trott *tt);

int trott_write_res(struct trott *tt, data_id fid);

#endif // TROTT_H
