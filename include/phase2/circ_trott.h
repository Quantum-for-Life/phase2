#ifndef CIRC_TROTT_H
#define CIRC_TROTT_H

#ifdef __cplusplus
extern "C" {
#endif

struct circ_trott_data {
	double delta;
	size_t nsteps;
};

struct circ_trott {
	struct circ circ;

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
