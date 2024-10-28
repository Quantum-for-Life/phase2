#ifndef CIRC_TROTT_H
#define CIRC_TROTT_H

#ifdef __cplusplus
extern "C" {
#endif

struct circ_trott_data {
	double delta;
	size_t nsteps;
};

struct circ_trott_res {
	double delta;
	_Complex double *steps;
	size_t nsteps;
};

#ifdef __cplusplus
}
#endif

#endif // CIRC_TROTT_H
