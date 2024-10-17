#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"

#ifdef __cplusplus
extern "C" {
#endif

struct circ_hamil {
	size_t nqb;
	size_t nterms;
	double *cfs; /* coefficients */
	struct paulis *ops;
};

int circ_hamil_init(struct circ_hamil *h, size_t nterm);

/* Convenience function. Use instead of circ_hamil_init */
int circ_hamil_init_from_file(struct circ_hamil *h, data_id fid);

void circ_hamil_destroy(struct circ_hamil *h);

struct circ_multidet {
	struct {
		uint64_t idx;
		double cf[2];
	} *dets;
	size_t ndets;
};

int circ_multidet_init(struct circ_multidet *md, size_t ndet);

/* Convenience function. Use instead of circ_multidet_init */
int circ_multidet_init_from_file(struct circ_multidet *md, data_id fid);

void circ_multidet_destroy(struct circ_multidet *md);

/* Circuit: trott */
struct circ_trott_data {
	struct circ_hamil hamil;
	struct circ_multidet multidet;

	double delta;

	double *tsteps[2]; /* Trotter steps */
	size_t ntsteps;
};

int circ_trott_data_init(struct circ_trott_data *cd, size_t ntstep);

/* Convenience function. Use instead of circ_trott_data_init. */
int circ_trott_data_init_from_file(
	struct circ_trott_data *cd, size_t num_steps, data_id fid);

void circ_trott_data_destroy(struct circ_trott_data *cd);

int circ_trott_simulate(const struct circ_trott_data *cd);

/* Circuit: qdrift */
struct circ_qdrift_data {
	struct circ_hamil hamil;
	struct circ_multidet multidet;

	double step_size;
	size_t depth;

	double *samples[2];
	size_t nsamples;
};

int circ_qdrift_data_init(struct circ_qdrift_data *cd, data_id fid);

void circ_qdrift_data_destroy(struct circ_qdrift_data *cd);

int circ_qdrift_simulate(const struct circ_qdrift_data *cd);

#ifdef __cplusplus
}
#endif

#endif // CIRC_H
