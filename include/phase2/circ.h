#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/paulis.h"

#ifdef __cplusplus
extern "C" {
#endif

struct circ_hamil {
	size_t nqb;
	struct {
		struct paulis op;
		double cf;
	} *terms;
	size_t nterms;
};

struct circ_muldet {
	struct {
		uint64_t idx;
		_Complex double cf;
	} *dets;
	size_t ndets;
};

struct circ {
	struct circ_hamil hamil;
	struct circ_muldet muldet;

	void *data;
	void *res;
};

int circ_init(struct circ *c, data_id fid, void *data);

int circ_simulate(struct circ *c);

int circ_res_init(struct circ *c);

void circ_res_destroy(struct circ *c);

int circ_res_write(struct circ *c, data_id fid);

void circ_destroy(struct circ *c);

#ifdef __cplusplus
}
#endif

#endif // CIRC_H
