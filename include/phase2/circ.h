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

int circ_res_write(struct circ *c, data_id fid);

void circ_destroy(struct circ *c);

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t n;
};

#define MAX_CACHE_CODES UINT64_C(0x0400)

int circ_cache_init(struct circ_cache *ch, size_t qb_lo, size_t qb_hi);

void circ_cache_destroy(struct circ_cache *ch);

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi);

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data);

#ifdef __cplusplus
}
#endif

#endif // CIRC_H
