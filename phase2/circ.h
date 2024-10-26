#ifndef CIRC_IMPL_H
#define CIRC_IMPL_H

#include "phase2/circ.h"

#ifdef __cplusplus
extern "C" {
#endif

int circ_res_init(struct circ *c);

void circ_res_destroy(struct circ *c);

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t n;
};

int circ_cache_init(struct circ_cache *ch, size_t qb_lo, size_t qb_hi);

void circ_cache_destroy(struct circ_cache *ch);

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi);

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data);

#ifdef __cplusplus
}
#endif

#endif // CIRC_IMPL_H
