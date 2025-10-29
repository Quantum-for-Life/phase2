#ifndef CIRC_CACHE
#define CIRC_CACHE

#include <stddef.h>

#include "phase2/paulis.h"

#define CIRC_CACHE_CODES_MAX UINT64_C(0x0400)

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t len;
};

struct circ_cache *circ_cache_new(uint32_t qb_lo, uint32_t qb_hi);
void circ_cache_free(struct circ_cache *ch);
int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi);
void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data);

#endif /* CIRC_CACHE */
