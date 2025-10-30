#ifndef CIRC_CACHE
#define CIRC_CACHE

#include <stddef.h>

#include "phase2/paulis.h"

int circ_cache_init(int hi, int lo);

int circ_cache_insert(struct paulis pa, double phi);

typedef void (*circ_cache_op)(
	struct paulis, const struct paulis *, double *, size_t, void *);

void circ_cache_flush(circ_cache_op fn, void *data);

size_t circ_cache_len(void);

#endif /* CIRC_CACHE */
