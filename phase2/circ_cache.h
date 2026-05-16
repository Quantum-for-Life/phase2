#ifndef CIRC_CACHE
#define CIRC_CACHE

#include <stddef.h>

#include "phase2/paulis.h"

/*
 * Pauli-rotation batch cache.  An instance accumulates
 * Hamiltonian terms that share a single hi-qubit Pauli
 * code; flushing emits one MPI exchange + multi-rotation
 * for all batched terms.
 *
 * The full struct definition lives in circ_cache.c;
 * callers (currently only phase2/circ.c) hold an opaque
 * pointer and go through the API below.
 */
struct circ_cache;

/* Allocate a new cache for the given (hi, lo) qubit
 * partition.  Returns NULL on allocation failure. */
struct circ_cache *circ_cache_init(int hi, int lo);

void circ_cache_free(struct circ_cache *c);

/* Try to insert (pa, phi).  Returns 0 on success, -1 if
 * the cache is full or `pa`'s hi-code differs from the
 * current batch -- caller must flush then retry. */
int circ_cache_insert(struct circ_cache *c,
	struct paulis pa, double phi);

/* Flush the cache: invoke fn with the accumulated hi
 * code, lo codes, phis, term count, and caller-supplied
 * data; then reset to empty.  Calling on an empty cache
 * is a no-op. */
void circ_cache_flush(struct circ_cache *c,
	void (*fn)(struct paulis hi, const struct paulis *lo,
		double *phis, size_t n, void *data),
	void *data);

size_t circ_cache_len(const struct circ_cache *c);

#endif /* CIRC_CACHE */
