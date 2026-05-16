/*
 * circ_cache.c - Pauli rotation batch cache.
 *
 * Hamiltonian terms that share the same hi-qubit Pauli code
 * can reuse a single MPI amplitude exchange: the expensive
 * all-to-all communication depends only on the hi part.
 * This cache collects consecutive terms with identical hi
 * codes and batches their lo-qubit rotations into one
 * exchange + multi-rotation call.
 *
 * When a term with a different hi code arrives, or the cache
 * is full (CACHE_MAX = 1024 terms), the caller must flush
 * the cache (triggering the MPI exchange and rotations) and
 * then re-insert.  Overflow returns -1 to signal this.
 *
 * Each cache is its own instance (no file-static state),
 * so a process can hold several active circ contexts in
 * parallel -- useful for Python plugins and multi-circuit
 * drivers.
 */
#define LOG_SUBSYS "cache"

#include "c23_compat.h"
#include <stdlib.h>

#include "log.h"

#include "circ_cache.h"

#define CACHE_MAX UINT64_C(0x0400)

struct circ_cache {
	int qb_hi, qb_lo;
	struct paulis pa_hi;
	struct paulis pa_lo[CACHE_MAX];
	double phis[CACHE_MAX];
	size_t ch_len;
};

struct circ_cache *circ_cache_init(int hi, int lo)
{
	struct circ_cache *c = calloc(1, sizeof *c);
	if (!c)
		return nullptr;
	c->qb_hi = hi;
	c->qb_lo = lo;
	return c;
}

void circ_cache_free(struct circ_cache *c)
{
	free(c);
}

int circ_cache_insert(struct circ_cache *c,
	struct paulis pa, double phi)
{
	struct paulis lo, hi;
	paulis_split(pa, c->qb_lo, c->qb_hi, &lo, &hi);
	if (c->ch_len == 0) {
		c->pa_hi = hi;
		c->pa_lo[0] = lo;
		c->phis[0] = phi;
		c->ch_len = 1;
		return 0;
	}

	if (c->ch_len < CACHE_MAX && paulis_eq(c->pa_hi, hi)) {
		const size_t k = c->ch_len++;
		c->pa_lo[k] = lo;
		c->phis[k] = phi;
		return 0;
	}

	return -1;
}

void circ_cache_flush(struct circ_cache *c,
	circ_cache_op fn, void *data)
{
	if (c->ch_len > 0 && fn) {
		log_trace("flush: %zu terms (cache_max=%llu)", c->ch_len,
			(unsigned long long)CACHE_MAX);
		fn(c->pa_hi, c->pa_lo, c->phis, c->ch_len, data);
	}

	c->ch_len = 0;
}

size_t circ_cache_len(const struct circ_cache *c)
{
	return c->ch_len;
}
