#include "c23_compat.h"
#include <stdlib.h>

#include "circ_cache.h"

struct circ_cache *circ_cache_new(uint32_t qb_lo, uint32_t qb_hi)
{
	struct circ_cache *ch = malloc(sizeof *ch);
	if (!ch)
		goto err_alloc;

	struct paulis *lo =
		malloc(sizeof(struct paulis) * CIRC_CACHE_CODES_MAX);
	if (!lo)
		goto err_lo;
	double *angles = malloc(sizeof(double) * CIRC_CACHE_CODES_MAX);
	if (!angles)
		goto err_angles;

	ch->codes_lo = lo;
	ch->phis = angles;
	ch->qb_lo = qb_lo;
	ch->qb_hi = qb_hi;
	ch->len = 0;

	return ch;

	// free(phs);
err_angles:
	free(lo);
err_lo:
err_alloc:
	return nullptr;
}

void circ_cache_free(struct circ_cache *ch)
{
	if (ch) {
		free(ch->codes_lo);
		free(ch->phis);
		free(ch);
	}
}

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi)
{
	struct paulis lo, hi;
	paulis_split(code, ch->qb_lo, ch->qb_hi, &lo, &hi);
	if (ch->len == 0) {
		ch->code_hi = hi;
		ch->codes_lo[0] = lo;
		ch->phis[0] = phi;
		ch->len = 1;
		return 0;
	}

	if (ch->len < CIRC_CACHE_CODES_MAX && paulis_eq(ch->code_hi, hi)) {
		const size_t k = ch->len++;
		ch->codes_lo[k] = lo;
		ch->phis[k] = phi;
		return 0;
	}

	return -1;
}

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data)
{
	if (ch->len > 0 && op)
		op(ch->code_hi, ch->codes_lo, ch->phis, ch->len, data);
	ch->len = 0;
}
