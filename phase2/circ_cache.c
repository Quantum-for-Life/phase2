#include "c23_compat.h"
#include <stdlib.h>

#include "circ_cache.h"

#define CACHE_MAX UINT64_C(0x0400)

static int qb_hi, qb_lo;
static struct paulis pa_hi, pa_lo[CACHE_MAX];
static double phis[CACHE_MAX];
static size_t ch_len;

int circ_cache_init(int hi, int lo)
{
	qb_hi = hi;
	qb_lo = lo;
	ch_len = 0;

	return 0;
}

int circ_cache_insert(struct paulis pa, double phi)
{
	struct paulis lo, hi;
	paulis_split(pa, qb_lo, qb_hi, &lo, &hi);
	if (ch_len == 0) {
		pa_hi = hi;
		pa_lo[0] = lo;
		phis[0] = phi;
		ch_len = 1;
		return 0;
	}

	if (ch_len < CACHE_MAX && paulis_eq(pa_hi, hi)) {
		const size_t k = ch_len++;
		pa_lo[k] = lo;
		phis[k] = phi;
		return 0;
	}

	return -1;
}

void circ_cache_flush(circ_cache_op fn, void *data)
{
	if (ch_len > 0 && fn)
		fn(pa_hi, pa_lo, phis, ch_len, data);

	ch_len = 0;
}

size_t circ_cache_len(void)
{
	return ch_len;
}
