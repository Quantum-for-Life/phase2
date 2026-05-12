/*
 * combinations - lex-ascending k-subset enumerator.
 *
 * Public API:
 *   void combinations_init(struct combo *c, uint32_t n, uint32_t k);
 *   int  combinations_next(struct combo *c, uint32_t *out);
 *
 * Algorithm L from Knuth, TAOCP vol. 4A, sec. 7.2.1.3.
 *
 * Internal state stores the current tuple in `state[0..k)` plus
 * two sentinels (`state[k] = n`, `state[k+1] = 0`) used by the
 * step procedure.  The `done` flag tracks exhaustion.  Calling
 * `combinations_next` after the iterator is exhausted is
 * idempotent (no state change, returns 1).
 *
 * Algorithm reference: phase2/doc/state-prep.md "Slater-Condon
 * expansion algorithm".
 */

#include "c23_compat.h"
#include <stdint.h>
#include <string.h>

#include "phase2/combinations.h"

void combinations_init(struct combo *c, const uint32_t n, const uint32_t k)
{
	memset(c, 0, sizeof *c);
	c->n = n;
	c->k = k;
	c->done = 0;
	c->started = 0;

	if (k > n || k > COMBINATIONS_MAX_K) {
		c->done = 1;
		return;
	}

	for (uint32_t i = 0; i < k; i++)
		c->state[i] = i;
	c->state[k] = n;
	c->state[k + 1] = 0;
}

int combinations_next(struct combo *c, uint32_t *out)
{
	if (c->done)
		return 1;

	const uint32_t k = c->k;
	const uint32_t n = c->n;

	if (!c->started) {
		c->started = 1;
		for (uint32_t i = 0; i < k; i++)
			out[i] = c->state[i];

		if (k == 0) {
			c->done = 1;
			return 0;
		}
		if (k == n) {
			c->done = 1;
			return 0;
		}

		return 0;
	}

	/*
	 * Lex step (matches Python itertools.combinations):
	 * find the largest j such that state[j] < n - (k - 1 - j),
	 * i.e. the rightmost position that can still increase
	 * leaving room for k-1-j higher values.  Increment that
	 * position and cascade the suffix to consecutive values.
	 */
	int32_t j = (int32_t)k - 1;
	while (j >= 0 && c->state[j] == n - k + (uint32_t)j)
		j--;
	if (j < 0) {
		c->done = 1;
		return 1;
	}
	c->state[j]++;
	for (uint32_t i = (uint32_t)j + 1; i < k; i++)
		c->state[i] = c->state[i - 1] + 1;

	for (uint32_t i = 0; i < k; i++)
		out[i] = c->state[i];

	return 0;
}
