#ifndef PHASE2_COMBINATIONS_H
#define PHASE2_COMBINATIONS_H

#include <stdint.h>

/*
 * combinations - lex-ascending k-subset enumerator.
 *
 * Iterator over all k-subsets of {0, 1, ..., n-1} in
 * lexicographic order (Knuth, TAOCP vol. 4A, sec. 7.2.1.3,
 * algorithm L).
 *
 * Public API:
 *   void combinations_init(struct combo *c, uint32_t n, uint32_t k);
 *   int  combinations_next(struct combo *c, uint32_t *out);
 *
 * `out` must point to an array of at least k uint32_t values.
 * `combinations_next` writes the next subset into `out[0..k)` and
 * returns 0 (more subsets to come) or 1 (sequence exhausted).
 * After the iterator is exhausted, subsequent calls keep
 * returning 1 idempotently and leave `out` unchanged.
 *
 * Index i of a returned tuple equals its rank in the lex
 * sequence; this lets `state_prep_coeff_expand` index a
 * precomputed determinant table by tuple rank.
 *
 * Edge cases:
 *   k == 0 yields exactly one empty tuple, then done.
 *   k >  n yields zero tuples (immediate done).
 *
 * The maximum supported subset size is COMBINATIONS_MAX_K.
 *
 * Algorithm reference: phase2/doc/state-prep.md, "Slater-Condon
 * expansion algorithm".
 */

#define COMBINATIONS_MAX_K (32u)

struct combo {
	uint32_t n;
	uint32_t k;
	uint32_t state[COMBINATIONS_MAX_K + 2];
	int done;
	int started;
};

void combinations_init(struct combo *c, uint32_t n, uint32_t k);
int combinations_next(struct combo *c, uint32_t *out);

#endif /* PHASE2_COMBINATIONS_H */
