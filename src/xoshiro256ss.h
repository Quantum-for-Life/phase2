/*
 * A thin wrapper around the original implementation of Xoshiro256**
 * by Blackman and Vigna:
 *
 *   David Blackman and Sebastiano Vigna. 2021.
 *   Scrambled Linear Pseudorandom Number Generators.
 *   ACM Trans. Math. Softw. 47, 4,
 *   Article 36 (December 2021), 32 pages.
 *   https://doi.org/10.1145/3460772
 *
 * Their implementation uses a static internal state for the PRNG.
 * We modify original functions to accept a local state, so that the
 * implementation is thread-safe.
 */

#ifndef XOSHIRO256SS_H
#define XOSHIRO256SS_H

#include <stdint.h>

struct xoshiro256ss {
	uint64_t s[4];
};

void xoshiro256ss_init(struct xoshiro256ss *st, uint64_t seed);

uint64_t xoshiro256ss_next(struct xoshiro256ss *st);

void xoshiro256ss_jump(struct xoshiro256ss *st);

void xoshiro256ss_longjump(struct xoshiro256ss *st);

#endif /* XOSHIRO256SS_H */
