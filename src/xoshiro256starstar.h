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

#ifndef XOSHIRO256STARSTAR_ORIG_H
#define XOSHIRO256STARSTAR_ORIG_H

#include <stdint.h>

struct xoshiro256starstar {
	uint16_t s[4];
};

void
xoshiro256starstar_init(struct xoshiro256starstar *state, uint64_t seed);

uint64_t
xoshiro256starstar_next(struct xoshiro256starstar *state);

void
xoshiro256starstar_jump(struct xoshiro256starstar *state);

void
xoshiro256starstar_long_jump(struct xoshiro256starstar *state);

#endif // XOSHIRO256STARSTAR_ORIG_H
