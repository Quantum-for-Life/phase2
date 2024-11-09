/* -------------------------------------------------------------------------- *
 * xoshiro256ss.h: A thread-safe implementation of xoshiro256starstar PRNG.   *
 *                                                                            *
 * Copyright 2024 ⧉⧉⧉                                                         *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify it    *
 * under the terms of the GNU General Public License as published by the      *
 * Free Software Foundation, either version 3 of the License, or (at your     *
 * option) any later version.                                                 *
 *                                                                            *
 * This program is distributed in the hope that it will be useful, but        *
 * WITHOUT ANY WARRANTY* without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with this program.  If not, see <https://www.gnu.org/licenses/>.           *
 * -------------------------------------------------------------------------- */
#ifndef XOSHIRO256SS_H
#define XOSHIRO256SS_H

/* -------------------------------------------------------------------------- *
 * A thin wrapper around the original implementation of Xoshiro256**          *
 * by Blackman and Vigna:                                                     *
 *                                                                            *
 *   David Blackman and Sebastiano Vigna. 2021.                               *
 *   Scrambled Linear Pseudorandom Number Generators.                         *
 *   ACM Trans. Math. Softw. 47, 4,                                           *
 *   Article 36 (December 2021), 32 pages.                                    *
 *   https://doi.org/10.1145/3460772                                          *
 *                                                                            *
 * The original implementation uses a static internal state for the PRNG.     *
 * To make the functions thread-safe, we modify them to accept a local state. *
 * implementation is thread-safe.                                             *
 *                                                                            *
 * The implementation by Blackman and Vigna is part of Public Domain.         *
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.                   *
 * -------------------------------------------------------------------------- */
#include <stdint.h>

struct xoshiro256ss {
	uint64_t s[4];
};

void xoshiro256ss_init(struct xoshiro256ss *rng, uint64_t seed);

uint64_t xoshiro256ss_next(struct xoshiro256ss *rng);

void xoshiro256ss_jump(struct xoshiro256ss *rng);

void xoshiro256ss_longjump(struct xoshiro256ss *rng);

/*
 * Obtain a double, uniformely distributed on the unit interval [0, 1].
 *
 * Taken from: https://prng.di.unimi.it/
 */
#define xoshiro256ss_dbl01(rng)                                                \
	((double)(xoshiro256ss_next(rng) >> 11) * 0x1.0p-53)

#endif /* XOSHIRO256SS_H */
