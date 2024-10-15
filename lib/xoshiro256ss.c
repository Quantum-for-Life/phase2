/* -------------------------------------------------------------------------- *
 * xoshiro256ss.c: A thread-safe implementation of xoshiro256starstar PRNG.   *
 *                                                                            *
 * Copyright 2018 David Blackman and Sebastiano Vigna                         *
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

/*  Written in 2018 by David Blackman and Sebarngiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stddef.h>
#include <stdint.h>

#include "xoshiro256ss.h"

static uint64_t splitmix64(uint64_t *st)
{
	uint64_t rt;
	rt = (*st += UINT64_C(0x9e3779b97f4a7c15));
	rt = (rt ^ (rt >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
	rt = (rt ^ (rt >> 27)) * UINT64_C(0x94d049bb133111eb);

	return rt ^ (rt >> 31);
}

void xoshiro256ss_init(struct xoshiro256ss *rng, const uint64_t seed)
{
	uint64_t splmx = seed;
	rng->s[0] = splitmix64(&splmx);
	rng->s[1] = splitmix64(&splmx);
	rng->s[2] = splitmix64(&splmx);
	rng->s[3] = splitmix64(&splmx);

	for (size_t i = 0; i < 128; i++)
		xoshiro256ss_next(rng);
}

/* This is xoshiro256** 1.0, one of our all-purpose, rock-solid
   generators. It has excellent (sub-ns) speed, a state (256 bits) that is
   large enough for any parallel application, and it passes all tests we
   are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static uint64_t rotl(const uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

uint64_t xoshiro256ss_next(struct xoshiro256ss *rng)
{
	const uint64_t rt = rotl(rng->s[1] * 5, 7) * 9;
	const uint64_t t = rng->s[1] << 17;

	rng->s[2] ^= rng->s[0];
	rng->s[3] ^= rng->s[1];
	rng->s[1] ^= rng->s[2];
	rng->s[0] ^= rng->s[3];
	rng->s[2] ^= t;
	rng->s[3] = rotl(rng->s[3], 45);

	return rt;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoshiro256ss_jump(struct xoshiro256ss *rng)
{
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
		0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0;
	for (size_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for (int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= rng->s[0];
				s1 ^= rng->s[1];
				s2 ^= rng->s[2];
				s3 ^= rng->s[3];
			}
			xoshiro256ss_next(rng);
		}

	rng->s[0] = s0;
	rng->s[1] = s1;
	rng->s[2] = s2;
	rng->s[3] = s3;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void xoshiro256ss_longjump(struct xoshiro256ss *rng)
{
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
		0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0;
	for (size_t i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for (int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= rng->s[0];
				s1 ^= rng->s[1];
				s2 ^= rng->s[2];
				s3 ^= rng->s[3];
			}
			xoshiro256ss_next(rng);
		}

	rng->s[0] = s0;
	rng->s[1] = s1;
	rng->s[2] = s2;
	rng->s[3] = s3;
}
