/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stddef.h>
#include <stdint.h>

#include "xoshiro256ss.h"

static uint64_t splitmix64(uint64_t *st)
{
	uint64_t rt = (*st += UINT64_C(0x9E3779B97f4A7C15));
	rt	    = (rt ^ (rt >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	rt	    = (rt ^ (rt >> 27)) * UINT64_C(0x94D049BB133111EB);

	return rt ^ (rt >> 31);
}

void xoshiro256ss_init(struct xoshiro256ss *st, const uint64_t seed)
{
	uint64_t splmx = seed;
	st->s[0]       = splitmix64(&splmx);
	st->s[1]       = splitmix64(&splmx);
	st->s[2]       = splitmix64(&splmx);
	st->s[3]       = splitmix64(&splmx);

	for (size_t _ = 0; _ < 128; _++)
		xoshiro256ss_next(st);
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

uint64_t xoshiro256ss_next(struct xoshiro256ss *st)
{
	const uint64_t result = rotl(st->s[1] * 5, 7) * 9;
	const uint64_t t      = st->s[1] << 17;
	st->s[2] ^= st->s[0];
	st->s[3] ^= st->s[1];
	st->s[1] ^= st->s[2];
	st->s[0] ^= st->s[3];
	st->s[2] ^= t;
	st->s[3] = rotl(st->s[3], 45);

	return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoshiro256ss_jump(struct xoshiro256ss *st)
{
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
		0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0;
	for (size_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for (int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= st->s[0];
				s1 ^= st->s[1];
				s2 ^= st->s[2];
				s3 ^= st->s[3];
			}
			xoshiro256ss_next(st);
		}

	st->s[0] = s0;
	st->s[1] = s1;
	st->s[2] = s2;
	st->s[3] = s3;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void xoshiro256ss_longjump(struct xoshiro256ss *st)
{
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
		0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0;
	for (size_t i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for (int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= st->s[0];
				s1 ^= st->s[1];
				s2 ^= st->s[2];
				s3 ^= st->s[3];
			}
			xoshiro256ss_next(st);
		}

	st->s[0] = s0;
	st->s[1] = s1;
	st->s[2] = s2;
	st->s[3] = s3;
}
