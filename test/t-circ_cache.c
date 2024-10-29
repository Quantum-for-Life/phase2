#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0x77e8fe9b90caf912)
static struct world WD;

#define SEED UINT64_C(0x334b06fc8c7b40ea)
static struct xoshiro256ss RNG;

static double rand_double(void)
{
	const uint64_t x = xoshiro256ss_next(&RNG);

	return (x >> 11) * 0x1.0p-53;
}

static _Complex double rand_complex(void)
{
	const _Complex double z = CMPLX(rand_double(), rand_double());

	return z / cabs(z);
}

static int rand_pauli(void)
{
	const int x = xoshiro256ss_next(&RNG) % 4;

	return x;
}

static struct paulis rand_paulis(uint32_t qb)
{
	struct paulis ps = paulis_new();
	for (uint32_t k = 0; k < qb; k++)
		paulis_set(&ps, rand_pauli(), k);

	return ps;
}

void t_cache_00(void)
{
	struct paulis p = paulis_new();
	struct circ_cache ch;

	TEST_ASSERT(circ_cache_init(&ch, 1, 1) == 0, "can't init cache");
	TEST_ASSERT(ch.n == 0, "init size must be zero");
	TEST_ASSERT(circ_cache_insert(&ch, p, 0.0) == 0,
		"cannot insert first code");
	TEST_ASSERT(ch.n == 1, "cache size should be 1");
	TEST_ASSERT(
		circ_cache_insert(&ch, p, 0.1) == 0, "insert the same code");
	TEST_ASSERT(ch.n == 2, "cache size should be 2");
	paulis_set(&p, PAULI_X, 1);
	TEST_ASSERT(
		circ_cache_insert(&ch, p, 0.0) < 1, "insert different code");
	TEST_ASSERT(ch.n == 2, "cache size should be 2");
	circ_cache_flush(&ch, nullptr, nullptr);
	TEST_ASSERT(ch.n == 0, "cache size should be 0");
	TEST_ASSERT(
		circ_cache_insert(&ch, p, 0.0) == 0, "insert different code");
	TEST_ASSERT(ch.n == 1, "cache size should be 1");

	circ_cache_destroy(&ch);
}

void TEST_MAIN(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	t_cache_00();

	world_destroy();
}
