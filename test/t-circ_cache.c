#define LOG_SUBSYS "test"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "log.h"
#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "circ_cache.h"

#include "test.h"

#define WD_SEED UINT64_C(0x77e8fe9b90caf912)
static struct world_info WD;

#define SEED UINT64_C(0x334b06fc8c7b40ea)
static struct xoshiro256ss RNG;

/*
 * Basic insert/flush/len semantics: an empty cache
 * accepts a first code, accepts subsequent codes with
 * the same hi part, rejects a different hi part with
 * -1, and after flush is back to empty.
 */
void t_cache_basic(void)
{
	struct paulis p = paulis_new();

	struct circ_cache *c = circ_cache_init(1, 1);
	TEST_ASSERT(c, "can't init cache");
	TEST_ASSERT(circ_cache_len(c) == 0, "init size must be zero");
	TEST_ASSERT(circ_cache_insert(c, p, 0.0) == 0,
		"cannot insert first code");
	TEST_ASSERT(circ_cache_len(c), "cache size should be 1");
	TEST_ASSERT(circ_cache_insert(c, p, 0.1) == 0, "insert the same code");
	TEST_ASSERT(circ_cache_len(c) == 2, "cache size should be 2");
	paulis_set(&p, PAULI_X, 1);
	TEST_ASSERT(circ_cache_insert(c, p, 0.0) < 0, "insert different code");
	TEST_ASSERT(circ_cache_len(c) == 2, "cache size should be 2");
	circ_cache_flush(c, nullptr, nullptr);
	TEST_ASSERT(circ_cache_len(c) == 0, "cache size should be 0");
	TEST_ASSERT(circ_cache_insert(c, p, 0.0) == 0,
		"insert after flush should succeed");
	TEST_ASSERT(circ_cache_len(c) == 1, "cache size should be 1");
	circ_cache_free(c);
}

/*
 * Boundary: CACHE_MAX inserts (all matching hi)
 * succeed; the next insert returns -1; a flush
 * empties the cache and the rejected insert
 * succeeds.  The exact cap (1024) is private to
 * circ_cache.c, so the test discovers it by
 * pushing until the first -1.
 */
#define CACHE_MAX_EXPECTED 1024
void t_cache_overflow(void)
{
	struct paulis p = paulis_new();
	struct circ_cache *c = circ_cache_init(1, 1);
	TEST_ASSERT(c, "can't init cache");

	size_t inserted = 0;
	while (circ_cache_insert(c, p, (double)inserted) == 0)
		inserted++;
	TEST_EQ(inserted, (size_t)CACHE_MAX_EXPECTED);
	TEST_EQ(circ_cache_len(c), (size_t)CACHE_MAX_EXPECTED);

	/* One more insert: refused. */
	TEST_ASSERT(circ_cache_insert(c, p, -1.0) < 0,
		"full cache must reject insert");
	TEST_EQ(circ_cache_len(c), (size_t)CACHE_MAX_EXPECTED);

	/* Flush -> empty -> next insert succeeds. */
	circ_cache_flush(c, nullptr, nullptr);
	TEST_EQ(circ_cache_len(c), (size_t)0);
	TEST_ASSERT(circ_cache_insert(c, p, -1.0) == 0,
		"insert after flush must succeed");

	circ_cache_free(c);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	t_cache_basic();
	t_cache_overflow();

	world_free();
}
