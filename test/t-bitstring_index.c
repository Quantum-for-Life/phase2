#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xa11b7c4422bb7d1c)

static inline uint64_t occ_pair_to_idx(const uint32_t *occ_a, uint32_t na,
	const uint32_t *occ_b, uint32_t nb, uint32_t n_sites)
{
	uint64_t idx = 0;
	for (uint32_t i = 0; i < na; i++)
		idx |= (UINT64_C(1) << occ_a[i]);
	for (uint32_t i = 0; i < nb; i++)
		idx |= (UINT64_C(1) << (n_sites + occ_b[i]));
	return idx;
}

static inline uint64_t drop_two_bits(uint64_t idx, uint32_t n_sites)
{
	uint64_t lo = (idx >> 1) & ((UINT64_C(1) << (n_sites - 1)) - 1);
	uint64_t hi = idx >> (n_sites + 1);
	return lo | (hi << (n_sites - 1));
}

static void t_pack_unpack(uint32_t n_sites)
{
	const uint32_t occ_a[3] = { 0, 1, 2 };
	const uint32_t occ_b[3] = { 0, 1, 2 };
	const uint64_t idx =
		occ_pair_to_idx(occ_a, 2, occ_b, 1, n_sites);
	const uint64_t expect = (UINT64_C(1) << 0) | (UINT64_C(1) << 1)
				| (UINT64_C(1) << (n_sites + 0));
	TEST_EQ(idx, expect);
}

static void t_empty(void)
{
	const uint64_t idx = occ_pair_to_idx(NULL, 0, NULL, 0, 4);
	TEST_EQ(idx, UINT64_C(0));
}

static void t_drop_round_trip(void)
{
	const uint32_t n_sites = 8;
	const uint32_t occ_a[3] = { 0, 3, 5 };
	const uint32_t occ_b[3] = { 0, 2, 6 };
	const uint64_t full = occ_pair_to_idx(occ_a, 3, occ_b, 3, n_sites);
	const uint64_t mask = UINT64_C(1) | (UINT64_C(1) << n_sites);
	TEST_ASSERT((full & mask) == mask,
		"bit 0 and bit n_sites should be set, idx=0x%lx",
		(unsigned long)full);

	const uint64_t dropped = drop_two_bits(full, n_sites);
	const uint32_t total_bits = 2 * (n_sites - 1);
	TEST_ASSERT((dropped >> total_bits) == 0,
		"dropped should fit in 2*(n_sites-1) bits");

	uint64_t reconstructed = 0;
	const uint64_t low = dropped & ((UINT64_C(1) << (n_sites - 1)) - 1);
	const uint64_t high = dropped >> (n_sites - 1);
	reconstructed = (low << 1) | UINT64_C(1)
			| (high << (n_sites + 1))
			| (UINT64_C(1) << n_sites);
	TEST_ASSERT(reconstructed == full,
		"reconstruction mismatch: got 0x%lx vs 0x%lx",
		(unsigned long)reconstructed, (unsigned long)full);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_empty();
	t_pack_unpack(4);
	t_pack_unpack(8);
	t_pack_unpack(14);
	t_pack_unpack(18);

	t_drop_round_trip();

	world_free();
}
