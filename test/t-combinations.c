#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/combinations.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x4d12c3ee9988b4a7)

static uint64_t binom(uint32_t n, uint32_t k)
{
	if (k > n)
		return 0;
	if (k > n - k)
		k = n - k;
	uint64_t num = 1, den = 1;
	for (uint32_t i = 1; i <= k; i++) {
		num *= (n - k + i);
		den *= i;
	}
	return num / den;
}

static void t_lex_order(uint32_t n, uint32_t k)
{
	struct combo it;
	combinations_init(&it, n, k);

	const uint64_t expected = binom(n, k);
	uint64_t count = 0;
	uint32_t prev[COMBINATIONS_MAX_K];
	uint32_t buf[COMBINATIONS_MAX_K];
	int has_prev = 0;

	while (combinations_next(&it, buf) == 0) {
		for (uint32_t i = 1; i < k; i++)
			TEST_ASSERT(buf[i] > buf[i - 1],
				"non-ascending tuple at count=%lu",
				(unsigned long)count);
		if (has_prev) {
			int lt = 0;
			for (uint32_t i = 0; i < k; i++) {
				if (buf[i] > prev[i]) {
					lt = 1;
					break;
				}
				if (buf[i] < prev[i])
					break;
			}
			TEST_ASSERT(lt, "lex order broken at count=%lu",
				(unsigned long)count);
		}
		for (uint32_t i = 0; i < k; i++)
			prev[i] = buf[i];
		has_prev = 1;
		count++;
	}
	TEST_ASSERT(count == expected,
		"n=%u k=%u: count=%lu expected=%lu", n, k,
		(unsigned long)count, (unsigned long)expected);

	for (int x = 0; x < 3; x++)
		TEST_ASSERT(combinations_next(&it, buf) == 1,
			"post-done idempotent");
}

static void t_k_zero(void)
{
	struct combo it;
	combinations_init(&it, 5, 0);
	uint32_t buf[1] = { 0xdeadbeef };
	TEST_EQ(combinations_next(&it, buf), 0);
	TEST_EQ(combinations_next(&it, buf), 1);
	TEST_EQ(combinations_next(&it, buf), 1);
}

static void t_k_gt_n(void)
{
	struct combo it;
	combinations_init(&it, 3, 5);
	uint32_t buf[5];
	TEST_EQ(combinations_next(&it, buf), 1);
	TEST_EQ(combinations_next(&it, buf), 1);
}

static void t_k_eq_n(void)
{
	struct combo it;
	combinations_init(&it, 4, 4);
	uint32_t buf[4] = { 0 };
	TEST_EQ(combinations_next(&it, buf), 0);
	for (uint32_t i = 0; i < 4; i++)
		TEST_EQ(buf[i], i);
	TEST_EQ(combinations_next(&it, buf), 1);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_k_zero();
	t_k_gt_n();
	t_k_eq_n();

	t_lex_order(4, 2);
	t_lex_order(6, 3);
	t_lex_order(8, 4);
	t_lex_order(10, 5);
	t_lex_order(17, 8);
	t_lex_order(18, 9);

	world_free();
}
