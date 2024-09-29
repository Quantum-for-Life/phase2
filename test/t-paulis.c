#include "paulis.h"
#include "xoshiro256ss.h"

#include "test.h"

#define SEED (0x2a44fe101UL)
static struct xoshiro256ss RNG;

#define WIDTH (64)

void test_paulis_new(void)
{
	struct paulis ps = paulis_new();

	for (size_t k = 0; k < WIDTH; k++)
		TEST_ASSERT(paulis_get(ps, k) == PAULI_I,
			"k=%zu, should be PAULI_I", k);
}

void test_paulis_getset(size_t tag)
{
	enum pauli_op op[WIDTH];
	struct paulis ps = paulis_new();

	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, op[k] = (int)(xoshiro256ss_next(&RNG) % 4), k);

	for (size_t k = 0; k < WIDTH; k++) {
		enum pauli_op p = paulis_get(ps, k);
		TEST_ASSERT(p == op[k],
			"[%zu] k=%zu, pauli=%d, expected=%d", tag, k, p, op[k]);
	}
}

void test_paulis_eq(size_t tag)
{
	enum pauli_op op;
	struct paulis ps1, ps2;

	for (size_t k = 0; k < WIDTH; k++) {
		op = (int)(xoshiro256ss_next(&RNG) % 4);
		paulis_set(&ps1, op, k);
		paulis_set(&ps2, op, k);
	}
	TEST_ASSERT(paulis_eq(ps1, ps2), "[%zu] should be equal", tag);

	paulis_set(&ps1, op == PAULI_I ? PAULI_X : PAULI_I, WIDTH - 1);
	TEST_ASSERT(!paulis_eq(ps1, ps2), "[%zu] should not be equal", tag);
}


void test_paulis_shr(size_t n)
{
	enum pauli_op op;
	struct paulis ps1, ps2;

	ps1 = ps2 = paulis_new();
	for (size_t k = 0; k < WIDTH; k++) {
		op = (int)(xoshiro256ss_next(&RNG) % 3 + 1); /* no PAULI_I */
		paulis_set(&ps1, op, k);
		if (k >= n)
			paulis_set(&ps2, op, k - n);
	}
	if (n > 0)
		TEST_ASSERT(!paulis_eq(ps1, ps2),
			"n=%zu, should not be equal", n);

	paulis_shr(&ps1, n);
	TEST_ASSERT(paulis_eq(ps1, ps2), "n=%zu, should be equal", n);
}

void TEST_MAIN(void)
{
	xoshiro256ss_init(&RNG, SEED);

	test_paulis_new();
	for (size_t n = 0; n < 100; n++) {
		test_paulis_eq(n);
		test_paulis_getset(n);
		test_paulis_shr(n);
	}
}
