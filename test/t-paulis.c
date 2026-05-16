#include "c23_compat.h"
#include <complex.h>
#include <stdint.h>

#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0x92ee46fade7a5742)

#define SEED (0x2a44fe101UL)
static struct xoshiro256ss RNG;

#define WIDTH (64)

static void t_paulis_new(void)
{
	struct paulis ps = paulis_new();

	for (size_t k = 0; k < WIDTH; k++)
		TEST_ASSERT(paulis_get(ps, k) == PAULI_I,
			"k=%zu, should be PAULI_I", k);
}

static void t_paulis_getset(size_t tag)
{
	enum pauli_op op[WIDTH];
	struct paulis ps = paulis_new();

	for (size_t k = 0; k < WIDTH; k++) {
		op[k] = (enum pauli_op)(xoshiro256ss_next(&RNG) % 4);
		paulis_set(&ps, op[k], k);
	}

	for (size_t k = 0; k < WIDTH; k++) {
		enum pauli_op p = paulis_get(ps, k);
		TEST_ASSERT(p == op[k], "[%zu] k=%zu, pauli=%d, expected=%d",
			tag, k, p, op[k]);
	}
}

static void t_paulis_eq(size_t tag)
{
	enum pauli_op op = PAULI_I;
	struct paulis ps1, ps2;

	for (size_t k = 0; k < WIDTH; k++) {
		op = (enum pauli_op)(xoshiro256ss_next(&RNG) % 4);
		paulis_set(&ps1, op, k);
		paulis_set(&ps2, op, k);
	}
	TEST_ASSERT(paulis_eq(ps1, ps2), "[%zu] should be equal", tag);

	paulis_set(&ps1, op == PAULI_I ? PAULI_X : PAULI_I, WIDTH - 1);
	TEST_ASSERT(!paulis_eq(ps1, ps2), "[%zu] should not be equal", tag);
}

static void t_paulis_shr(size_t n)
{
	enum pauli_op op;
	struct paulis ps1, ps2;

	ps1 = ps2 = paulis_new();
	for (size_t k = 0; k < WIDTH; k++) {
		op = (enum pauli_op)(xoshiro256ss_next(&RNG) % 3 + 1); /* no I */
		paulis_set(&ps1, op, k);
		if (k >= n)
			paulis_set(&ps2, op, k - n);
	}
	if (n > 0)
		TEST_ASSERT(
			!paulis_eq(ps1, ps2), "n=%zu, should not be equal", n);

	paulis_shr(&ps1, n);
	TEST_ASSERT(paulis_eq(ps1, ps2), "n=%zu, should be equal", n);
}

static void t_paulis_effect_00(void)
{
	_Complex double z = 11.3;
	struct paulis ps = paulis_new();

	for (size_t i = 0; i < 99; i++) {
		uint64_t x, y = xoshiro256ss_next(&RNG);
		x = paulis_effect(ps, y, &z);
		TEST_ASSERT(x == y, "[%zu] x=%lu, y=%lu", i, x, y);
		TEST_ASSERT(
			z == 11.3, "[%zu] z = %f + %fi", i, creal(z), cimag(z));
	}
}

static void t_paulis_effect_i1(void)
{
	_Complex double z;
	struct paulis ps = paulis_new();
	paulis_set(&ps, PAULI_Y, 0);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x00, &z), 0x01);
	TEST_EQ(z, I);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x01, &z), 0x00);
	TEST_EQ(z, -I);
}

static void t_paulis_effect_i2(void)
{
	_Complex double z;
	struct paulis ps = paulis_new();
	paulis_set(&ps, PAULI_Y, 0);
	paulis_set(&ps, PAULI_Y, 1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x00, &z), 0x03);
	TEST_EQ(z, -1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x01, &z), 0x02);
	TEST_EQ(z, 1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x02, &z), 0x01);
	TEST_EQ(z, 1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x03, &z), 0x00);
	TEST_EQ(z, -1);
}

static void t_paulis_effect_01(void)
{
	_Complex double z;
	struct paulis ps = paulis_new();

	paulis_set(&ps, PAULI_X, 0);
	paulis_set(&ps, PAULI_X, 1);
	paulis_set(&ps, PAULI_X, 3);

	TEST_EQ(paulis_effect(ps, 0x00, nullptr), 0x0B); /* 0xB = 0b1011 */
	TEST_EQ(paulis_effect(ps, 0x01, nullptr), 0x0A);
	TEST_EQ(paulis_effect(ps, 0x02, nullptr), 0x09);
	TEST_EQ(paulis_effect(ps, 0x03, nullptr), 0x08);

	paulis_set(&ps, PAULI_Z, 1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x00, &z), 0x09); /* 0x09 = 0x1001 */
	TEST_EQ(z, 1.0);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x01, &z), 0x08);
	TEST_EQ(z, 1.0);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x02, &z), 0x0B);
	TEST_EQ(z, -1.0);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x03, &z), 0x0A);
	TEST_EQ(z, -1.0);

	paulis_set(&ps, PAULI_Y, 1);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x00, &z), 0x0B); /* 0x09 = 0x1001 */
	TEST_EQ(z, I);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x01, &z), 0x0A);
	TEST_EQ(z, I);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x02, &z), 0x09);
	TEST_EQ(z, -I);

	z = 1.0;
	TEST_EQ(paulis_effect(ps, 0x03, &z), 0x08);
	TEST_EQ(z, -I);
}

static void t_paulis_effect_02(size_t tag)
{
	uint64_t x, y, y_exp, kk;
	_Complex double z, z_exp;
	enum pauli_op op;
	struct paulis ps = paulis_new();

	x = xoshiro256ss_next(&RNG);
	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, (enum pauli_op)(xoshiro256ss_next(&RNG) % 4), k);

	y_exp = 0;
	z_exp = z = 1.0;
	for (size_t k = 0; k < WIDTH; k++) {
		op = paulis_get(ps, k);
		kk = UINT64_C(1) << k;

		switch (op) {
		case PAULI_I:
			y_exp |= x & kk;
			break;
		case PAULI_X:
			y_exp |= ~x & kk;
			break;
		case PAULI_Y:
			y_exp |= ~x & kk;
			if (x & kk)
				z_exp *= -I;
			else
				z_exp *= I;
			break;
		case PAULI_Z:
			y_exp |= x & kk;
			if (x & kk)
				z_exp *= -1;
			break;
		}
	}

	y = paulis_effect(ps, x, &z);

	TEST_ASSERT(y == y_exp, "[%zu] y=0x%lx, y_exp=0x%lx", tag, y, y_exp);
	TEST_ASSERT(z == z_exp, "[%zu] z=%f+%fi, z_exp=%f+%fi", tag, creal(z),
		cimag(z), creal(z_exp), cimag(z_exp));
}

static void t_paulis_split_01(size_t tag)
{
	uint32_t lo, hi;
	struct paulis ps_lo, ps_hi, ps = paulis_new();
	enum pauli_op op_lo, op_hi, op_ex;

	lo = xoshiro256ss_next(&RNG) % (WIDTH - 1);
	hi = xoshiro256ss_next(&RNG) % (WIDTH - lo);
	TEST_ASSERT(lo + hi <= WIDTH, "wrong test params");

	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, (enum pauli_op)(xoshiro256ss_next(&RNG) % 4), k);

	paulis_split(ps, lo, hi, &ps_lo, &ps_hi);

	for (uint32_t k = 0; k < lo; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);
		op_ex = paulis_get(ps, k);

		TEST_ASSERT(op_lo == op_ex,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d, op_ex=%d", tag,
			lo, hi, k, op_lo, op_ex);
		TEST_ASSERT(op_hi == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_hi=%d", tag, lo, hi, k,
			op_hi);
	}
	for (uint32_t k = lo; k < hi; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);
		op_ex = paulis_get(ps, k);

		TEST_ASSERT(op_lo == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d", tag, lo, hi, k,
			op_lo);
		TEST_ASSERT(op_hi == op_ex,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d, op_ex=%d", tag,
			lo, hi, k, op_hi, op_ex);
	}
	for (uint32_t k = lo + hi; k < WIDTH; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);

		TEST_ASSERT(op_lo == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d", tag, lo, hi, k,
			op_lo);
		TEST_ASSERT(op_hi == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_hi=%d", tag, lo, hi, k,
			op_hi);
	}
}

static void t_paulis_cmp_00()
{
	const struct paulis a = paulis_new();

	struct paulis b = paulis_new();
	paulis_set(&b, PAULI_X, 1);

	TEST_ASSERT(paulis_cmp(a, a) == 0, "cmp with itself");
	TEST_ASSERT(paulis_cmp(b, b) == 0, "cmp with itself");

	TEST_ASSERT(paulis_cmp(a, b) == -1, "cmp with I");
	TEST_ASSERT(paulis_cmp(b, a) == 1, "cmp with I");
}

static int cmp_explicit(const struct paulis a, const struct paulis b)
{
	for (uint32_t n = 0; n < WIDTH; n++) {
		const int x = paulis_get(a, WIDTH - n - 1);
		const int y = paulis_get(b, WIDTH - n - 1);
		if (x < y)
			return -1;
		if (x > y)
			return 1;
	}

	return 0;
}

static void t_paulis_cmp_01(size_t tag)
{
	struct paulis a = paulis_new();
	struct paulis b = paulis_new();

	for (uint32_t n = 0; n < WIDTH; n++) {
		const int x = xoshiro256ss_next(&RNG) % 4;
		const int y = xoshiro256ss_next(&RNG) % 4;
		paulis_set(&a, x, n);
		paulis_set(&b, y, n);
	}

	const int res = paulis_cmp(a, b);
	const int exp = cmp_explicit(a, b);
	TEST_ASSERT(res == exp, "[%zu] res=%d, exp=%d", tag, res, exp);
}

static void t_paulis_cmp_02_eq(size_t tag)
{
	struct paulis a = paulis_new();
	struct paulis b = paulis_new();

	for (uint32_t n = 0; n < WIDTH; n++) {
		const int x = xoshiro256ss_next(&RNG) % 4;
		paulis_set(&a, x, n);
		paulis_set(&b, x, n);
	}
	TEST_ASSERT(paulis_eq(a, b), "[%zu]", tag);

	const int res = paulis_cmp(a, b);
	const int exp = cmp_explicit(a, b);
	TEST_ASSERT(res == exp, "[%zu] res=%d, exp=%d", tag, res, exp);
}

/*
 * Top-bit boundary: paulis_set / paulis_get and the
 * shift helpers must work at qubit 63 (the highest
 * bit of the 64-bit word).  Off-by-one in the mask
 * (e.g. `<< n` when n == 63 vs n == 64) would
 * silently wrap.
 */
static void t_paulis_n63(void)
{
	const enum pauli_op ops[] = { PAULI_I, PAULI_X, PAULI_Y, PAULI_Z };

	for (size_t k = 0; k < 4; k++) {
		struct paulis ps = paulis_new();
		paulis_set(&ps, ops[k], 63);
		TEST_EQ(paulis_get(ps, 63), ops[k]);

		/* Every other qubit must stay I. */
		for (uint32_t n = 0; n < 63; n++)
			TEST_EQ(paulis_get(ps, n), PAULI_I);
	}

	/* shl/shr across the top bit: a Pauli at qubit
	 * 0 shifted left by 63 lands at qubit 63 and
	 * round-trips on shr. */
	struct paulis ps = paulis_new();
	paulis_set(&ps, PAULI_Y, 0);
	paulis_shl(&ps, 63);
	TEST_EQ(paulis_get(ps, 63), PAULI_Y);
	TEST_EQ(paulis_get(ps, 0), PAULI_I);
	paulis_shr(&ps, 63);
	TEST_EQ(paulis_get(ps, 0), PAULI_Y);
	TEST_EQ(paulis_get(ps, 63), PAULI_I);
}

/*
 * paulis_cmp transitivity: pick three random Pauli
 * strings, sort them via cmp, verify the implied
 * order a <= b <= c yields cmp(a, c) <= 0 and
 * cmp(c, a) >= 0.  Run N iterations.
 */
static void t_paulis_cmp_transitive(size_t tag)
{
	struct paulis p[3];
	for (size_t i = 0; i < 3; i++) {
		p[i] = paulis_new();
		for (size_t k = 0; k < WIDTH; k++)
			paulis_set(&p[i],
				(enum pauli_op)(xoshiro256ss_next(&RNG) % 4),
				k);
	}

	/* Tiny three-element sort by cmp. */
	if (paulis_cmp(p[0], p[1]) > 0) {
		struct paulis t = p[0];
		p[0] = p[1];
		p[1] = t;
	}
	if (paulis_cmp(p[1], p[2]) > 0) {
		struct paulis t = p[1];
		p[1] = p[2];
		p[2] = t;
	}
	if (paulis_cmp(p[0], p[1]) > 0) {
		struct paulis t = p[0];
		p[0] = p[1];
		p[1] = t;
	}

	/* Pairwise: with the sort done, cmp(lo, hi) <= 0
	 * for every (i < j), and antisymmetric reverse. */
	for (size_t i = 0; i < 3; i++)
		for (size_t j = i + 1; j < 3; j++) {
			const int ij = paulis_cmp(p[i], p[j]);
			const int ji = paulis_cmp(p[j], p[i]);
			TEST_ASSERT(ij <= 0,
				"[%zu] cmp(p[%zu], p[%zu]) = %d "
				"after sort", tag, i, j, ij);
			TEST_ASSERT(ji >= 0,
				"[%zu] cmp(p[%zu], p[%zu]) = %d "
				"after sort", tag, j, i, ji);
		}
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	xoshiro256ss_init(&RNG, SEED);

	t_paulis_new();
	for (size_t n = 0; n < 999; n++) {
		t_paulis_eq(n);
		t_paulis_getset(n);
	}

	for (size_t n = 0; n < WIDTH; n++)
		t_paulis_shr(n);

	t_paulis_effect_00();
	t_paulis_effect_i1();
	t_paulis_effect_i2();
	t_paulis_effect_01();
	for (size_t n = 0; n < 999; n++)
		t_paulis_effect_02(n);

	for (size_t n = 0; n < 999; n++)
		t_paulis_split_01(n);

	t_paulis_cmp_00();
	for (size_t n = 0; n < 999; n++)
		t_paulis_cmp_01(n);
	for (size_t n = 0; n < 999; n++)
		t_paulis_cmp_02_eq(n);
	for (size_t n = 0; n < 999; n++)
		t_paulis_cmp_transitive(n);

	t_paulis_n63();

	world_free();
}
