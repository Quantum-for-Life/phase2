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
	pauli_op_t op[WIDTH];
	struct paulis ps = paulis_new();

	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, op[k] = (int)(xoshiro256ss_next(&RNG) % 4), k);

	for (size_t k = 0; k < WIDTH; k++) {
		pauli_op_t p = paulis_get(ps, k);
		TEST_ASSERT(p == op[k],
			"[%zu] k=%zu, pauli=%d, expected=%d", tag, k, p, op[k]);
	}
}

static void t_paulis_eq(size_t tag)
{
	pauli_op_t op;
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


static void t_paulis_shr(size_t n)
{
	pauli_op_t op;
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

static void t_paulis_effect_00(void)
{
	_Complex double z = 11.3;
	struct paulis ps = paulis_new();

	for (size_t i = 0; i < 99; i++) {
		uint64_t x, y = xoshiro256ss_next(&RNG);
		x = paulis_effect(ps, y, &z);
		TEST_ASSERT(x == y,
			"[%zu] x=%lu, y=%lu", i, x, y);
		TEST_ASSERT(z == 11.3, "[%zu] z = %f + %fi",
			i, creal(z), cimag(z));
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

	TEST_EQ(paulis_effect(ps, 0x00, (void *)0), 0x0B); /* 0xB = 0b1011 */
	TEST_EQ(paulis_effect(ps, 0x01, (void *)0), 0x0A);
	TEST_EQ(paulis_effect(ps, 0x02, (void *)0), 0x09);
	TEST_EQ(paulis_effect(ps, 0x03, (void *)0), 0x08);


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
	pauli_op_t op;
	struct paulis ps = paulis_new();

	x = xoshiro256ss_next(&RNG);
	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, (int)(xoshiro256ss_next(&RNG) % 4), k);

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
	TEST_ASSERT(z == z_exp, "[%zu] z=%f+%fi, z_exp=%f+%fi", tag,
		creal(z), cimag(z), creal(z_exp), cimag(z_exp));

}

static void t_paulis_split_01(size_t tag)
{
	uint32_t lo, hi;
	struct paulis ps_lo, ps_hi, ps = paulis_new();
	pauli_op_t op_lo, op_hi, op_ex;

	lo = xoshiro256ss_next(&RNG) % (WIDTH - 1);
	hi = xoshiro256ss_next(&RNG) % (WIDTH - lo);
	TEST_ASSERT(lo + hi <= WIDTH, "wrong test params");

	for (size_t k = 0; k < WIDTH; k++)
		paulis_set(&ps, (int)(xoshiro256ss_next(&RNG) % 4), k);

	paulis_split(ps, lo, hi, &ps_lo, &ps_hi);

	for (uint32_t k = 0; k < lo; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);
		op_ex = paulis_get(ps, k);

		TEST_ASSERT(op_lo == op_ex,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d, op_ex=%d",
			tag, lo, hi, k, op_lo, op_ex);
		TEST_ASSERT(op_hi == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_hi=%d",
			tag, lo, hi, k, op_hi);
	}
	for (uint32_t k = lo; k < hi; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);
		op_ex = paulis_get(ps, k);

		TEST_ASSERT(op_lo == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d",
			tag, lo, hi, k, op_lo);
		TEST_ASSERT(op_hi == op_ex,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d, op_ex=%d",
			tag, lo, hi, k, op_hi, op_ex);
	}
	for (uint32_t k = lo + hi; k < WIDTH; k++) {
		op_lo = paulis_get(ps_lo, k);
		op_hi = paulis_get(ps_hi, k);

		TEST_ASSERT(op_lo == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_lo=%d",
			tag, lo, hi, k, op_lo);
		TEST_ASSERT(op_hi == PAULI_I,
			"[%zu], lo=%u, hi=%u, k=%u, op_hi=%d",
			tag, lo, hi, k, op_hi);
	}
}

static void TEST_MAIN(void)
{
	world_init((void *)0, (void *)0, WD_SEED);

	xoshiro256ss_init(&RNG, SEED);

	t_paulis_new();
	for (size_t n = 0; n < 9999; n++) {
		t_paulis_eq(n);
		t_paulis_getset(n);
	}
	
	for (size_t n = 0; n < WIDTH; n++)
		t_paulis_shr(n);

	t_paulis_effect_00();
	t_paulis_effect_i1();
	t_paulis_effect_i2();
	t_paulis_effect_01();
	for (size_t n = 0; n < 9999; n++)
		t_paulis_effect_02(n);

	for (size_t n = 0; n < 9999; n++)
		t_paulis_split_01(n);

	world_destroy();
}
