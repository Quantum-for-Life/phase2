#include <complex.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "paulis.h"

static int RET = 0;

#define TEST_FAIL(...)                                                         \
	{                                                                      \
		fprintf(stderr, "FAILED %s:%d \"", __FILE__, __LINE__);        \
		fprintf(stderr, __VA_ARGS__);                                  \
		fprintf(stderr, "\"\n");                                       \
		RET = -1;                                                      \
	}

#define NUM_PAULI_STRS (22)
static char *PAULI_STR[NUM_PAULI_STRS] = {
	"I",
	"X",
	"Y",
	"Z",
	"XX",
	"YX",
	"ZX",
	"XY",
	"YY",
	"ZY",
	"XZ",
	"YZ",
	"ZZ",
	"XYXXZIYYYZZY",
	"IZZIZZIZ",
	"YZIXXXXXXXXXXZ",
	"IIIIIIIIZIIIIIIIX",
	"YYYYZYZZYZZZZZIIIZZ",
	"XXYYXXYXXXY",
	"ZXYXZYXZYZXZYZXIIZX",
	"ZZXY",
	"XYIYXYZXYXZYXYXYXXIXZXZ",
};

static struct paulis paulis_fromstr(const char *s)
{
	struct paulis code = paulis_new();

	int i = 0;
	while (s[i] != '\0' && i < PAULI_MAX_WIDTH) {
		enum pauli_op op;

		switch (s[i]) {
		case 'I':
			op = PAULI_I;
			break;
		case 'X':
			op = PAULI_X;
			break;
		case 'Y':
			op = PAULI_Y;
			break;
		case 'Z':
			op = PAULI_Z;
			break;
		default:
			TEST_FAIL("ERROR: %s[%d] = %c", s, i, s[i]);
			exit(EXIT_FAILURE);
		}
		paulis_set(&code, op, i);

		i++;
	}

	return code;
}

void test_paulis_getset(int i)
{
	const char	   *pstr = PAULI_STR[i];
	const struct paulis code = paulis_fromstr(pstr);

	int n = 0;
	while (pstr[n] != '\0') {
		const enum pauli_op have = paulis_get(code, n);

		enum pauli_op should;
		switch (pstr[n]) {
		case 'I':
			should = PAULI_I;
			break;
		case 'X':
			should = PAULI_X;
			break;
		case 'Y':
			should = PAULI_Y;
			break;
		case 'Z':
			should = PAULI_Z;
			break;
		default:
			TEST_FAIL(
				"ERROR: PAULI_STR[%d][%d] = %c", i, n, pstr[n]);
			exit(EXIT_FAILURE);
		}

		if (have != should) {
			TEST_FAIL("pauli op = %d, but PAULI_STR[%d][%d] = %d",
				have, i, n, should);
			break;
		}
		n++;
	}

	for (; n < PAULI_MAX_WIDTH; n++) {
		const enum pauli_op have = paulis_get(code, n);

		if (have != PAULI_I) {
			TEST_FAIL("remaining pauli_op = %d, but should be 0",
				have);
			break;
		}
	}
}

void test_paulis_eff(const int str_idx)
{
	const char	   *pstr = PAULI_STR[str_idx];
	const struct paulis code = paulis_fromstr(pstr);

	u32 len = 0;
	while (pstr[len] != '\0')
		len++;

	const u64 num_amps = (u64)1 << len;
	for (u64 i = 0; i < num_amps; i++) {
		u64 j_should  = i;
		f64 z_have[2] = { 1.0, 0.0 }, z_should[2] = { 1.0, 0.0 };

		const u64 j_have = paulis_effect(code, i, &z_have);

		for (u32 k = 0; k < len; k++) {
			const enum pauli_op pauli = paulis_get(code, k);

			const u64 k_mask = (u64)1 << k;
			switch (pauli) {
			case PAULI_I:
				break;
			case PAULI_X:
				j_should ^= k_mask;
				break;
			case PAULI_Y:
				j_should ^= k_mask;
				const fl a  = z_should[0];
				z_should[0] = -z_should[1];
				z_should[1] = a;
				if ((i & k_mask) != 0) {
					z_should[0] *= -1;
					z_should[1] *= -1;
				}
				break;
			case PAULI_Z:
				if ((i & k_mask) != 0) {
					z_should[0] *= -1;
					z_should[1] *= -1;
				}
				break;
			default:
				break;
			}
		}

		if (j_have != j_should ||
			abs(z_have[0] - z_should[0]) > DBL_EPSILON ||
			abs(z_have[1] - z_should[1]) > DBL_EPSILON) {
			TEST_FAIL(
				"%s |%lu> = (%f+%fi) |%lu>, but should be (%f+%fi) |%lu>",
				pstr, i, z_have[0], z_have[1], j_have,
				z_should[0], z_should[1], j_should);
			break;
		}
	}
}

int main(void)
{
	for (int i = 0; i < NUM_PAULI_STRS; i++)
		test_paulis_getset(i);

	for (int i = 0; i < NUM_PAULI_STRS; i++)
		test_paulis_eff(i);

	return RET;
}
