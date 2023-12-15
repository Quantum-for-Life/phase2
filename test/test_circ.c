/**
 * A mock-up circuit for testing.
 */

#include <stdio.h>
#include <string.h>

#include "circ.h"

#include "test.h"

#define MOCK_CIRCUIT_NAME "mock"

struct mock_circ_data {
	int reset_val;
	int index;
	int values[3];
};

int mock_reset(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->index = 0;
	dat->reset_val = 777;

	return 0;
}

int mock_state_prep(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++] = 222;

	return 0;
}

int mock_routine(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++] = 333;

	return 0;
}

int mock_state_post(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++] = 444;

	return 0;
}

static struct circuit MOCK_CIRCUIT = { .name = MOCK_CIRCUIT_NAME,
				       .num_mea_qb = 1,
				       .num_sys_qb = 2,
				       .num_anc_qb = 3,
				       .reset = mock_reset,
				       .prepst = mock_state_prep,
				       .effect = mock_routine,
				       .measure = mock_state_post };

TEST(mock_circ_init, struct circuit *ct)
{
	struct mock_circ_data dat;
	struct circ *c = circ_init(ct, &dat);
	TEST_ASSERT(c != NULL, "cannot initialize circuit")
	TEST_ASSERT(memcmp(ct->name, MOCK_CIRCUIT_NAME, 4) == 0,
		    "wrong circuit passed")

	TEST_FINALIZE
	circ_destroy(c);
}
TEST_END

TEST(mock_circ_reset, struct circuit *ct)
{
	struct mock_circ_data dat = { .index = 0 };
	struct circ *c = circ_init(ct, &dat);

	TEST_ASSERT(c != NULL, "cannot initialize circuit")
	TEST_ASSERT(circ_reset(c) == 0, "reset")
	TEST_ASSERT(dat.reset_val == 777, "circuit not reset")

	TEST_FINALIZE
	circ_destroy(c);
}
TEST_END

TEST(mock_circ_simulate, struct circuit *ct)
{
	struct mock_circ_data dat = { .index = 0 };
	struct circ *c = circ_init(ct, &dat);

	TEST_ASSERT(c != NULL, "cannot initialize circuit")
	TEST_ASSERT(circ_simulate(c) == 0, "simulation error")
	TEST_ASSERT(dat.reset_val == 777, "reset value")
	TEST_ASSERT(dat.values[0] == 222, "state_prep value")
	TEST_ASSERT(dat.values[1] == 333, "effect value")
	TEST_ASSERT(dat.values[2] == 444, "measure value")

	TEST_FINALIZE
	circ_destroy(c);
}
TEST_END

TEST(mock_circ_suite, void)
{
	TEST_CASE(mock_circ_init(&MOCK_CIRCUIT))
	TEST_CASE(mock_circ_reset(&MOCK_CIRCUIT))
	TEST_CASE(mock_circ_simulate(&MOCK_CIRCUIT))
	TEST_CASE(mock_circ_simulate(&MOCK_CIRCUIT))

	TEST_FINALIZE
	circ_env_shutdown();
}
TEST_END

int main(void)
{
	return mock_circ_suite();
}
