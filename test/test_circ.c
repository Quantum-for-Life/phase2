#include <stdarg.h>
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

int mock_reset(const struct circ *c)
{
	struct mock_circ_data *dat = c->data;
	dat->reset_val = 777;

	return 0;
}

int mock_state_prep(const struct circ *c)
{
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 222;

	return 0;
}

int mock_routine(const struct circ *c)
{
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 333;

	return 0;
}

int mock_state_post(const struct circ *c)
{
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 444;

	return 0;
}

static struct circuit MOCK_CIRCUIT = { .name = MOCK_CIRCUIT_NAME,
				       .data = NULL,
				       .num_mea_qb = 1,
				       .num_sys_qb = 2,
				       .num_anc_qb = 3,
				       .reset = mock_reset,
				       .state_prep = mock_state_prep,
				       .routine = mock_routine,
				       .state_post = mock_state_post };

TEST(mock_circ_init, struct circ_env *env, struct circuit *ct)
{
	struct circ c;
	struct mock_circ_data dat;
	TEST_ASSERT(circ_init(&c, env, ct, &dat) == 0,
		    "cannot initialize circuit")
	TEST_ASSERT(memcmp(ct->name, MOCK_CIRCUIT_NAME, 4) == 0,
		    "wrong circuit passed")

	TEST_FINALIZE
	circ_destroy(&c);
}
TEST_END

TEST(mock_circ_reset, struct circ_env *env, struct circuit *ct)
{
	struct circ c;
	struct mock_circ_data dat = { .index = 0 };

	TEST_ASSERT(circ_init(&c, env, ct, &dat) == 0,
		    "cannot initialize circuit")
	TEST_ASSERT(circ_reset(&c) == 0, "reset")
	TEST_ASSERT(dat.reset_val == 777, "circuit not reset")

	TEST_FINALIZE
	circ_destroy(&c);
}
TEST_END

TEST(mock_circ_simulate, struct circ_env *env, struct circuit *ct)
{
	struct circ c;
	struct mock_circ_data dat = { .index = 0 };

	TEST_ASSERT(circ_init(&c, env, ct, &dat) == 0,
		    "cannot initialize circuit")
	TEST_ASSERT(circ_simulate(&c) == 0, "simulation error")
	TEST_ASSERT(dat.reset_val == 777, "reset value")
	TEST_ASSERT(dat.values[0] == 222, "state_prep value")
	TEST_ASSERT(dat.values[1] == 333, "routine value")
	TEST_ASSERT(dat.values[2] == 444, "state_post value")

	TEST_FINALIZE
	circ_destroy(&c);
}
TEST_END

TEST(mock_circ_suite, void)
{
	struct circ_env env;
	TEST_ASSERT(circ_env_init(&env) == 0,
		    "cannot initializa circ environment");

	TEST_CASE(mock_circ_init(&env, &MOCK_CIRCUIT))
	TEST_CASE(mock_circ_reset(&env, &MOCK_CIRCUIT))
	TEST_CASE(mock_circ_simulate(&env, &MOCK_CIRCUIT))

	TEST_FINALIZE
	circ_env_destroy(&env);
}
TEST_END

int main(void)
{
	return mock_circ_suite();
}
