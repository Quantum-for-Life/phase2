#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "circ.h"

#include "test.h"

#define MOCK_CIRCUIT_NAME "mock"

struct mock_circ_data {
	int index;
	int values[4];
};

int mock_reset(struct circ *c) {
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 111;

	return CIRC_OK;
}

int mock_state_prep(struct circ *c) {
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 222;

	return CIRC_OK;
}

int mock_routine(struct circ *c) {
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 333;

	return CIRC_OK;
}

int mock_state_post(struct circ *c) {
	struct mock_circ_data *dat = c->data;
	dat->values[dat->index++] = 444;

	return CIRC_OK;
}

static struct circuit MOCK_CIRCUIT = {
	.name = MOCK_CIRCUIT_NAME,
	.data = NULL,
	.num_mea_qb = 1,
	.num_sys_qb = 2,
	.num_anc_qb = 3,
	.reset = mock_reset,
	.state_prep = mock_state_prep,
	.routine = mock_routine,
	.state_post = mock_state_post
};

TEST(mock_circ_init, struct circ_env env, struct circuit ct)
	struct circ c;

	TEST_ASSERT(circ_init(&c, env, ct, NULL) == CIRC_OK,
	            "cannot initialize circuit")
	TEST_ASSERT(memcmp(ct.name, MOCK_CIRCUIT_NAME, 4) == 0,
	            "wrong circuit passed")

	TEST_FINALIZE
	circ_destroy(&c);
TEST_END


TEST(mock_circ_reset, struct circ_env env, struct circuit ct)
	struct circ c;
	struct mock_circ_data dat = {.index = 0};

	TEST_ASSERT(circ_init(&c, env, ct, &dat) == CIRC_OK,
	            "cannot initialize circuit")
	TEST_ASSERT(circ_reset(&c) == CIRC_OK, "reset")
	TEST_ASSERT(dat.values[0] == 111, "circuit not reset")

	TEST_FINALIZE
	circ_destroy(&c);
TEST_END


TEST(mock_circ_simulate, struct circ_env env, struct circuit ct)
	struct circ c;
	struct mock_circ_data dat = {.index = 0};

	TEST_ASSERT(circ_init(&c, env, ct, &dat) == CIRC_OK,
	            "cannot initialize circuit")
	TEST_ASSERT(circ_simulate(&c) == CIRC_OK, "simulation error")
	TEST_ASSERT(dat.values[0] == 111, "reset value")
	TEST_ASSERT(dat.values[1] == 222, "state_prep value")
	TEST_ASSERT(dat.values[2] == 333, "routine value")
	TEST_ASSERT(dat.values[3] == 444, "state_post value")

	TEST_FINALIZE
	circ_destroy(&c);
TEST_END

TEST(mock_circ_suite, void)
	struct circ_env env;
	TEST_ASSERT(circ_env_init(&env) == CIRC_OK,
	            "cannot initializa circ environment");

	TEST_CASE(mock_circ_init(env, MOCK_CIRCUIT))
	TEST_CASE(mock_circ_reset(env, MOCK_CIRCUIT))
	TEST_CASE(mock_circ_simulate(env, MOCK_CIRCUIT))

	TEST_FINALIZE
	circ_env_destroy(&env);
TEST_END


int main(void) {
	return mock_circ_suite();
}
