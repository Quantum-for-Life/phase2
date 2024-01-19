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
	dat->index		   = 0;
	dat->reset_val		   = 777;

	return 0;
}

int mock_state_prep(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++]  = 222;

	return 0;
}

int mock_routine(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++]  = 333;

	return 0;
}

int mock_state_post(struct circ *c)
{
	struct mock_circ_data *dat = circ_data(c);
	dat->values[dat->index++]  = 444;

	return 0;
}

static struct circuit MOCK_CIRCUIT = { .name = MOCK_CIRCUIT_NAME,
	.num_mea_qb			     = 1,
	.num_sys_qb			     = 2,
	.num_anc_qb			     = 3,
	.reset				     = mock_reset,
	.prepst				     = mock_state_prep,
	.effect				     = mock_routine,
	.measure			     = mock_state_post };

static int mock_circ_init(struct circuit *ct)
{
	struct mock_circ_data dat;
	struct circ	    *c = circ_create(ct, &dat);

	if (c == NULL) {
		TEST_FAIL("cannot initialize circuit");
		goto error;
	}
	if (memcmp(ct->name, MOCK_CIRCUIT_NAME, 4) != 0) {
		TEST_FAIL("wrong circuit passed");
		goto error;
	}

	circ_destroy(c);
	return 0;
error:
	circ_destroy(c);
	return -1;
}

static int mock_circ_reset(struct circuit *ct)
{
	struct mock_circ_data dat = { .index = 0 };
	struct circ	    *c	  = circ_create(ct, &dat);

	if (c == NULL) {
		TEST_FAIL("cannot initialize circuit");
		goto error;
	}
	if (circ_reset(c) != 0) {
		TEST_FAIL("reset");
		goto error;
	}
	if (dat.reset_val != 777) {
		TEST_FAIL("circuit not reset");
		goto error;
	}

	circ_destroy(c);
	return 0;
error:
	circ_destroy(c);
	return -1;
}

static int mock_circ_simulate(struct circuit *ct)
{
	struct mock_circ_data dat = { .index = 0 };
	struct circ	    *c	  = circ_create(ct, &dat);

	if (c == NULL) {
		TEST_FAIL("cannot initialize circuit");
		goto error;
	}
	if (circ_run(c) != 0) {
		TEST_FAIL("simulation error");
		goto error;
	}
	if (dat.reset_val != 777) {
		TEST_FAIL("reset value");
		goto error;
	}
	if (dat.values[0] != 222) {
		TEST_FAIL("state_prep value");
		goto error;
	}
	if (dat.values[1] != 333) {
		TEST_FAIL("effect value");
		goto error;
	}
	if (dat.values[2] != 444) {
		TEST_FAIL("measure value");
		goto error;
	}

	circ_destroy(c);
	return 0;
error:
	circ_destroy(c);
	return -1;
}

int main(void)
{
	if (mock_circ_init(&MOCK_CIRCUIT) != 0)
		goto error;
	if (mock_circ_reset(&MOCK_CIRCUIT) != 0)
		goto error;
	if (mock_circ_simulate(&MOCK_CIRCUIT) != 0)
		goto error;
	if (mock_circ_simulate(&MOCK_CIRCUIT) != 0)
		goto error;

	return 0;
error:
	return -1;
}
