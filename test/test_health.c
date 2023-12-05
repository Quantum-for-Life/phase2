/*
 * Test if we check initialize MPI environment and allocate memory.
 *
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "circ.h"

#include "test.h"

#define HEALTHCHECK_LABEL "healthcheck"

static struct circuit HC_CIRCUIT = { .name = HEALTHCHECK_LABEL,
				     .data = NULL,
				     .num_mea_qb = 7,
				     .num_sys_qb = 8,
				     .num_anc_qb = 9,
				     .reset = NULL,
				     .state_prep = NULL,
				     .routine = NULL,
				     .state_post = NULL };

TEST(healthcheck, struct circ_env *env)
{
	struct circ c;
	TEST_ASSERT(circ_init(&c, &env, &HC_CIRCUIT, NULL) == 0,
		    "cannot initialize circuit");

	TEST_ASSERT(c.ct->num_mea_qb == 7, " ")
	TEST_ASSERT(c.ct->num_sys_qb == 8, " ")
	TEST_ASSERT(c.ct->num_anc_qb == 9, " ")

	TEST_FINALIZE
	circ_destroy(&c);
}
TEST_END

int main(void)
{
	int rc = 0;

	struct circ_env env;
	circ_env_init(&env);

	rc |= healthcheck(&env);

	circ_env_destroy(&env);
	return rc;
}
