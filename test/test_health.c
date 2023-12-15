/*
 * Test if we check initialize MPI environment and allocate memory.
 *
 */

#include <stdio.h>
#include <string.h>

#include "circ.h"

#include "test.h"

#define HEALTHCHECK_LABEL "healthcheck"

static struct circuit HC_CIRCUIT = { .name = HEALTHCHECK_LABEL,
				     .num_mea_qb = 7,
				     .num_sys_qb = 8,
				     .num_anc_qb = 9,
				     .reset = NULL,
				     .prepst = NULL,
				     .effect = NULL,
				     .measure = NULL };

TEST(healthcheck, void)
{
	struct circ *c = circ_create(&HC_CIRCUIT, NULL);
	TEST_ASSERT(c != NULL, "cannot initialize circuit");

	TEST_FINALIZE
	circ_destroy(c);
	circ_shutdown();
}
TEST_END

int main(void)
{
	return healthcheck();
}
