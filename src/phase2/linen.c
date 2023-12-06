/** circ: linen.
 *
 * Simple circ implementation for testing.
 */

#include "circ.h"
#include "linen.h"

int linen_reset(struct circ *c)
{
	(void)c;

	return 0;
}

int linen_state_prep(struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->state_prep_value = ct_dat->state_prep_value;

	return 0;
}

int linen_state_post(struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->state_post_value = ct_dat->state_post_value;

	return 0;
}

int linen_routine(struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->routine_value = ct_dat->routine_value;

	return 0;
}

void linen_circuit_init(struct circuit *ct, struct linen_circuit_data *ct_dat)
{
	ct->name = LINEN_NAME;
	ct->data = (void *)ct_dat;
	ct->num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB;
	ct->num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB;
	ct->num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB;
	ct->reset = linen_reset;
	ct->state_prep = linen_state_prep;
	ct->routine = linen_routine;
	ct->state_post = linen_state_post;
}

int linen_simulate(struct circ_env env)
{
	circ_env_report(&env);

	struct circuit ct;
	struct linen_circuit_data ct_dat = { .state_prep_value = 1,
					     .routine_value = 22,
					     .state_post_value = 333 };
	linen_circuit_init(&ct, &ct_dat);

	struct circ c;
	struct linen_circ_data circ_dat;
	if (circ_init(&c, &env, &ct, &circ_dat) < 0)
		return -1;
	circ_report(&c);
	circ_simulate(&c);
	circ_destroy(&c);

	return 0;
}
