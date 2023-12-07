/** circ: linen.
 *
 * Simple circ implementation for testing.
 */

#include "circ.h"
#include "linen.h"

int linen_reset(const struct circ *c)
{
	(void)c;

	return 0;
}

int linen_prepst(const struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->pass_prepst = ct_dat->val_prepst;

	return 0;
}

int linen_effect(const struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->pass_effect = ct_dat->val_effect;

	return 0;
}

int linen_measure(const struct circ *c)
{
	const struct linen_circuit_data *ct_dat = c->ct->data;
	struct linen_circ_data *dat = c->data;
	dat->pass_measure = ct_dat->val_measure;

	return 0;
}

void linen_circuit_init(struct circuit *ct, struct linen_circuit_data *ct_dat)
{
	ct->name = LINEN_NAME;
	ct->data = ct_dat;
	ct->num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB;
	ct->num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB;
	ct->num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB;
	ct->reset = (circ_op)linen_reset;
	ct->prepst = (circ_op)linen_prepst;
	ct->effect = (circ_op)linen_effect;
	ct->measure = (circ_op)linen_measure;
}

int linen_simulate(struct circ_env *env)
{
	circ_env_report(env);

	struct circuit ct;
	struct linen_circuit_data ct_dat = { .val_prepst = 1,
					     .val_effect = 22,
					     .val_measure = 333 };
	linen_circuit_init(&ct, &ct_dat);

	struct circ c;
	struct linen_circ_data circ_dat;
	if (circ_init(&c, env, &ct, &circ_dat) < 0)
		return -1;
	circ_reset(&c);
	circ_report(&c);
	circ_simulate(&c);
	circ_destroy(&c);

	return 0;
}
