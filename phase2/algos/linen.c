/** circ: linen.
 *
 * Simple circ implementation for testing.
 */

#include "circ.h"
#include "algos/linen.h"

int linen_reset(struct circ *c)
{
	(void)c;

	return 0;
}

int linen_prepst(struct circ *c)
{
	struct linen_circ_data *dat = circ_data(c);
	dat->pass_prepst = dat->val_prepst;

	return 0;
}

int linen_effect(struct circ *c)
{
	struct linen_circ_data *dat = circ_data(c);
	dat->pass_effect = dat->val_effect;

	return 0;
}

int linen_measure(struct circ *c)
{
	struct linen_circ_data *dat = circ_data(c);
	dat->pass_measure = dat->val_measure;

	return 0;
}

void linen_circuit_init(struct circuit *ct, struct linen_circ_data *data)
{
	ct->name = LINEN_NAME;
	ct->num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB;
	ct->num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB;
	ct->num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB;
	ct->reset = linen_reset;
	ct->prepst = linen_prepst;
	ct->effect = linen_effect;
	ct->measure = linen_measure;
}

int linen_simulate(void)
{
	struct circuit ct;
	struct linen_circ_data cdat = { .val_prepst = 1,
					.val_effect = 22,
					.val_measure = 333 };
	linen_circuit_init(&ct, &cdat);

	struct linen_circ_data circ_dat;
	struct circ *c = circ_create(&ct, &circ_dat);
	if (!c)
		return -1;
	circ_reset(c);
	circ_report(c);
	circ_simulate(c);
	circ_destroy(c);

	return 0;
}
