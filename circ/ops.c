#include <complex.h>

#include "qreg.h"

#include "circ.h"
#include "circ_private.h"

void circ_ops_hadamard(struct circ *c, qbid qb)
{
	//hadamard(c->quest_qureg, qb);
}

void circ_ops_sgate(struct circ *c, qbid qb)
{
	//sGate(c->quest_qureg, qb);
}


void circ_ops_sgate_conj(struct circ *c, qbid qb)
{
	//sGate(c->quest_qureg, qb);
	//pauliZ(c->quest_qureg, qb);
}

double circ_ops_prob0(struct circ *c, qbid qb)
{
	return 0.0; //calcProbOfOutcome(c->quest_qureg, qb, 0);
}

void circ_ops_blankstate(struct circ *c)
{
	qreg_zerostate(&c->quest_qureg);
}

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	double	  amp_re    = creal(amp);
	double	  amp_im    = cimag(amp);
	long long start_ind = idx << circ_num_meaqb(c);
	//setAmps(c->quest_qureg, start_ind, &amp_re, &amp_im, 1);
}

void circ_ops_crotpauli(struct circ *c, int *paulis, double angle)
{
	const size_t num_mea_qb = circ_num_meaqb(c);
	const size_t num_sys_qb = circ_num_sysqb(c);
	for (size_t i = 0; i < num_mea_qb + num_sys_qb; i++) {
		c->qb[i] = i;
	}
	int *qb_ctl = c->qb;
	int *qb_trg = c->qb + num_mea_qb;

	//multiControlledMultiRotatePauli(c->quest_qureg, qb_ctl, num_mea_qb,
	//	qb_trg, (enum pauliOpType *)paulis, num_sys_qb, -2.0 * angle);
}
