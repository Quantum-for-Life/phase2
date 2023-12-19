#include <complex.h>

#include "QuEST.h"

#include "circ.h"
#include "circ_private.h"

void
capi_hadamard(struct circ *c, qbid qb)
{
	Qureg qureg = circ_intl_quest_qureg(c);
	hadamard(qureg, qb);
}

void
capi_sgate(struct circ *c, qbid qb)
{
	Qureg qureg = circ_intl_quest_qureg(c);
	sGate(qureg, qb);
}

double
capi_prob0(struct circ *c, qbid qb)
{
	Qureg qureg = circ_intl_quest_qureg(c);
	return calcProbOfOutcome(qureg, qb, 0);
}

void
capi_blankstate(struct circ *c)
{
	Qureg qureg = circ_intl_quest_qureg(c);
	initBlankState(qureg);
}

void
capi_set_sysamp(struct circ *c, size_t idx, _Complex double amp)
{
	Qureg	  qureg	    = circ_intl_quest_qureg(c);
	double	  amp_re    = creal(amp);
	double	  amp_im    = cimag(amp);
	long long start_ind = idx << circ_num_meaqb(c);
	setAmps(qureg, start_ind, &amp_re, &amp_im, 1);
}

void
capi_ctl_rotate_pauli(struct circ *c, int *paulis, double angle)
{
	Qureg	     qureg	= circ_intl_quest_qureg(c);
	int	    *qb		= circ_intl_get_qb(c);
	const size_t num_mea_qb = circ_num_meaqb(c);
	const size_t num_sys_qb = circ_num_sysqb(c);
	for (size_t i = 0; i < num_mea_qb + num_sys_qb; i++) {
		qb[i] = i;
	}
	int *qb_ctl = qb;
	int *qb_trg = qb + num_mea_qb;

	multiControlledMultiRotatePauli(qureg, qb_ctl, num_mea_qb, qb_trg,
		(enum pauliOpType *)paulis, num_sys_qb, -2.0 * angle);
}