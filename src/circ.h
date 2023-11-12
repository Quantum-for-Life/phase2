/** Circuit interface.
 */

#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "QuEST.h"


typedef enum {
    CIRC_OK,
    CIRC_ERR,
} circ_result;

/** circuit environment for a run on MPI cluster.
 *
 * To be created and destroyed _only once_
 * during the experiment.
 */
typedef struct circ_env_ circ_env;

/** circuit instance.
 */
typedef struct circ_ circ;

/** circuit specification.
 */
typedef struct circuit_ {
    const char *name;
    void *data;

    size_t num_mea_qb;
    size_t num_sys_qb;
    size_t num_anc_qb;

    /** Reset the circuit.
     *
     * On top of standard resetting the Qureg and mea_cl.
     * This can be NULL.
     */
    circ_result (*reset)(circ *);

    /** circuit specification divided into 4 steps.
     *
     * The order of execution of these steps in circ_simulate() is always:
     *   1. state_prep()
     *   2. routine()
     *   3. state_post()
     *   4. measure()
     *
     * The cicuit is _not_ reset before the call to circ_simulate().
     *
     * The data pointer passed to these function is the same stored in
     * circ.data, and the same as passed to circ_create_circuit().
     *
     * There pointers can be NULL, meaning the step is not specified.
     */
    circ_result (*state_prep)(circ *, void *);

    circ_result (*routine)(circ *, void *);

    circ_result (*state_post)(circ *, void *);

} circuit;


circ_env *circ_create_env();

void circ_destroy_env(circ_env *);

void circ_report_env(circ_env *);

circ *circ_create(circuit, circ_env *, void *);

void circ_destroy(circ *);

int *circ_mea_cl(circ *);

double *circ_mea_cl_prob(circ *);

int *circ_mea_qb(circ *);

int *circ_sys_qb(circ *);

int *circ_anc_qb(circ *);

size_t circ_num_mea_qb(circ *);

size_t circ_num_sys_qb(circ *);

size_t circ_num_anc_qb(circ *);

size_t circ_num_tot_qb(circ *);

const char *circ_name(circ *);

void circ_report(circ *);

void *circ_circuit_data(circ *);

Qureg circ_qureg(circ *);

circ_result circ_reset(circ *);

circ_result circ_simulate(circ *);

#endif //PHASE2_CIRC_H