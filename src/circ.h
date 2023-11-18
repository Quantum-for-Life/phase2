/** Circuit interface.
 */

#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "QuEST.h"

typedef struct circ_env_ circ_env;
typedef struct circuit_ circuit;
typedef struct circ_ circ;

enum {
    circ_ok,
    circ_err,
};

/** circuit environment for a run on MPI cluster.
 *
 * To be created and destroyed _only once_
 * during the experiment.
 */

struct circ_env_ {
    QuESTEnv *quest_env;
};

/** circuit specification.
 */
struct circuit_ {
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
    int (*reset)(circ *);

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
    int (*state_prep)(circ *, void *);

    int (*routine)(circ *, void *);

    int (*state_post)(circ *, void *);
};


/** circuit instance.
 */
struct circ_ {
    circuit ct;
    void *data;

    circ_env env;
    Qureg *qureg;

    int *mea_cl;
    double *mea_cl_prob;

    int *mea_qb;
    int *sys_qb;
    int *anc_qb;
};

int circ_env_init(circ_env *);

void circ_env_destroy(circ_env *);

void circ_env_report(circ_env);

int circ_init(circ *, circuit, circ_env, void *);

void circ_destroy(circ *);

void circ_report(circ);

int circ_reset(circ *);

int circ_simulate(circ *);

#endif //PHASE2_CIRC_H