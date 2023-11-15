/** Circuit interface.
 */

#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "QuEST.h"

typedef struct circ_env_ circ_env;
typedef struct circuit_ circuit;
typedef struct circ_ circ;

typedef enum {
    CIRC_OK,
    CIRC_ERR,
} circ_result;

/** circuit environment for a run on MPI cluster.
 *
 * To be created and destroyed _only once_
 * during the experiment.
 */

struct circ_env_ {
    QuESTEnv quest_env;
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
    circ_result (*reset)(circ);

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
    circ_result (*state_prep)(circ, void *);

    circ_result (*routine)(circ, void *);

    circ_result (*state_post)(circ, void *);
};


/** circuit instance.
 */
struct circ_ {
    circuit ct;
    /** Arbitrary data passed during initialization.
     *
     * The user is responsible to manage the memory pointed to.
     */
    void *data;

    circ_env env;
    Qureg qureg;

    /** Subregisters:
     *
     *    mea (measurement)
     *       mea_cl - classical register
     *       mea_qb - quantum register
     *    sys (system)
     *    anc (ancilla)
     */
    int *mea_cl;
    double *mea_cl_prob;

    int *mea_qb;
    int *sys_qb;
    int *anc_qb;

    size_t simul_counter;
};

circ_result circ_env_init(circ_env *);

void circ_env_drop(circ_env);

void circ_env_report(circ_env);

circ_result circ_init(circ *, circuit, circ_env, void *);

void circ_drop(circ);

void circ_report(circ);

circ_result circ_reset(circ);

circ_result circ_simulate(circ);

#endif //PHASE2_CIRC_H