/** Circuit interface.
 */

#ifndef CIRC_H
#define CIRC_H

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

typedef struct circuit_data_ {
    PauliHamil hamil;
    void* data;
} circuit_data;

typedef struct circ_sample_ {
    double time;
    int outcome;
    void* data;
} circ_sample;

/** circuit specification.
 */
typedef struct circuit_ {
    const char *name;
    circuit_data data;

    size_t num_mea_cl;
    size_t num_mea_qb;
    size_t num_sys_qb;
    size_t num_anc_qb;

    /** Reset the circuit.
     *
     * If the pointer is NULL, this will only zero the mea_cl bits.
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
    circ_result (*state_prep)(circ *, circ_sample*);

    circ_result (*routine)(circ *, circ_sample*);

    circ_result (*state_post)(circ *, circ_sample*);

    circ_result (*measure)(circ *, circ_sample*);
} circuit;

circ_env *circ_create_env();

void circ_destroy_env(circ_env *);

void circ_report_env(circ_env *);

circ *circ_create(circuit, circ_env *, circ_sample*);

void circ_destroy(circ *);

int *circ_mea_cl(circ *);

int *circ_mea_qb(circ *);

int *circ_sys_qb(circ *);

int *circ_anc_qb(circ *);

size_t circ_num_mea_cl(circ *);

size_t circ_num_mea_qb(circ *);

size_t circ_num_sys_qb(circ *);

size_t circ_num_anc_qb(circ *);

size_t circ_num_tot_qb(circ *);

const char *circ_name(circ *);

void circ_report(circ *);

circuit_data circ_circuit_data(circ *);

Qureg circ_qureg(circ *);

circ_result circ_reset(circ *);

circ_result circ_simulate(circ *);

#endif // CIRC_H