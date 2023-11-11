#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "../src/circ.h"

#include "mock_circ.h"


struct test_circ_fixture {
    circuit *ct;
    circ_env *env;
    void *data;
};

void case_qb(void *data, circ *c) {
    (void) data;
    size_t num_tot_qb = circ_num_tot_qb(c);
    int *idx_qb = malloc(sizeof(int) * num_tot_qb);
    assert(idx_qb);
    for (size_t i = 0; i < num_tot_qb; i++) {
        idx_qb[i] = i;
    }

    assert(
            memcmp(circ_mea_qb(c),
                   idx_qb,
                   circ_num_mea_qb(c)
            ) == 0
    );
    assert(
            memcmp(circ_sys_qb(c),
                   idx_qb + circ_num_mea_qb(c),
                   circ_num_sys_qb(c)
            ) == 0
    );
    assert(
            memcmp(circ_anc_qb(c),
                   idx_qb + circ_num_mea_qb(c) + circ_num_sys_qb(c),
                   circ_num_anc_qb(c)
            ) == 0);

    free(idx_qb);
}

void case_num_qb(void *data, circ *c) {
    struct test_circ_fixture *fix = (struct test_circ_fixture *) data;
    circuit *ct = fix->ct;

    assert(circ_num_mea_cl(c) == ct->num_mea_cl);
    assert(circ_num_mea_qb(c) == ct->num_mea_qb);
    assert(circ_num_sys_qb(c) == ct->num_sys_qb);
    assert(circ_num_anc_qb(c) == ct->num_anc_qb);
    assert(circ_num_tot_qb(c) ==
           ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb);
}

void case_name(void *data, circ *c) {
    struct test_circ_fixture *fix = (struct test_circ_fixture *) data;
    circuit *ct = fix->ct;
    assert(strcmp(circ_name(c), ct->name) == 0);
}

void case_simulate(void *data, circ *c) {
    struct test_circ_fixture *fix = (struct test_circ_fixture *) data;
    assert(circ_simulate(c) == CIRC_OK);
    struct mock_circ_sample *s = (struct mock_circ_sample *) fix->data;
    assert(s->result == s->input + 1);
}

static void test_circ_run(void(*test_case)(void *data, circ *), void *data) {
    struct test_circ_fixture *fixture = (struct test_circ_fixture *) data;
    circuit ct = *fixture->ct;
    circ_env *env = fixture->env;
    void *circ_data = fixture->data;

    circ *c = circ_create(ct, env, circ_data);
    assert(c != NULL);

    test_case(data, c);

    circ_destroy(c);
}

void test_circ_mock_circ_inst(circuit ct, circ_env *env) {
    struct test_circ_fixture fix = {.ct = &ct, .env = env, .data = NULL};
    test_circ_run(case_qb, &fix);
    test_circ_run(case_num_qb, &fix);
    test_circ_run(case_name, &fix);

    struct mock_circ_sample s = {.input = 0};
    fix.data = &s;
    test_circ_run(case_simulate, &fix);
    s.input = 1;
    test_circ_run(case_simulate, &fix);
    s.input = 2;
    test_circ_run(case_simulate, &fix);
}

void test_circ_mock_circ(circ_env *env, int const *num_qb, char const *name) {
    circuit ct = mock_circ_circuit;
    ct.num_mea_cl = num_qb[0];
    ct.num_mea_qb = num_qb[1];
    ct.num_sys_qb = num_qb[2];
    ct.num_anc_qb = num_qb[3];
    ct.name = name;
    test_circ_mock_circ_inst(ct, env);
}

void test_circ_unit01(circ_env *env) {
    test_circ_mock_circ(env, (int[]) {0, 1, 2, 3}, "mock1");

    test_circ_mock_circ(env, (int[]) {3, 1, 0, 0}, "mock2");
    test_circ_mock_circ(env, (int[]) {3, 0, 1, 0}, "mock2");
    test_circ_mock_circ(env, (int[]) {3, 0, 0, 1}, "mock2");

    test_circ_mock_circ(env, (int[]) {5, 5, 5, 5}, "mock3");
    test_circ_mock_circ(env, (int[]) {7, 7, 7, 7}, "mock3");
}

int main(int argc, char **argv) {
    (void) argc;
    (void) argv;

    circ_env *env = circ_create_env();

    test_circ_unit01(env);

    circ_destroy_env(env);
    return EXIT_SUCCESS;
}
