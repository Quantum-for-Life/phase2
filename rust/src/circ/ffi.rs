use std::ffi::{
    c_char,
    c_double,
    c_int,
    c_void,
};

type QuESTEnv = c_void;
type Qureg = c_void;

#[derive(Debug, Copy, Clone, PartialEq)]
#[repr(C)]
pub(crate) enum circ_result {
    CIRC_OK,
    CIRC_ERR,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct circ_env {
    quest_env: *mut QuESTEnv,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct circuit {
    name: *const c_char,

    pub(crate) data: *mut c_void,

    num_mea_qb: usize,
    num_sys_qb: usize,
    num_anc_qb: usize,

    reset: extern "C" fn(c: *mut circ) -> c_int,

    state_prep: extern "C" fn(c: *mut circ) -> c_int,
    routine:    extern "C" fn(c: *mut circ) -> c_int,
    state_post: extern "C" fn(c: *mut circ) -> c_int,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct circ {
    env:  circ_env,
    ct:   circuit,
    data: *mut c_void,

    qureg: *mut Qureg,

    mea_cl:      *mut c_int,
    mea_cl_prob: *mut c_double,
    mea_qb:      *mut c_int,
    sys_qb:      *mut c_int,
    anc_qb:      *mut c_int,
}

#[link(name = "phase2")]
extern "C" {

    pub(crate) fn circ_env_init(env: *mut circ_env) -> circ_result;

    pub(crate) fn circ_env_destroy(env: *mut circ_env);

    pub(crate) fn circ_env_report(env: *const circ_env);

    pub(crate) fn circ_init(
        c: *mut circ,
        env: circ_env,
        ct: circuit,
        data: *mut c_void,
    ) -> circ_result;

    pub(crate) fn circ_destroy(c: *mut circ);

    pub(crate) fn circ_report(c: *const circ);

    pub(crate) fn circ_reset(c: *mut circ) -> circ_result;

    pub(crate) fn circ_simulate(c: *mut circ) -> circ_result;

}

#[derive(Debug)]
#[repr(C)]
pub(crate) struct linen_circuit_data {
    pub(crate) state_prep_value: c_int,
    pub(crate) routine_value:    c_int,
    pub(crate) state_post_value: c_int,
}

#[derive(Debug)]
#[repr(C)]
pub(crate) struct linen_circ_data {
    pub(crate) state_prep_value: c_int,
    pub(crate) routine_value:    c_int,
    pub(crate) state_post_value: c_int,
}

#[link(name = "phase2")]
extern "C" {
    pub static linen_circuit: circuit;
}
