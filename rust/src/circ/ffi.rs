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
    pub(crate) name: *const c_char,
    pub(crate) data: *mut c_void,

    pub(crate) num_mea_qb: usize,
    pub(crate) num_sys_qb: usize,
    pub(crate) num_anc_qb: usize,

    reset: extern "C" fn(c: circ) -> circ_result,

    state_prep: extern "C" fn(c: circ, data: *mut c_void) -> circ_result,
    routine:    extern "C" fn(c: circ, data: *mut c_void) -> circ_result,
    state_post: extern "C" fn(c: circ, data: *mut c_void) -> circ_result,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct circ {
    ct:   circuit,
    data: *mut c_void,

    env:   circ_env,
    qureg: *mut Qureg,

    mea_cl:      *mut c_int,
    mea_cl_prob: *mut c_double,

    mea_qb: *mut c_int,
    sys_qb: *mut c_int,
    anc_qb: *mut c_int,
}

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

extern "C" {
    pub static linen_circuit: circuit;
}
