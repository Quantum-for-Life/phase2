#![feature(c_size_t)]
extern crate core;

use core::ffi::c_size_t;
use std::{
    ffi::{
        c_char,
        c_int,
        c_void,
        CStr,
        CString,
    },
    ptr::null_mut,
};

#[allow(non_camel_case_types)]
type circ_env = c_void;

#[allow(non_camel_case_types)]
type circ = c_void;

#[repr(C)]
#[allow(non_camel_case_types)]
enum circ_result {
    CIRC_OK,
    CIRC_ERR,
}

#[repr(C)]
#[allow(non_camel_case_types)]
struct circuit {
    name: *const c_char,
    data: *mut c_void,

    num_mea_cl: c_size_t,
    num_mea_qb: c_size_t,
    num_sys_qb: c_size_t,
    num_anc_qb: c_size_t,

    reset: extern "C" fn(*mut circ) -> circ_result,

    state_prep: extern "C" fn(*mut circ, *mut c_void) -> circ_result,
    routine:    extern "C" fn(*mut circ, *mut c_void) -> circ_result,
    state_post: extern "C" fn(*mut circ, *mut c_void) -> circ_result,
    measure:    extern "C" fn(*mut circ, *mut c_void) -> circ_result,
}

// #[link(name = "QuEST")]
#[link(name = "circ")]
#[allow(non_snake_case)]
extern "C" {
    fn circ_create_env() -> *mut circ_env;
    fn circ_destroy_env(env: *mut circ_env);
    fn circ_report_env(env: *mut circ_env);

    fn circ_create(
        ct: circuit,
        env: *mut circ_env,
        data: *mut c_void,
    ) -> *mut circ;
    fn circ_destroy(c: *mut circ);

    fn circ_mea_cl(c: *mut circ) -> *mut c_int;
    fn circ_mea_qb(c: *mut circ) -> *mut c_int;
    fn circ_sys_qb(c: *mut circ) -> *mut c_int;
    fn circ_anc_qb(c: *mut circ) -> *mut c_int;

    fn circ_num_mea_cl(c: *mut circ) -> c_size_t;
    fn circ_num_mea_qb(c: *mut circ) -> c_size_t;
    fn circ_num_sys_qb(c: *mut circ) -> c_size_t;
    fn circ_num_anc_qb(c: *mut circ) -> c_size_t;
    fn circ_num_tot_qb(c: *mut circ) -> c_size_t;

    fn circ_name(c: *mut circ) -> *const c_char;

    fn circ_report(c: *mut circ);

    fn circ_circuit_data(c: *mut circ) -> *mut c_void;

    // Qureg circ_qureg(circ *);

    fn circ_reset(c: *mut circ) -> circ_result;

    fn circ_simulate(c: *mut circ) -> circ_result;
}

extern "C" fn ttt_reset(c: *mut circ) -> circ_result {
    circ_result::CIRC_OK
}

extern "C" fn ttt_state_prep(
    c: *mut circ,
    data: *mut c_void,
) -> circ_result {
    circ_result::CIRC_OK
}

extern "C" fn ttt_routine(
    c: *mut circ,
    data: *mut c_void,
) -> circ_result {
    circ_result::CIRC_OK
}

extern "C" fn ttt_state_post(
    c: *mut circ,
    data: *mut c_void,
) -> circ_result {
    circ_result::CIRC_OK
}

extern "C" fn ttt_measure(
    c: *mut circ,
    data: *mut c_void,
) -> circ_result {
    circ_result::CIRC_OK
}

fn main() {
    println!("Hello, world!");

    let name = CString::new("ttt").unwrap();
    let ct = circuit {
        name:       name.as_ptr(),
        data:       null_mut(),
        num_mea_cl: 3,
        num_mea_qb: 3,
        num_sys_qb: 3,
        num_anc_qb: 0,
        reset:      ttt_reset,
        state_prep: ttt_state_prep,
        routine:    ttt_routine,
        state_post: ttt_state_post,
        measure:    ttt_measure,
    };

    unsafe {
        let env = circ_create_env();
        circ_report_env(env);

        let c = circ_create(ct, env, null_mut());
        circ_report(c);

        let slice = CStr::from_ptr(circ_name(c));
        println!("string returned: {}", slice.to_str().unwrap());

        circ_simulate(c);

        circ_destroy(c);
        circ_destroy_env(env);
    }
}
