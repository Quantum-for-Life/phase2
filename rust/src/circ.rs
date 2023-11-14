#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::{
    mem,
    mem::MaybeUninit,
    ptr,
};

use crate::circ::ffi::{
    circ_env,
    circ_env_init,
    circ_report,
    circ_result,
};

pub(crate) mod ffi {
    use std::ffi::{
        c_char,
        c_double,
        c_int,
        c_void,
    };

    use crate::quest_sys;

    #[derive(Debug, Copy, Clone, PartialEq)]
    #[repr(C)]
    pub(crate) enum circ_result {
        CIRC_OK,
        CIRC_ERR,
    }

    #[derive(Debug, Copy, Clone)]
    #[repr(C)]
    pub(crate) struct circ_env {
        quest_env: quest_sys::QuESTEnv,
    }

    #[derive(Debug, Copy, Clone)]
    #[repr(C)]
    pub(crate) struct circuit {
        pub(crate) name: *const c_char,
        pub(crate) data: *mut c_void,

        pub(crate) num_mea_qb: usize,
        pub(crate) num_sys_qb: usize,
        pub(crate) num_anc_qb: usize,

        reset: extern "C" fn(c: circ),

        state_prep: extern "C" fn(c: circ, data: *mut c_void),
        routine:    extern "C" fn(c: circ, data: *mut c_void),
        state_post: extern "C" fn(c: circ, data: *mut c_void),
    }

    #[derive(Debug, Copy, Clone)]
    #[repr(C)]
    pub(crate) struct circ {
        ct:   circuit,
        data: *mut c_void,

        env:   circ_env,
        qureg: quest_sys::Qureg,

        mea_cl:      *mut c_int,
        mea_cl_prob: *mut c_double,

        mea_qb: *mut c_int,
        sys_qb: *mut c_int,
        anc_qb: *mut c_int,

        simul_counter: usize,
    }

    extern "C" {

        pub(crate) fn circ_env_init(env: *mut circ_env) -> circ_result;

        pub(crate) fn circ_env_drop(env: circ_env);

        pub(crate) fn circ_env_report(env: circ_env);

        pub(crate) fn circ_init(
            c: *mut circ,
            ct: circuit,
            env: circ_env,
            data: *mut c_void,
        ) -> circ_result;

        pub(crate) fn circ_drop(c: circ);

        pub(crate) fn circ_num_tot_qb(c: circ) -> usize;
        pub(crate) fn circ_report(c: circ);

        pub(crate) fn circ_reset(c: circ) -> circ_result;

        pub(crate) fn circ_simulate(c: circ) -> circ_result;

    }
}

#[derive(Debug)]
pub enum Error {
    Init { msg: String },
}

pub struct CircEnv(ffi::circ_env);

impl CircEnv {
    pub fn try_new() -> Result<Self, Error> {
        let mut env = MaybeUninit::uninit();
        let env_init = unsafe {
            (circ_env_init(env.as_mut_ptr()) == circ_result::CIRC_OK)
                .then(|| env.assume_init())
                .ok_or(Error::Init {
                    msg: "cannot initialize environment".to_string(),
                })
        }?;

        Ok(Self(env_init))
    }

    pub fn report(&self) {
        unsafe {
            ffi::circ_env_report(self.0);
        }
    }
}

impl Drop for CircEnv {
    fn drop(&mut self) {
        unsafe { ffi::circ_env_drop(self.0) }
    }
}

pub struct Circuit<T> {
    pub(crate) circuit: ffi::circuit,
    pub(crate) data:    T,
}

pub struct Circ<'a, C, T> {
    env:  &'a CircEnv,
    ct:   &'a Circuit<C>,
    circ: ffi::circ,
    data: T,
}

impl<'a, C, T> Circ<'a, C, T> {
    pub fn try_new(
        ct: &'a Circuit<C>,
        env: &'a CircEnv,
        data: T,
    ) -> Result<Self, Error> {
        let mut circ_uninit = MaybeUninit::uninit();
        let mut data = data;
        let data_ptr: *mut T = &mut data;

        unsafe {
            ffi::circ_init(
                circ_uninit.as_mut_ptr(),
                ct.circuit,
                env.0,
                mem::transmute(data_ptr),
            )
        };

        let circ = unsafe { circ_uninit.assume_init() };

        Ok(Self {
            env,
            ct,
            circ,
            data,
        })
    }

    pub fn report(&self) {
        unsafe {
            circ_report(self.circ);
        }
    }
}

impl<'a, C, T> Drop for Circ<'a, C, T> {
    fn drop(&mut self) {
        unsafe { ffi::circ_drop(self.circ) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_circ_env() {
        let env = CircEnv::try_new().unwrap();
        env.report();
    }
}
