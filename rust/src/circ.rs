#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::{
    mem,
    mem::MaybeUninit,
};

use crate::circ::ffi::{
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

    use crate::{
        quest_sys,
        quest_sys::Qureg,
    };

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
        qureg: Qureg,

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

        pub(crate) fn circ_simulate(c: *mut circ) -> circ_result;

    }
}

#[derive(Debug)]
pub enum Error {
    Init { msg: String },
}

pub struct CircEnv {
    env: ffi::circ_env,
}

impl CircEnv {
    pub fn try_new() -> Result<Self, Error> {
        let mut env_uninit = MaybeUninit::uninit();

        let env = (unsafe { circ_env_init(env_uninit.as_mut_ptr()) }
            == circ_result::CIRC_OK)
            .then(|| unsafe { env_uninit.assume_init() })
            .ok_or(Error::Init {
                msg: "cannot initialize environment".to_string(),
            })?;

        Ok(Self {
            env,
        })
    }

    pub fn report(&self) {
        unsafe {
            ffi::circ_env_report(self.env);
        }
    }
}

impl Drop for CircEnv {
    fn drop(&mut self) {
        unsafe { ffi::circ_env_drop(self.env) }
    }
}

#[derive(Debug)]
pub struct Circuit<T> {
    pub(crate) ct:   ffi::circuit,
    pub(crate) data: T,
}

pub struct Circ<'a, C, T> {
    circ: ffi::circ,
    env:  &'a CircEnv,
    ct:   &'a Circuit<C>,
    data: T,
}

impl<'a, C, T> Circ<'a, C, T> {
    pub fn try_new<Ct>(
        ct: &'a Ct,
        env: &'a CircEnv,
        data: T,
    ) -> Result<Self, Error>
    where
        Ct: AsRef<Circuit<C>>,
    {
        let mut data = data;
        let data_ptr: *mut T = &mut data;

        let mut circ_uninit = MaybeUninit::uninit();
        let circ = (unsafe {
            ffi::circ_init(
                circ_uninit.as_mut_ptr(),
                ct.as_ref().ct,
                env.env,
                mem::transmute(data_ptr),
            )
        } == circ_result::CIRC_OK)
            .then(|| unsafe { circ_uninit.assume_init() })
            .ok_or(Error::Init {
                msg: "circ initialization".to_string(),
            })?;

        Ok(Self {
            env,
            ct: ct.as_ref(),
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
    use crate::linen::LinenCircuit;

    #[test]
    fn report_circ_env() {
        let env = CircEnv::try_new().unwrap();
        env.report();
    }

    #[test]
    fn circ_init() {
        let env = CircEnv::try_new().unwrap();
        let ct = LinenCircuit::new(0);
        let c = Circ::try_new(&ct, &env, 1).unwrap();
    }
}
