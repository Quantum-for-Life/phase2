use std::{
    ffi::CStr,
    mem,
    ops::DerefMut,
    str::Utf8Error,
};

use crate::circ::ffi::{
    circ_name,
    circ_report,
    circ_reset,
    circ_simulate,
};

pub(crate) mod ffi {
    use std::ffi::{
        c_char,
        c_int,
        c_void,
    };

    #[derive(Debug, Copy, Clone)]
    #[repr(C)]
    #[allow(non_camel_case_types)]
    pub(crate) enum circ_result {
        CIRC_OK,
        CIRC_ERR,
    }

    #[allow(non_camel_case_types)]
    pub(crate) type circ_env = c_void;

    #[allow(non_camel_case_types)]
    pub(crate) type circ = c_void;

    #[derive(Debug, Copy, Clone)]
    #[repr(C)]
    #[allow(non_camel_case_types)]
    pub(crate) struct circuit {
        pub(crate) name: *const c_char,
        pub(crate) data: *mut c_void,

        pub(crate) num_mea_cl: usize,
        pub(crate) num_mea_qb: usize,
        pub(crate) num_sys_qb: usize,
        pub(crate) num_anc_qb: usize,

        reset: extern "C" fn(*mut circ) -> circ_result,

        state_prep: extern "C" fn(*mut circ, *mut c_void) -> circ_result,
        routine:    extern "C" fn(*mut circ, *mut c_void) -> circ_result,
        state_post: extern "C" fn(*mut circ, *mut c_void) -> circ_result,
        measure:    extern "C" fn(*mut circ, *mut c_void) -> circ_result,
    }

    #[link(name = "QuEST")]
    #[link(name = "circ")]
    #[allow(non_snake_case)]
    extern "C" {
        pub(crate) fn circ_create_env() -> *mut circ_env;
        pub(crate) fn circ_destroy_env(env: *mut circ_env);
        pub(crate) fn circ_report_env(env: *mut circ_env);

        pub(crate) fn circ_create(
            ct: circuit,
            env: *const circ_env,
            data: *mut c_void,
        ) -> *mut circ;
        pub(crate) fn circ_destroy(c: *mut circ);

        pub(crate) fn circ_mea_cl(c: *mut circ) -> *mut c_int;
        pub(crate) fn circ_mea_qb(c: *mut circ) -> *mut c_int;
        pub(crate) fn circ_sys_qb(c: *mut circ) -> *mut c_int;
        pub(crate) fn circ_anc_qb(c: *mut circ) -> *mut c_int;

        pub(crate) fn circ_num_mea_cl(c: *mut circ) -> usize;
        pub(crate) fn circ_num_mea_qb(c: *mut circ) -> usize;
        pub(crate) fn circ_num_sys_qb(c: *mut circ) -> usize;
        pub(crate) fn circ_num_anc_qb(c: *mut circ) -> usize;
        pub(crate) fn circ_num_tot_qb(c: *mut circ) -> usize;

        pub(crate) fn circ_name(c: *mut circ) -> *const c_char;

        pub(crate) fn circ_report(c: *mut circ);

        pub(crate) fn circ_circuit_data(c: *mut circ) -> *mut c_void;

        // Qureg circ_qureg(circ *);

        pub(crate) fn circ_reset(c: *mut circ) -> circ_result;

        pub(crate) fn circ_simulate(c: *mut circ) -> circ_result;
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Error {
    Reset,
    Simulate,
}

pub struct CircEnv(*mut ffi::circ_env);

impl CircEnv {
    pub fn new() -> Self {
        Self(unsafe { ffi::circ_create_env() })
    }

    pub fn report(&self) {
        unsafe { ffi::circ_report_env(self.0) };
    }
}

impl Drop for CircEnv {
    fn drop(&mut self) {
        unsafe { ffi::circ_destroy_env(self.0) }
    }
}

pub struct Circ<'e, T> {
    circ: *mut ffi::circ,
    env:  &'e CircEnv,
    data: T,
}

impl<'e, T> Circ<'e, T> {
    pub fn new<K>(
        ct: Circuit<K>,
        env: &'e CircEnv,
        data: T,
    ) -> Self {
        let mut data = data;
        let circ = {
            let data_ptr: *mut T = &mut data;
            unsafe { ffi::circ_create(ct.ct, env.0, mem::transmute(data_ptr)) }
        };

        Self {
            circ,
            env,
            data,
        }
    }

    pub fn name(&self) -> Result<&str, Utf8Error> {
        let name_ptr = unsafe { circ_name(self.circ) };
        unsafe { CStr::from_ptr(name_ptr) }.to_str()
    }

    pub fn report(&self) {
        unsafe { circ_report(self.circ) }
    }

    pub fn data(&self) -> &T {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut T {
        &mut self.data
    }

    pub fn reset(&mut self) -> Result<(), Error> {
        match unsafe { circ_reset(self.circ) } {
            ffi::circ_result::CIRC_OK => Ok(()),
            ffi::circ_result::CIRC_ERR => Err(Error::Reset),
        }
    }

    pub fn simulate(&mut self) -> Result<(), Error> {
        match unsafe { circ_simulate(self.circ) } {
            ffi::circ_result::CIRC_OK => Ok(()),
            ffi::circ_result::CIRC_ERR => Err(Error::Simulate),
        }
    }
}

impl<'e, T> Drop for Circ<'e, T> {
    fn drop(&mut self) {
        unsafe {
            ffi::circ_destroy(self.circ);
        }
    }
}

#[derive(Copy, Clone)]
pub struct Circuit<T> {
    ct:   ffi::circuit,
    data: T,
}

impl<T> Circuit<T> {
    pub(crate) fn with_ffi_circuit_factory(
        ct: ffi::circuit,
        data: T,
    ) -> Self {
        let mut ct = ct;
        let mut data = data;
        let data_ptr: *mut T = &mut data;
        ct.data = unsafe { mem::transmute(data_ptr) };

        Self {
            ct,
            data,
        }
    }

    pub fn data(&self) -> &T {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut T {
        &mut self.data
    }
}
