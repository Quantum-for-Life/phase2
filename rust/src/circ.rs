#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::{
    mem,
    mem::MaybeUninit,
};

mod ffi;
pub mod linen;
pub mod rayon;

use crate::circ::ffi::{
    circ_env_init,
    circ_report,
    circ_result,
};

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
            ffi::circ_env_report(&self.env as *const _);
        }
    }
}

impl Drop for CircEnv {
    fn drop(&mut self) {
        let env_ptr = &mut self.env as *mut _;
        unsafe { ffi::circ_env_destroy(env_ptr) }
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
        env: &'a CircEnv,
        ct: &'a Ct,
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
                env.env,
                ct.as_ref().ct,
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
            circ_report(&self.circ as *const _);
        }
    }
}
impl<'a, C, T> Drop for Circ<'a, C, T> {
    fn drop(&mut self) {
        let circ_ptr = &mut self.circ as *mut _;
        unsafe { ffi::circ_destroy(circ_ptr) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circ::linen::LinenCircuit;

    #[test]
    fn report_circ_env() {
        let env = CircEnv::try_new().unwrap();
        env.report();
    }

    #[test]
    fn circ_init() {
        let env = CircEnv::try_new().unwrap();
        let ct = LinenCircuit::new(0);
        let c = Circ::try_new(&env, &ct, 1).unwrap();
    }
}
