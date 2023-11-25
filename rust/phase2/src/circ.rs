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
    Simulate,
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
pub struct Circuit<'a, T> {
    #[allow(dead_code)]
    ct:   ffi::circuit,
    #[allow(dead_code)]
    data: &'a mut T,
}

pub struct Circ<'a, 'b: 'a, C, T> {
    circ: ffi::circ,
    #[allow(dead_code)]
    env:  &'a CircEnv,
    #[allow(dead_code)]
    ct:   &'a Circuit<'b, C>,
    #[allow(dead_code)]
    data: &'a mut T,
}

impl<'a, 'b: 'a, C, T> Circ<'a, 'b, C, T> {
    pub fn try_new<Ct>(
        env: &'a CircEnv,
        ct: &'a Ct,
        data: &'a mut T,
    ) -> Result<Self, Error>
    where
        Ct: AsRef<Circuit<'b, C>>,
    {
        let data_ptr = data as *mut _;
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

    pub fn simulate(&mut self) -> Result<(), Error> {
        let circ_ptr = &mut self.circ as *mut _;
        match unsafe { ffi::circ_simulate(circ_ptr) } {
            circ_result::CIRC_OK => Ok(()),
            circ_result::CIRC_ERR => Err(Error::Simulate),
        }
    }
}
impl<'a, 'b: 'a, C, T> Drop for Circ<'a, 'b, C, T> {
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
        CircEnv::try_new().unwrap().report();
    }

    #[test]
    fn linen_circ_01() {
        let env = CircEnv::try_new().unwrap();
        let mut ct_data = ffi::linen_circuit_data {
            state_prep_value: 0,
            routine_value:    11,
            state_post_value: 222,
        };
        let ct = LinenCircuit::new(&mut ct_data);
        assert_eq!(ct.0.data.state_prep_value, 0);
        assert_eq!(ct.0.data.routine_value, 11);
        assert_eq!(ct.0.data.state_post_value, 222);

        let mut dat = ffi::linen_circ_data {
            state_prep_value: -1,
            routine_value:    -1,
            state_post_value: -1,
        };
        let mut c = Circ::try_new(&env, &ct, &mut dat).unwrap();
        c.simulate().unwrap();
        drop(c);

        assert_eq!(dat.state_prep_value, 0);
        assert_eq!(dat.routine_value, 11);
        assert_eq!(dat.state_post_value, 222);
    }
}
