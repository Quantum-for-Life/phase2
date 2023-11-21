use std::mem;

use crate::circ::{
    ffi,
    Circuit,
};

#[derive(Debug)]
pub struct LinenCircuit<T>(Circuit<T>);

impl<T> LinenCircuit<T> {
    pub fn new(data: T) -> Self {
        let mut data = data;
        let data_ptr = &mut data as *mut _;
        let mut circuit = unsafe { ffi::linen_circuit };
        unsafe {
            circuit.data = mem::transmute(data_ptr);
        }
        Self(Circuit {
            ct: circuit,
            data,
        })
    }
}

impl<T> AsRef<Circuit<T>> for LinenCircuit<T> {
    fn as_ref(&self) -> &Circuit<T> {
        &self.0
    }
}
