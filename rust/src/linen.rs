use std::mem;

use crate::circ::Circuit;

mod ffi {
    use std::ffi::c_void;

    use crate::circ::ffi::circuit;

    extern "C" {
        pub(crate) fn linen_circuit_factory(data: *mut c_void) -> circuit;

    }
}

#[derive(Debug)]
pub struct LinenCircuit<T>(Circuit<T>);

impl<T> LinenCircuit<T> {
    pub fn new(data: T) -> Self {
        let mut data = data;
        let data_ptr: *mut T = &mut data;
        let circuit =
            unsafe { ffi::linen_circuit_factory(mem::transmute(data_ptr)) };
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
