use std::mem;

use crate::circ::{
    ffi,
    Circuit,
};

#[derive(Debug)]
pub struct LinenCircuit<'a>(pub(crate) Circuit<'a, ffi::linen_circuit_data>);

impl<'a> LinenCircuit<'a> {
    pub fn new(data: &'a mut ffi::linen_circuit_data) -> Self {
        let data_ptr = data as *mut _;
        let mut circuit = unsafe { ffi::LINEN_CIRCUIT_TEMPLATE };
        unsafe {
            circuit.data = mem::transmute(data_ptr);
        }
        Self(Circuit {
            ct: circuit,
            data,
        })
    }
}

impl<'a> AsRef<Circuit<'a, ffi::linen_circuit_data>> for LinenCircuit<'a> {
    fn as_ref(&self) -> &Circuit<'a, ffi::linen_circuit_data> {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circ::ffi::linen_circuit_data;

    #[test]
    fn linen_circuit_init() {
        let mut ct_data = linen_circuit_data {
            state_prep_value: 0,
            routine_value:    11,
            state_post_value: 222,
        };
        let ct = LinenCircuit::new(&mut ct_data);

        assert_eq!(ct.0.data.state_prep_value, 0);
        assert_eq!(ct.0.data.routine_value, 11);
        assert_eq!(ct.0.data.state_post_value, 222);

        let data_ptr: *mut ffi::linen_circuit_data =
            unsafe { mem::transmute(ct.0.ct.data) };

        let data: (i32, i32, i32) = unsafe {
            (
                (*data_ptr).state_prep_value,
                (*data_ptr).routine_value,
                (*data_ptr).state_post_value,
            )
        };

        assert_eq!(data.0, 0);
        assert_eq!(data.1, 11);
        assert_eq!(data.2, 222);
    }
}
