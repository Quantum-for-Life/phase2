use std::{
    mem,
    mem::MaybeUninit,
};

use crate::circ::{
    ffi::circuit,
    Circuit,
};

mod ffi {
    use std::ffi::c_void;

    use crate::circ::ffi::circuit;

    extern "C" {
        pub(crate) fn linen_circuit_factory(data: *mut c_void) -> circuit;

    }
}

pub(crate) fn linen_circuit_factory<T>(data: T) -> Circuit<T> {
    let mut data = data;
    let data_ptr: *mut T = &mut data;
    let circuit =
        unsafe { ffi::linen_circuit_factory(mem::transmute(data_ptr)) };
    Circuit {
        circuit,
        data,
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::circ::{
        Circ,
        CircEnv,
    };

    #[test]
    fn linen_circuit() {
        let env = CircEnv::try_new().unwrap();
        let ct = linen_circuit_factory(0);
        let c = Circ::try_new(&ct, &env, ()).unwrap();
        c.report();
    }
}
