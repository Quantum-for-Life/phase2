use std::{
    mem,
    mem::MaybeUninit,
};

use crate::circ::Circuit;

mod ffi {
    use std::ffi::c_void;

    use crate::circ::ffi::circuit;

    extern "C" {
        pub(crate) fn linen_circuit_init(
            ct: *mut circuit,
            data: *mut c_void,
        );

    }
}

pub(crate) fn linen_circuit_init<T>(data: T) -> Circuit<T> {
    let mut data = data;
    let data_ptr: *mut T = &mut data;
    let mut circuit_uninit = MaybeUninit::uninit();
    unsafe {
        ffi::linen_circuit_init(
            circuit_uninit.as_mut_ptr(),
            mem::transmute(data_ptr),
        )
    };
    Circuit {
        circuit: unsafe { circuit_uninit.assume_init() },
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
        let ct = linen_circuit_init(0);
        let c = Circ::try_new(&ct, &env, ()).unwrap();
        c.report();
    }
}
