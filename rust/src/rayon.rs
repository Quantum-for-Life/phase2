use crate::circ::Circuit;

mod ffi {
    use crate::circ::ffi::*;

    #[link(name = "phase2")]
    extern "C" {
        pub(crate) static rayon_circuit: circuit;
    }
}

pub fn rayon_circuit<T>(data: T) -> Circuit<T> {
    Circuit::with_ffi_circuit_factory(unsafe { ffi::rayon_circuit }, data)
}
