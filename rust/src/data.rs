use std::{
    ffi::CString,
    mem::MaybeUninit,
    path::Path,
    ptr::slice_from_raw_parts,
};

use crate::data::ffi::data_result::DATA_OK;

mod ffi;

#[derive(Debug)]
pub enum Error {
    FileOpen,
    FileRead,
}

pub struct Empty<T>(T);

pub struct DataHandle(ffi::dataid_t);

impl DataHandle {
    pub fn open(path: &Path) -> Result<Self, Error> {
        let path_str = path.to_str().ok_or(Error::FileOpen)?;
        let filename = CString::new(path_str).map_err(|_| Error::FileOpen)?;
        let file_id = unsafe { ffi::data_file_open(filename.as_ptr()) };
        if file_id == ffi::DATA_INVALID_OBJID {
            return Err(Error::FileOpen);
        }

        Ok(Self(file_id))
    }
}

impl Drop for DataHandle {
    fn drop(&mut self) {
        unsafe {
            ffi::data_file_close(self.0);
        }
    }
}

#[derive(Debug)]
pub struct PauliHamil(ffi::data_pauli_hamil);

impl PauliHamil {
    pub fn new() -> Empty<Self> {
        let mut dat_uninit = MaybeUninit::uninit();
        let dat = unsafe {
            ffi::data_pauli_hamil_init(dat_uninit.as_mut_ptr());
            dat_uninit.assume_init()
        };
        Empty(Self(dat))
    }

    pub fn num_qubits(&self) -> usize {
        self.0.num_qubits
    }

    pub fn num_terms(&self) -> usize {
        self.0.num_terms
    }

    pub fn coeffs(&self) -> &[f64] {
        let slice_ptr = slice_from_raw_parts(self.0.coeffs, self.0.num_terms);
        unsafe { &*slice_ptr }
    }

    pub fn paulis(&self) -> &[u8] {
        let slice_ptr = slice_from_raw_parts(
            self.0.paulis,
            self.0.num_terms * self.0.num_qubits,
        );
        unsafe { &*slice_ptr }
    }
}

impl Empty<PauliHamil> {
    pub fn read(
        mut self,
        handle: &mut DataHandle,
    ) -> Result<PauliHamil, Error> {
        let res = unsafe {
            ffi::data_pauli_hamil_read(&mut self.0 .0 as *mut _, handle.0)
        };
        (res == ffi::data_result::DATA_OK)
            .then_some(self.0)
            .ok_or(Error::FileRead)
    }
}

impl Drop for PauliHamil {
    fn drop(&mut self) {
        unsafe { ffi::data_pauli_hamil_destroy(&mut self.0 as *mut _) };
    }
}

#[cfg(test)]
mod tests {
    use std::{
        iter::zip,
        path::PathBuf,
    };

    use super::*;

    const MARGIN: f64 = 1e-4;

    #[test]
    fn data_file_open() {
        let data_dir = PathBuf::from("./dat");
        let filename = data_dir.join("./simul_H2_2.h5");
        let handle = DataHandle::open(&filename).unwrap();
    }

    #[test]
    fn data_pauli_hamil_read() {
        let data_dir = PathBuf::from("./dat");
        let filename = data_dir.join("./simul_H2_2.h5");
        let mut handle = DataHandle::open(&filename).unwrap();

        let data = PauliHamil::new().read(&mut handle).unwrap();
        assert_eq!(data.num_qubits(), 4);
        assert_eq!(data.num_terms(), 15);

        let expected_coeffs = [
            -1.16395, 0.298454, 0.00407158, 0.298454, 0.00407158, 0.0749901,
            0.163535, 0.0857172, 0.0107271, 0.0107271, 0.0107271, 0.0107271,
            0.0857172, 0.073949, 0.0749901,
        ];
        for (coeff, exp_coeff) in zip(data.coeffs(), expected_coeffs) {
            assert!(f64::abs(coeff - exp_coeff) < MARGIN);
        }

        let expected_paulis = [
            0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 3, 3,
            0, 0, 3, 0, 3, 0, 3, 0, 0, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2,
            1, 1, 1, 1, 0, 3, 3, 0, 0, 3, 0, 3, 0, 0, 3, 3,
        ];
        assert_eq!(data.paulis(), &expected_paulis);
    }
}
