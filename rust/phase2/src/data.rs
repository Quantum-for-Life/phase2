use std::{
    ffi::CString,
    fmt::{
        Display,
        Formatter,
    },
    mem::{
        ManuallyDrop,
        MaybeUninit,
    },
    path::Path,
    ptr::slice_from_raw_parts,
};

use crate::data::ffi::data_state_prep_multidet_init;

mod ffi;

#[derive(Debug)]
pub enum Error {
    FileOpen,
    FileRead,
    FileWrite,
}

impl Display for Error {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            Self::FileOpen => write!(f, "file open error"),
            Self::FileRead => write!(f, "file read error"),
            Self::FileWrite => write!(f, "file write error"),
        }
    }
}

impl std::error::Error for Error {}

#[derive(Debug, PartialEq)]
pub struct Empty<T>(T);

pub(crate) trait DataId: Sized {
    fn is_valid(&self) -> bool;

    fn valid_then<T>(
        self,
        f: impl FnOnce(Self) -> T,
    ) -> Option<T> {
        self.is_valid().then(|| f(self))
    }
}

impl DataId for ffi::dataid_t {
    fn is_valid(&self) -> bool {
        *self != ffi::DATA_INVALID_OBJID
    }
}

pub struct Handle(ffi::dataid_t);

impl Handle {
    pub fn open(path: &Path) -> Result<Self, Error> {
        let path_str = path.to_str().ok_or(Error::FileOpen)?;
        let filename = CString::new(path_str).map_err(|_| Error::FileOpen)?;
        unsafe { ffi::data_file_open(filename.as_ptr()) }
            .valid_then(Self)
            .ok_or(Error::FileOpen)
    }
}

impl Drop for Handle {
    fn drop(&mut self) {
        unsafe {
            ffi::data_file_close(self.0);
        }
    }
}

#[derive(Debug)]
pub struct StatePrep(ffi::data_state_prep);

impl StatePrep {
    pub fn new() -> Empty<Self> {
        let mut state_prep_uninit = MaybeUninit::uninit();
        unsafe { ffi::data_state_prep_init(state_prep_uninit.as_mut_ptr()) };
        let state_prep = unsafe { state_prep_uninit.assume_init() };

        Empty(Self(state_prep))
    }

    pub fn for_multidet(
        &self,
        f: impl FnOnce(Option<&MultiDet>),
    ) {
        let multidet_ptr = self.0.multidet;
        let multidet = (!multidet_ptr.is_null())
            .then(|| ManuallyDrop::new(MultiDet(unsafe { *multidet_ptr })));

        f(multidet.as_deref());
    }
}

impl Empty<StatePrep> {
    pub fn read(
        mut self,
        handle: &mut Handle,
    ) -> Result<StatePrep, Error> {
        unsafe { ffi::data_state_prep_read(&mut self.0 .0 as *mut _, handle.0) }
            .is_data_ok()
            .then_some(self.0)
            .ok_or(Error::FileRead)
    }
}

impl Drop for StatePrep {
    fn drop(&mut self) {
        unsafe { ffi::data_state_prep_destroy(&mut self.0 as *mut _) };
    }
}

#[derive(Debug)]
pub struct MultiDet(ffi::data_state_prep_multidet);

impl MultiDet {
    pub fn new() -> Empty<Self> {
        let mut dat_uninit = MaybeUninit::uninit();
        unsafe { data_state_prep_multidet_init(dat_uninit.as_mut_ptr()) };
        let dat = unsafe { dat_uninit.assume_init() };

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

    pub fn dets(&self) -> &[u8] {
        let slice_ptr = slice_from_raw_parts(
            self.0.dets,
            self.0.num_terms * self.0.num_qubits,
        );
        unsafe { &*slice_ptr }
    }
}

impl Drop for MultiDet {
    fn drop(&mut self) {
        unsafe { ffi::data_state_prep_multidet_destroy(&mut self.0 as *mut _) };
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

    pub fn norm(&self) -> f64 {
        self.0.norm
    }
}

impl Empty<PauliHamil> {
    pub fn read(
        mut self,
        handle: &mut Handle,
    ) -> Result<PauliHamil, Error> {
        unsafe {
            ffi::data_pauli_hamil_read(&mut self.0 .0 as *mut _, handle.0)
        }
        .is_data_ok()
        .then_some(self.0)
        .ok_or(Error::FileRead)
    }
}

impl Drop for PauliHamil {
    fn drop(&mut self) {
        unsafe { ffi::data_pauli_hamil_destroy(&mut self.0 as *mut _) };
    }
}

#[derive(Debug)]
pub struct TimeSeries(ffi::data_time_series);

impl TimeSeries {
    pub fn new() -> Empty<TimeSeries> {
        let mut dat_uninit = MaybeUninit::uninit();
        unsafe { ffi::data_time_series_init(dat_uninit.as_mut_ptr()) };
        Empty(Self(unsafe { dat_uninit.assume_init() }))
    }

    pub fn write(
        &mut self,
        handle: &mut Handle,
    ) -> Result<(), Error> {
        unsafe { ffi::data_time_series_write(&mut self.0 as *mut _, handle.0) }
            .is_data_ok()
            .then_some(())
            .ok_or(Error::FileWrite)
    }

    pub fn num_steps(&self) -> usize {
        self.0.num_steps
    }

    pub fn times(&self) -> &[f64] {
        let slice_ptr = slice_from_raw_parts(self.0.times, self.0.num_steps);
        unsafe { &*slice_ptr }
    }

    pub fn values(&self) -> &[f64] {
        let slice_ptr =
            slice_from_raw_parts(self.0.values, self.0.num_steps * 2);
        unsafe { &*slice_ptr }
    }
}

impl Drop for TimeSeries {
    fn drop(&mut self) {
        unsafe { ffi::data_time_series_destroy(&mut self.0 as *mut _) };
    }
}

impl Empty<TimeSeries> {
    pub fn read(
        mut self,
        handle: &mut Handle,
    ) -> Result<TimeSeries, Error> {
        unsafe {
            ffi::data_time_series_read(&mut self.0 .0 as *mut _, handle.0)
        }
        .is_data_ok()
        .then_some(self.0)
        .ok_or(Error::FileRead)
    }
}

#[cfg(test)]
mod tests {
    use std::{
        iter::zip,
        path::PathBuf,
    };

    use super::*;

    const TEST_DAT_DIR: &str = "../dat";
    const MARGIN: f64 = 1e-4;

    #[test]
    fn data_file_open() {
        let data_dir = PathBuf::from(TEST_DAT_DIR);
        let filename = data_dir.join("./simul_H2_2.h5");
        Handle::open(&filename).unwrap();
    }

    #[test]
    fn data_pauli_hamil_read() {
        let data_dir = PathBuf::from(TEST_DAT_DIR);
        let filename = data_dir.join("./simul_H2_2.h5");
        let mut handle = Handle::open(&filename).unwrap();

        let data = PauliHamil::new().read(&mut handle).unwrap();
        assert_eq!(data.num_qubits(), 4);
        assert_eq!(data.num_terms(), 15);
        assert!(f64::abs(data.norm() - 2.370806) < MARGIN);

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

    #[test]
    fn data_time_series_read() {
        let data_dir = PathBuf::from(TEST_DAT_DIR);
        let filename = data_dir.join("./simul_H2_2.h5");
        let mut handle = Handle::open(&filename).unwrap();

        let data = TimeSeries::new().read(&mut handle).unwrap();
        assert_eq!(data.num_steps(), 111);

        let expected_times: Vec<_> = (0..111).map(f64::from).collect();
        for (time, exp_time) in zip(data.times(), expected_times) {
            assert!(f64::abs(time - exp_time) < MARGIN);
        }
        for v in data.values() {
            assert!(v.is_nan());
        }
    }

    #[test]
    fn data_state_prep_read() {
        let data_dir = PathBuf::from(TEST_DAT_DIR);
        let filename = data_dir.join("./simul_H2_2.h5");
        let mut handle = Handle::open(&filename).unwrap();

        let state_prep = StatePrep::new().read(&mut handle).unwrap();
        state_prep.for_multidet(|m| {
            if let Some(multidet) = m {
                assert_eq!(multidet.num_qubits(), 4);
                assert_eq!(multidet.num_terms(), 1);

                assert!(f64::abs(multidet.coeffs()[0] - 1.0) < f64::EPSILON);
                assert_eq!(multidet.dets(), &[1, 0, 1, 0]);
            }
        });
    }
}
