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

impl DataId for ffi::data_id {
    fn is_valid(&self) -> bool {
        *self != ffi::DATA_INVALID_FID
    }
}

pub struct Handle(ffi::data_id);

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
pub struct Data(ffi::data);

impl Data {
    pub fn new() -> Empty<Self> {
        let mut dat_uninit = MaybeUninit::uninit();
        unsafe { ffi::data_init(dat_uninit.as_mut_ptr()) };
        let data = unsafe { dat_uninit.assume_init() };
        Empty(Self(data))
    }

    pub fn for_state_prep(
        &self,
        f: impl FnOnce(&StatePrep),
    ) {
        let state_prep = ManuallyDrop::new(StatePrep(self.0.state_prep));
        f(&state_prep);
    }

    pub fn for_pauli_hamil(
        &self,
        f: impl FnOnce(&PauliHamil),
    ) {
        let pauli_hamil = ManuallyDrop::new(PauliHamil(self.0.pauli_hamil));
        f(&pauli_hamil);
    }

    pub fn for_time_series(
        &self,
        f: impl FnOnce(&TimeSeries),
    ) {
        let time_series = ManuallyDrop::new(TimeSeries(self.0.time_series));
        f(&time_series);
    }
}

impl Drop for Data {
    fn drop(&mut self) {
        unsafe { ffi::data_destroy(&mut self.0 as *mut _) };
    }
}

impl Empty<Data> {
    pub fn parse(
        self,
        fid: &mut Handle,
    ) -> Result<Data, Error> {
        let mut dat = self.0;
        unsafe { ffi::data_parse(&mut dat.0 as *mut _, fid.0) }
            .is_data_ok()
            .then_some(dat)
            .ok_or(Error::FileRead)
    }
}

#[derive(Debug)]
pub struct MultiDet(ffi::data_state_prep_multidet);

impl MultiDet {
    pub fn num_qubits(&self) -> usize {
        self.0.num_qubits
    }

    pub fn num_terms(&self) -> usize {
        self.0.num_terms
    }

    pub fn coeffs(&self) -> &[f64] {
        let slice_ptr =
            slice_from_raw_parts(self.0.coeffs as *const _, self.0.num_terms);
        unsafe { &*slice_ptr }
    }

    pub fn dets(&self) -> &[u8] {
        let slice_ptr = slice_from_raw_parts(
            self.0.dets as *const _,
            self.0.num_terms * self.num_qubits(),
        );
        unsafe { &*slice_ptr }
    }
}

#[derive(Debug)]
pub struct StatePrep(ffi::data_state_prep);

impl StatePrep {
    pub fn for_multidet(
        &self,
        f: impl FnOnce(&MultiDet),
    ) {
        let multidet = ManuallyDrop::new(MultiDet(self.0.multidet));
        f(&multidet);
    }
}

#[derive(Debug)]
pub struct PauliHamil(ffi::data_pauli_hamil);

impl PauliHamil {
    pub fn num_qubits(&self) -> usize {
        self.0.num_qubits
    }

    pub fn num_terms(&self) -> usize {
        self.0.num_terms
    }

    pub fn coeffs(&self) -> &[f64] {
        let slice_ptr =
            slice_from_raw_parts(self.0.coeffs as *const _, self.0.num_terms);
        unsafe { &*slice_ptr }
    }

    pub fn paulis(&self) -> &[u8] {
        let slice_ptr = slice_from_raw_parts(
            self.0.paulis as *const _,
            self.0.num_terms * self.0.num_qubits,
        );
        unsafe { &*slice_ptr }
    }

    pub fn norm(&self) -> f64 {
        self.0.norm
    }
}

#[derive(Debug)]
pub struct TimeSeries(ffi::data_time_series);

impl TimeSeries {
    pub fn num_steps(&self) -> usize {
        self.0.num_steps
    }

    pub fn times(&self) -> &[f64] {
        let slice_ptr =
            slice_from_raw_parts(self.0.times as *const _, self.0.num_steps);
        unsafe { &*slice_ptr }
    }

    pub fn values(&self) -> &[f64] {
        let slice_ptr = slice_from_raw_parts(
            self.0.values as *const _,
            self.0.num_steps * 2,
        );
        unsafe { &*slice_ptr }
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
    fn data_parse_state_prep() {
        Data::new()
            .parse(
                &mut Handle::open(
                    &PathBuf::from(TEST_DAT_DIR).join("./simul_H2_2.h5"),
                )
                .unwrap(),
            )
            .unwrap()
            .for_state_prep(|sp| {
                sp.for_multidet(|md| {
                    assert_eq!(md.num_qubits(), 4);
                    assert_eq!(md.num_terms(), 1);
                    assert_eq!(md.coeffs(), &[1.0]);
                    assert_eq!(md.dets(), &[1, 0, 1, 0]);
                });
            })
    }

    #[test]
    fn data_parse_pauli_hamil() {
        let expected_paulis = [
            0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 3, 3,
            0, 0, 3, 0, 3, 0, 3, 0, 0, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2,
            1, 1, 1, 1, 0, 3, 3, 0, 0, 3, 0, 3, 0, 0, 3, 3,
        ];
        let expected_coeffs = [
            -1.16395, 0.298454, 0.00407158, 0.298454, 0.00407158, 0.0749901,
            0.163535, 0.0857172, 0.0107271, 0.0107271, 0.0107271, 0.0107271,
            0.0857172, 0.073949, 0.0749901,
        ];

        Data::new()
            .parse(
                &mut Handle::open(
                    &PathBuf::from(TEST_DAT_DIR).join("./simul_H2_2.h5"),
                )
                .unwrap(),
            )
            .unwrap()
            .for_pauli_hamil(|ph| {
                assert_eq!(ph.num_qubits(), 4);
                assert_eq!(ph.num_terms(), 15);
                assert_eq!(ph.paulis(), &expected_paulis);
                for (coeff, exp_coeff) in zip(ph.coeffs(), expected_coeffs) {
                    assert!(f64::abs(coeff - exp_coeff) < MARGIN);
                }
            })
    }

    #[test]
    fn data_parse_time_series() {
        let expected_times: Vec<_> = (0..100).map(f64::from).collect();

        Data::new()
            .parse(
                &mut Handle::open(
                    &PathBuf::from(TEST_DAT_DIR).join("./simul_H2_2.h5"),
                )
                .unwrap(),
            )
            .unwrap()
            .for_time_series(|ts| {
                assert_eq!(ts.num_steps(), 100);

                for (time, exp_time) in zip(ts.times(), expected_times) {
                    assert!(f64::abs(time - exp_time) < MARGIN);
                }

                for v in ts.values() {
                    assert!(v.is_nan());
                }
            })
    }
}
