#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::ffi::{
    c_char,
    c_double,
    c_uchar,
};

pub(crate) const DATA_INVALID_FID: i64 = -1;

#[derive(Debug, PartialEq)]
#[repr(C)]
pub(crate) enum data_result {
    DATA_ERR = -1,
    DATA_OK = 0,
}

impl data_result {
    pub(crate) fn is_data_ok(&self) -> bool {
        *self == Self::DATA_OK
    }
}

pub(crate) type data_id = i64;

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_state_prep_multidet {
    pub(crate) num_qubits: usize,
    pub(crate) num_terms:  usize,
    pub(crate) coeffs:     *mut c_double,
    pub(crate) dets:       *mut c_uchar,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_state_prep {
    pub(crate) multidet: data_state_prep_multidet,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_pauli_hamil {
    pub(crate) num_qubits: usize,
    pub(crate) num_terms:  usize,
    pub(crate) coeffs:     *mut c_double,
    pub(crate) paulis:     *mut c_uchar,
    pub(crate) norm:       c_double,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_time_series {
    pub(crate) num_steps: usize,
    pub(crate) times:     *mut c_double,
    pub(crate) values:    *mut c_double,
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data {
    pub(crate) state_prep:  data_state_prep,
    pub(crate) pauli_hamil: data_pauli_hamil,
    pub(crate) time_series: data_time_series,
}

extern "C" {
    pub(crate) fn data_file_open(filename: *const c_char) -> data_id;
    pub(crate) fn data_file_close(file_id: data_id);

    pub(crate) fn data_init(dat: *mut data);

    pub(crate) fn data_destroy(dat: *mut data);

    pub(crate) fn data_parse(
        dat: *mut data,
        obj_id: data_id,
    ) -> data_result;

    pub(crate) fn data_time_series_write(
        dat: *mut data_time_series,
        file_id: data_id,
    ) -> data_result;
}
