#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::ffi::{
    c_char,
    c_double,
    c_uchar,
};

pub(crate) const DATA_INVALID_OBJID: i64 = -1;

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

pub(crate) type dataid_t = i64;

extern "C" {
    pub(crate) fn data_file_open(filename: *const c_char) -> dataid_t;
    pub(crate) fn data_file_close(file_id: dataid_t);
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_state_prep {
    pub(crate) multidet: *mut data_state_prep_multidet,
}

extern "C" {
    pub(crate) fn data_state_prep_init(dat: *mut data_state_prep);

    pub(crate) fn data_state_prep_destroy(dat: *mut data_state_prep);

    pub(crate) fn data_state_prep_parse(
        dat: *mut data_state_prep,
        obj_id: dataid_t,
    ) -> data_result;
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_state_prep_multidet {
    pub(crate) num_qubits: usize,
    pub(crate) num_terms:  usize,
    pub(crate) coeffs:     *mut c_double,
    pub(crate) dets:       *mut c_uchar,
}

extern "C" {
    pub(crate) fn data_state_prep_multidet_init(
        dat: *mut data_state_prep_multidet
    );

    pub(crate) fn data_state_prep_multidet_destroy(
        dat: *mut data_state_prep_multidet
    );

    pub(crate) fn data_state_prep_multidet_parse(
        dat: *mut data_state_prep_multidet,
        obj_id: dataid_t,
    ) -> data_result;
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

extern "C" {
    pub(crate) fn data_pauli_hamil_init(dat: *mut data_pauli_hamil);

    pub(crate) fn data_pauli_hamil_destroy(dat: *mut data_pauli_hamil);

    pub(crate) fn data_pauli_hamil_parse(
        dat: *mut data_pauli_hamil,
        file_id: dataid_t,
    ) -> data_result;
}

#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub(crate) struct data_time_series {
    pub(crate) num_steps: usize,
    pub(crate) times:     *mut c_double,
    pub(crate) values:    *mut c_double,
}

extern "C" {
    pub(crate) fn data_time_series_init(dat: *mut data_time_series);

    pub(crate) fn data_time_series_destroy(dat: *mut data_time_series);

    pub(crate) fn data_time_series_parse(
        dat: *mut data_time_series,
        file_id: dataid_t,
    ) -> data_result;

    pub(crate) fn data_time_series_write(
        dat: *mut data_time_series,
        file_id: dataid_t,
    ) -> data_result;
}
