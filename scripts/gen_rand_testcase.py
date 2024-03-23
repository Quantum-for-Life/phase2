#!/bin/env python

import argparse
import json
import math
import random
import uuid
from math import sqrt

import h5py
import numpy as np
import scipy as sp

FILENAME_STUB = "case"

PAULI_TABLE = {
    0: [[1, 0], [0, 1]],
    1: [[0, 1], [1, 0]],
    2: [[0, -1j], [1j, 0]],
    3: [[1, 0], [0, -1]],
}

NUM_STEPS = 128


def pauli_string_to_matrix(paulis: list):
    num_qubits = len(paulis)
    size = 2 ** num_qubits
    matrix = np.ndarray(shape=(size, size), dtype=np.cdouble)
    for i in range(0, size):
        for j in range(0, size):
            elem = 1.0
            for k in range(0, num_qubits):
                ti = i >> k & 1
                tj = j >> k & 1
                elem *= PAULI_TABLE.get(paulis[k])[ti][tj]
            matrix[i][j] = elem

    return matrix


def pauli_hamil_to_matrix(coeffs, paulis):
    num_qubits = len(paulis[0])
    num_terms = len(coeffs)
    size = 2 ** num_qubits
    matrix = np.zeros(shape=(size, size), dtype=np.cdouble)
    for i in range(0, num_terms):
        matrix += coeffs[i] * pauli_string_to_matrix(paulis[i])

    return matrix


def multidet_to_vector(coeffs, indices, num_qubits):
    size = 2 ** num_qubits
    vec = np.zeros(shape=(size,), dtype=np.cdouble)
    for i, c in zip(indices, coeffs):
        vec[i] = c
    return vec


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-n", "--num-qubits", required=True, type=int, help="Number of qubits"
    )
    parser.add_argument(
        "-d",
        "--num-dets",
        required=True,
        type=int,
        help="Number of slater " "determinant",
    )
    parser.add_argument(
        "-t", "--num-terms", required=True, type=int,
        help="Number of hamiltonian terms"
    )
    parser.add_argument(
        "--many", type=int,
        help="Generate many files (parameters as max " "values)"
    )
    parser.add_argument("--compute", action="store_true")

    return parser.parse_args()


class Case:
    def __init__(self, filename_prefix: str, num_qubits: int):
        self.id_str = str(uuid.uuid4())
        self.filename_prefix = filename_prefix + "-" + self.id_str[:8]
        self.num_qubits = num_qubits
        self.pauli_hamil = {}
        self.multidet = {}
        self.trotter_steps = {}

    def generate_pauli_hamil(self, num_terms: int):
        self.pauli_hamil["num_terms"] = num_terms
        terms = set()
        while len(terms) < num_terms:
            term = tuple(
                random.randrange(0, 4) for _ in range(0, self.num_qubits))
            if term != tuple(0 for _ in range(0, self.num_qubits)):
                terms.add(term)

        self.pauli_hamil["paulis"] = np.ndarray(
            shape=(num_terms, self.num_qubits), dtype="u1"
        )
        paulis = [list(x) for x in terms]
        self.pauli_hamil["paulis"][...] = paulis
        coeffs = np.array(
            [(random.random() - 0.5) * 2.0 * 0.1 for _ in
             range(0, num_terms)],
            dtype="d"
        )
        self.pauli_hamil["normalization"] = 1.0
        self.pauli_hamil["coeffs"] = coeffs
        self.pauli_hamil["offset"] = 0.0

    def generate_multidet(self, num_dets: int):
        assert num_dets <= 2 ** self.num_qubits
        self.multidet["num_dets"] = num_dets
        coeffs = [random.random() + 1j * random.random() for _ in
                  range(num_dets)]
        norm = sqrt(sum(abs(x * x.conjugate()) for x in coeffs))
        coeffs_norm = np.array([x / norm for x in coeffs], dtype=np.cdouble)
        self.multidet["coeffs"] = coeffs_norm
        dets = set()
        while len(dets) < num_dets:
            dets.add(random.randrange(0, 2 ** self.num_qubits))
        dets = list(dets)
        self.multidet["indices"] = dets
        self.multidet["dets"] = np.ndarray(
            shape=(num_dets, self.num_qubits), dtype="u1"
        )
        self.multidet["dets"][...] = [
            [dets[i] >> j & 1 for j in range(0, self.num_qubits)]
            for i in range(0, num_dets)
        ]

    def compute_values(self):
        trotter_steps = []
        state = multidet_to_vector(self.multidet["coeffs"],
                                   self.multidet["indices"],
                                   self.num_qubits)
        #print(state)
        state_work = state
        for step in range(NUM_STEPS):
            for k, pauli_str in enumerate(self.pauli_hamil["paulis"]):
                #print(f'step: {step}, Pauli {k}: {pauli_str}, state_work[0]:'
                #      f' {state_work[0]}')
                matrix = pauli_string_to_matrix(pauli_str)
                coeff = self.pauli_hamil["coeffs"][k]
                state_work = math.cos(coeff) * state_work + 1j * math.sin(
                    coeff) * matrix.dot(state_work)
            trotter_steps.append(np.vdot(state, state_work))
        self.trotter_steps["values_comp"] = trotter_steps

    def write_h5file(self, filename=None):

        if not filename:
            filename = self.filename_prefix + ".h5"

        with h5py.File(filename, "w") as f:
            f.attrs["uuid"] = self.id_str
            ph = self.pauli_hamil
            h5_ph = f.create_group("pauli_hamil")
            h5_ph.create_dataset("coeffs", shape=(ph["num_terms"],), dtype="d")[
                ...
            ] = ph["coeffs"]
            h5_ph.create_dataset(
                "paulis", shape=(ph["num_terms"], self.num_qubits), dtype="u1"
            )[...] = ph["paulis"]
            h5_ph.attrs["normalization"] = ph["normalization"]
            h5_ph.attrs["offset"] = ph["offset"]

            md = self.multidet
            h5_md = f.create_group("state_prep/multidet")
            h5_coeffs = h5_md.create_dataset(
                "coeffs", shape=(md["num_dets"], 2), dtype="d"
            )
            coeffs = self.multidet["coeffs"]
            h5_coeffs[...] = [[z.real, z.imag] for z in coeffs]
            h5_md.create_dataset(
                "dets", shape=(md["num_dets"], self.num_qubits), dtype="u1"
            )[...] = md["dets"]
            grp = f.create_group("trotter_steps")
            grp.attrs["time_factor"] = 1.0

    def write_h5file_solved(self):
        filename = self.filename_prefix + ".h5_solved"
        self.write_h5file(filename)
        with h5py.File(filename, 'a') as f:
            h5_ts = f["trotter_steps"]
            h5_ts_values = h5_ts.create_dataset("values", shape=(NUM_STEPS, 2),
                                                dtype="d")
            h5_ts_values[...] = [[z.real, z.imag] for z in
                                 self.trotter_steps["values_comp"]]


if __name__ == "__main__":
    args = parse_arguments()

    case = Case(FILENAME_STUB, args.num_qubits)
    case.generate_pauli_hamil(args.num_terms)
    case.generate_multidet(args.num_dets)
    case.write_h5file()
    if args.compute:
        print("compute trotter steps")
        case.compute_values()
        case.write_h5file_solved()
