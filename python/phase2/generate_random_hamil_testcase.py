#!/bin/env python

import argparse
import json
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


def pauli_string_to_matrix(paulis: list):
    num_qubits = len(paulis)
    size = 2**num_qubits
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
    size = 2**num_qubits
    matrix = np.zeros(shape=(size, size), dtype=np.cdouble)
    for i in range(0, num_terms):
        matrix += coeffs[i] * pauli_string_to_matrix(paulis[i])

    return matrix


def multidet_to_vector(coeffs, indices, num_qubits):
    size = 2**num_qubits
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
        "-s", "--num-steps", required=True, type=int, help="Time series (num_steps)"
    )
    parser.add_argument(
        "--time-max", type=float, help="Max time value (num_steps if not " "specified)"
    )
    parser.add_argument(
        "-d",
        "--num-dets",
        required=True,
        type=int,
        help="Number of slater " "determinant",
    )
    parser.add_argument(
        "-t", "--num-terms", required=True, type=int, help="Number of hamiltonian terms"
    )
    parser.add_argument(
        "--many", type=int, help="Generate many files (parameters as max " "values)"
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
        self.time_series = {}

    def generate_pauli_hamil(self, num_terms: int):
        self.pauli_hamil["num_terms"] = num_terms
        terms = set()
        while len(terms) < num_terms:
            term = tuple(random.randrange(0, 4) for _ in range(0, self.num_qubits))
            if term != tuple(0 for _ in range(0, self.num_qubits)):
                terms.add(term)

        self.pauli_hamil["paulis"] = np.ndarray(
            shape=(num_terms, self.num_qubits), dtype="u1"
        )
        paulis = [list(x) for x in terms]
        self.pauli_hamil["paulis"][...] = paulis
        coeffs = np.array(
            [(random.random() - 0.5) * 2.0 for _ in range(0, num_terms)], dtype="d"
        )
        self.pauli_hamil["normalization"] = 1.0 / sum((abs(x) for x in coeffs))
        self.pauli_hamil["coeffs"] = coeffs
        self.pauli_hamil["offset"] = (random.random() - 0.5) * 21.0

        # self.pauli_hamil["eivecs"] = [
        #     [[v[i].real, v[i].imag] for i in range(0, len(v))] for v in eivecs
        # ]

    def generate_multidet(self, num_dets: int):
        assert num_dets <= 2**self.num_qubits
        self.multidet["num_dets"] = num_dets
        coeffs = [random.random() + 1j * random.random() for _ in range(num_dets)]
        norm = sqrt(sum(abs(x * x.conjugate()) for x in coeffs))
        coeffs_norm = np.array([x / norm for x in coeffs], dtype=np.cdouble)
        self.multidet["coeffs"] = coeffs_norm
        dets = set()
        while len(dets) < num_dets:
            dets.add(random.randrange(0, 2**self.num_qubits))
        dets = list(dets)
        self.multidet["indices"] = dets
        self.multidet["dets"] = np.ndarray(
            shape=(num_dets, self.num_qubits), dtype="u1"
        )
        self.multidet["dets"][...] = [
            [dets[i] >> j & 1 for j in range(0, self.num_qubits)]
            for i in range(0, num_dets)
        ]

    def generate_time_series(self, num_steps: int, time_max: float | None):
        self.time_series["num_steps"] = num_steps
        self.time_series["times"] = np.ndarray(shape=(num_steps,), dtype="d")
        self.time_series["values"] = np.ndarray(shape=(num_steps, 2), dtype="d")
        for i in range(0, num_steps):
            tmax = num_steps
            if time_max:
                tmax = time_max
            self.time_series["times"][i] = int(random.random() * tmax)
            self.time_series["values"][i][0] = np.nan
            self.time_series["values"][i][1] = np.nan

    def compute_values(self):
        matrix = pauli_hamil_to_matrix(
            self.pauli_hamil["coeffs"], self.pauli_hamil["paulis"]
        )
        self.pauli_hamil["matrix"] = matrix
        eigs, eivecs = np.linalg.eig(matrix)
        for e in eigs:
            assert abs(e.imag) < 10e-10
        eigs_real = [e.real for e in eigs]
        eigs_real.sort()
        self.pauli_hamil["eigs"] = eigs_real
        multidet = multidet_to_vector(
            self.multidet["coeffs"], self.multidet["indices"], self.num_qubits
        )
        hamil = self.pauli_hamil["matrix"] * self.pauli_hamil["normalization"]
        values = []
        for t in self.time_series["times"]:
            unitary = sp.linalg.expm(-1j * t * hamil)
            tmp = np.dot(unitary, multidet)
            values.append(np.dot(multidet.conj().T, tmp))
        self.time_series["values_comp"] = values

    def write_h5file(self):
        with h5py.File(self.filename_prefix + ".h5", "w") as f:
            f.attrs["uuid"] = self.id_str
            ph = self.pauli_hamil
            h5_ph = f.create_group("pauli_hamil")
            h5_ph.create_dataset("coeffs", shape=(ph["num_terms"],), dtype="d")[
                ...
            ] = ph["coeffs"]
            h5_ph.create_dataset(
                "paulis", shape=(ph["num_terms"], self.num_qubits), dtype="d"
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

            ts = self.time_series
            h5_ts = f.create_group("time_series")
            h5_ts.create_dataset("times", shape=(ts["num_steps"],), dtype="d")[
                ...
            ] = self.time_series["times"]
            h5_ts.create_dataset("values", shape=(ts["num_steps"], 2), dtype="d")[
                ...
            ] = self.time_series["values"]

    def write_info(self, compute=False):
        case_info = {
            "filename": self.filename_prefix + ".h5",
            "num_qubits": self.num_qubits,
            "uuid": self.id_str,
            "pauli_hamil": {
                "offset": self.pauli_hamil["offset"],
                "normalization": self.pauli_hamil["normalization"],
                "coeffs": list(self.pauli_hamil["coeffs"]),
                "paulis": [
                    [
                        int(self.pauli_hamil["paulis"][i][j])
                        for j in range(0, self.num_qubits)
                    ]
                    for i in range(0, self.pauli_hamil["num_terms"])
                ],
            },
            "multidet": {
                "coeffs": [[z.real, z.imag] for z in self.multidet["coeffs"]],
                "indices": list(self.multidet["indices"]),
            },
            "time_series": {
                "times": list(self.time_series["times"]),
            },
        }
        if compute:
            case_info["times_series"]["values"] = [
                [z.real, z.imag] for z in self.time_series["values_comp"]
            ]
            case_info["pauli_hamil"]["eigs"] = list(self.pauli_hamil["eigs"])
        with open(self.filename_prefix + ".json", "w") as f:
            json_str = json.dumps(case_info, indent=4, allow_nan=True)
            f.write(json_str)


if __name__ == "__main__":
    args = parse_arguments()

    case = Case(FILENAME_STUB, args.num_qubits)
    case.generate_pauli_hamil(args.num_terms)
    case.generate_multidet(args.num_dets)
    case.generate_time_series(args.num_steps, args.time_max)
    if args.compute:
        print("compute times series values")
        case.compute_values()
    case.write_h5file()
    case.write_info(compute=args.compute)
