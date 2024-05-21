#!/usr/bin/env python

import argparse
import uuid

import h5py
import qiskit_nature
from qiskit_nature.second_q.formats.fcidump import FCIDump
from qiskit_nature.second_q.formats.fcidump_translator import fcidump_to_problem
from qiskit_nature.second_q.mappers.jordan_wigner_mapper import \
    JordanWignerMapper
from qiskit_nature.second_q.problems import ElectronicStructureProblem

import numpy as np

# Suppress deprecation warnings in fcidump_parse_fermionic_op() below
qiskit_nature.settings.use_symmetry_reduced_integrals = True
qiskit_nature.settings.use_pauli_sum_op = False

PAULI_TABLE = {"I": 0, "X": 1, "Y": 2, "Z": 3}


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="parse_fcidump",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="FCIDUMP input file")
    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("--sort-terms", action="store_true",
                        help="Sort Hamiltonian terms for faster computation")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def fcidump_parse_fermionic_op(
        filename: str, verbose: bool
) -> ElectronicStructureProblem:
    fcidump = FCIDump.from_file(filename)
    if verbose:
        print(f"Number of orbitals: {fcidump.num_orbitals}")

    return fcidump_to_problem(fcidump)


class H5Output:
    def __init__(self, problem: ElectronicStructureProblem,
                 sort_terms=False):
        fermionic_op = problem.hamiltonian.second_q_op()
        mapper = JordanWignerMapper()
        qubit_jw_op = mapper.map(fermionic_op)

        self.num_qubits = qubit_jw_op.num_qubits
        self.num_sum_terms = len(qubit_jw_op.coeffs)

        self.labels = []
        self.coeffs = np.ndarray((self.num_sum_terms,), dtype="d")
        self.paulis = np.ndarray((self.num_sum_terms, self.num_qubits), dtype="u1")
        for i in range(0, self.num_sum_terms):
            self.coeffs[i] = qubit_jw_op.coeffs[i].real
            self.labels.append(qubit_jw_op.paulis[i].to_label())
            # Qiskit uses LE convention for numbering qubits. We reverse the string.
            for j, p in enumerate(qubit_jw_op.paulis[i].to_label()[::-1]):
                self.paulis[i][j] = PAULI_TABLE[p]

        if sort_terms:
            idx = np.argsort(self.labels)
            self.paulis = np.array(self.paulis)[idx]
            self.coeffs = np.array(self.coeffs)[idx]

        # offset term
        offset_idx = None
        for i in range(self.num_sum_terms):
            if all(self.paulis[i][j] == 0 for j in range(self.num_qubits)):
                offset_idx = i
        self.offset = 0.0
        if offset_idx is not None:
            self.offset = self.coeffs[offset_idx] + problem.nuclear_repulsion_energy
            self.coeffs = np.delete(self.coeffs, offset_idx)
            self.paulis = np.delete(self.paulis, 0, offset_idx)
            self.num_sum_terms -= 1

    def write_h5file(self, outfile: str):
        with h5py.File(outfile, "w") as f:
            f.attrs["uuid"] = str(uuid.uuid4())
            grp = f.create_group("pauli_hamil")
            #pos_offst = min(self.coeffs)
            #if pos_offst < 0:

            #    pos_offst = -pos_offst
            #else:
            #    pos_offst = 0
            # self.coeffs = [x + pos_offst for x in self.coeffs]
            dset_coeffs = grp.create_dataset("coeffs", (self.num_sum_terms,), dtype="d")
            dset_coeffs[...] = self.coeffs
            dset_paulis = grp.create_dataset(
                "paulis", (self.num_sum_terms, self.num_qubits), dtype=np.dtype("u1")
            )
            dset_paulis[...] = self.paulis
            norm = np.ndarray((1,), dtype="d")
            norm[0] = 1.0 / sum(abs(c) for c in self.coeffs)
            grp.attrs["normalization"] = norm[0]
            grp.attrs["offset"] = self.offset
            #grp.attrs["positive_shift"] = pos_offst


if __name__ == "__main__":
    args = parse_arguments()
    fermionic_op = fcidump_parse_fermionic_op(args.filename,
                                              verbose=args.verbose)
    h5out = H5Output(fermionic_op,
                     sort_terms=args.sort_terms)

    if args.output:
        h5out.write_h5file(args.output)
