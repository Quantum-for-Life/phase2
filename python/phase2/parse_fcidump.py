#!/usr/bin/env python

import argparse

import h5py
import qiskit_nature
from qiskit_nature.second_q.formats.fcidump import FCIDump
from qiskit_nature.second_q.formats.fcidump_translator import (
    fcidump_to_problem,
)
from qiskit_nature.second_q.mappers.jordan_wigner_mapper import (
    JordanWignerMapper,
)
from qiskit_nature.second_q.operators import FermionicOp

# Suppress deprecation warnings in fcidump_parse_fermionic_op() below
qiskit_nature.settings.use_symmetry_reduced_integrals = True
qiskit_nature.settings.use_pauli_sum_op = False

MARGIN: float = 10e-10


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="jw_map",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="FCIDUMP input file")
    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def fcidump_parse_fermionic_op(filename: str, verbose: bool) -> FermionicOp:
    fcidump = FCIDump.from_file(filename)
    if verbose:
        print(f"Number of orbitals: {fcidump.num_orbitals}")

    problem = fcidump_to_problem(fcidump)
    return problem.hamiltonian.second_q_op()


def h5_output(fermionic_op: FermionicOp, outfile: str):
    pauli_table = {"I": 0, "X": 1, "Y": 2, "Z": 3}

    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    import numpy as np
    # eigs = [x.real for x in np.linalg.eig(qubit_jw_op.to_matrix())[0]]
    # eigs.sort()
    # print(eigs)

    num_qubits = qubit_jw_op.num_qubits
    num_sum_terms = len(qubit_jw_op.coeffs)
    coeffs = []
    paulis_shape = (num_sum_terms, num_qubits)
    paulis = np.ndarray(paulis_shape, dtype='f')
    for i in range(0, num_sum_terms):
        coeffs.append(qubit_jw_op.coeffs[i].real)
        # NOTE: Qiskit uses LE convention for numbering qubits. We reverse
        # the string
        for j, p in enumerate(qubit_jw_op.paulis[i].to_label()[::-1]):
            paulis[i][j] = pauli_table[p]

    with h5py.File(outfile, "w") as f:
        grp = f.create_group("pauli_hamil")
        dset_coeffs = grp.create_dataset("coeffs", (num_sum_terms,), dtype='d')
        dset_coeffs[...] = coeffs
        dset_paulis = grp.create_dataset("paulis",
                                         paulis_shape,
                                         dtype=np.dtype('u1'))
        dset_paulis[...] = paulis
        grp.attrs["normalization"] = sum(abs(c) for c in coeffs)


if __name__ == "__main__":
    args = parse_arguments()
    fermionic_op = fcidump_parse_fermionic_op(
        args.filename, verbose=args.verbose
    )

    if args.output:
        h5_output(fermionic_op, args.output)
