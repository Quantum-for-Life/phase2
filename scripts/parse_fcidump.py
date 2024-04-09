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
    parser.add_argument("--sort-terms", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def fcidump_parse_fermionic_op(
        filename: str, verbose: bool
) -> ElectronicStructureProblem:
    fcidump = FCIDump.from_file(filename)
    if verbose:
        print(f"Number of orbitals: {fcidump.num_orbitals}")

    return fcidump_to_problem(fcidump)


def h5_output(problem: ElectronicStructureProblem, outfile: str,
              sort_terms=False):
    pauli_table = {"I": 0, "X": 1, "Y": 2, "Z": 3}

    fermionic_op = problem.hamiltonian.second_q_op()
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    import numpy as np

    num_qubits = qubit_jw_op.num_qubits
    num_sum_terms = len(qubit_jw_op.coeffs)

    labels = []
    coeffs = np.ndarray((num_sum_terms,), dtype="d")
    paulis = np.ndarray((num_sum_terms, num_qubits), dtype="u1")
    for i in range(0, num_sum_terms):
        coeffs[i] = qubit_jw_op.coeffs[i].real
        labels.append(qubit_jw_op.paulis[i].to_label())
        # NOTE: Qiskit uses LE convention for numbering qubits. We reverse
        # the string
        for j, p in enumerate(qubit_jw_op.paulis[i].to_label()[::-1]):
            paulis[i][j] = pauli_table[p]

    if sort_terms:
        idx = np.argsort(labels)
        paulis = np.array(paulis)[idx]
        coeffs = np.array(coeffs)[idx]

    # offset term
    offset_idx = None
    for i in range(num_sum_terms):
        if all(paulis[i][j] == 0 for j in range(num_qubits)):
            offset_idx = i
    offset = 0.0
    if offset_idx is not None:
        offset = coeffs[offset_idx] + problem.nuclear_repulsion_energy
        coeffs = np.delete(coeffs, offset_idx)
        paulis = np.delete(paulis, 0, offset_idx)
        num_sum_terms -= 1

    with h5py.File(outfile, "w") as f:
        f.attrs["uuid"] = str(uuid.uuid4())
        grp = f.create_group("pauli_hamil")
        dset_coeffs = grp.create_dataset("coeffs", (num_sum_terms,), dtype="d")
        dset_coeffs[...] = coeffs
        dset_paulis = grp.create_dataset(
            "paulis", (num_sum_terms, num_qubits), dtype=np.dtype("u1")
        )
        dset_paulis[...] = paulis
        norm = np.ndarray((1,), dtype="d")
        norm[0] = 1.0 / sum(abs(c) for c in coeffs)
        grp.attrs["normalization"] = norm[0]
        grp.attrs["offset"] = offset

        h5_md = f.create_group("state_prep/multidet")
        h5_coeffs = h5_md.create_dataset(
            "coeffs", shape=(1, 2), dtype="d"
        )
        coeffs = [1 + 0j]
        det0 = []
        for i in range(int(num_qubits / 2)):
            det0.append(1)
            det0.append(0)
        h5_coeffs[...] = [[z.real, z.imag] for z in coeffs]
        h5_md.create_dataset(
            "dets", shape=(1, num_qubits), dtype="u1"
        )[...] = [det0]

        grp = f.create_group("trotter_steps")
        grp.attrs["time_factor"] = 1.0


if __name__ == "__main__":
    args = parse_arguments()
    fermionic_op = fcidump_parse_fermionic_op(args.filename,
                                              verbose=args.verbose)

    if args.output:
        h5_output(fermionic_op, args.output, sort_terms=args.sort_terms)
