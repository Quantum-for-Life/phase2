#!/usr/bin/env python

import argparse
import uuid

import h5py

import numpy as np

from pyscf import gto, scf
from pyscf.mcscf import CASCI
from resource_estimate.hamiltonian_conversion import (
    cas_hamiltonian_coeffs,
    get_pauli_op,
    make_openfermion_QubitOperator,
    pauli_op_trotter_unitary,
)
from resource_estimate.pauli import lexicographic_ordering
from resource_estimate.state_mapping import get_hf_qubit_bitstring


## default value d=0.740848 Ang = 1.4 Bohr and basis set STO-6G
## as in PRX QUANTUM 2, 030305 (2021), Lee et all (tensor hypercontraction paper)

def make_hydrogenchain_xyz(n, d=0.740848):
    xyz = ""
    for i in range(n):
        xyz = xyz + f"H{i}  {i * d}  0   0\n"
    xyz.rstrip("\n")
    return xyz


def get_exact_energy_Hchain(n_atom):
    xyz = make_hydrogenchain_xyz(n_atom)
    mol = gto.M(atom=xyz, basis="sto6g")
    hf = scf.RHF(mol).newton()
    hf.run()
    norb_act = nelec_act = n_atom
    casci = CASCI(hf, norb_act, nelec_act)  # FCI
    e_tot = casci.kernel()[0]
    return e_tot


def setup_integrals_Hchain(n_atom):
    xyz = make_hydrogenchain_xyz(n_atom)
    mol = gto.M(atom=xyz, basis="sto6g")
    hf = scf.RHF(mol).newton()
    hf.run()
    norb_act = nelec_act = n_atom
    h0, h1e_h2o, h2e_h2o = cas_hamiltonian_coeffs(hf, norb_act, nelec_act)
    return h0, h1e_h2o, h2e_h2o


def get_trotter_op_Hchain(n_atom, dt, mapping="symm-bk", order=1):
    if n_atom % 2 != 0:
        raise ValueError("Need even number of atoms")
    nelec = n_atom
    h0, h1e, h2e = setup_integrals_Hchain(n_atom)
    hamiltonian = make_openfermion_QubitOperator(h0, h1e, h2e, mapping, nelec)
    hamiltonian = lexicographic_ordering(hamiltonian)
    pauli_op = get_pauli_op(hamiltonian)
    return pauli_op_trotter_unitary(dt, pauli_op, order)


PAULI_TABLE = {"I": 0, "X": 1, "Y": 2, "Z": 3}


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="parse_fcidump",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("--num-atoms", help="num_atoms")
    parser.add_argument("--delta", help="delta")
    parser.add_argument("--order", help="trotter_order")
    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


class H5Output:
    def __init__(self, pauli_op):

        labels = [op[0] for op in pauli_op]
        coeffs = [op[1] for op in pauli_op]
        self.num_qubits = len(labels[0])
        self.num_sum_terms = len(labels)

        self.labels = []
        self.coeffs = np.ndarray((self.num_sum_terms,), dtype="d")
        self.paulis = np.ndarray((self.num_sum_terms, self.num_qubits), dtype="u1")
        for i in range(0, self.num_sum_terms):
            self.coeffs[i] = coeffs[i]
            self.labels = labels
            # Qiskit uses LE convention for numbering qubits. We reverse the string.
            for j, p in enumerate(labels[i][::-1]):
                self.paulis[i][j] = PAULI_TABLE[p]

        # offset term
        offset_idx = None
        for i in range(self.num_sum_terms):
            if all(self.paulis[i][j] == 0 for j in range(self.num_qubits)):
                offset_idx = i
        self.offset = 0.0
        if offset_idx is not None:
            self.offset = self.coeffs[offset_idx]
            self.coeffs = np.delete(self.coeffs, offset_idx)
            self.paulis = np.delete(self.paulis, 0, offset_idx)
            self.num_sum_terms -= 1

        self.norm = 1.0 / sum(abs(c) for c in self.coeffs)

    def write_h5file(self, outfile: str):
        with h5py.File(outfile, "w") as f:
            f.attrs["uuid"] = str(uuid.uuid4())
            grp = f.create_group("pauli_hamil")

            dset_coeffs = grp.create_dataset("coeffs", (self.num_sum_terms,), dtype="d")
            dset_coeffs[...] = self.coeffs
            dset_paulis = grp.create_dataset(
                "paulis", (self.num_sum_terms, self.num_qubits), dtype=np.dtype("u1")
            )
            dset_paulis[...] = self.paulis

            grp.attrs["normalization"] = self.norm
            grp.attrs["offset"] = self.offset


PAULI_LOOKUP = {
    'I': [[1, 0], [0, 1]],
    'X': [[0, 1], [1, 0]],
    'Y': [[0, -1j], [1j, 0]],
    'Z': [[1, 0], [0, -1]],
}


def paulistr_to_matrix(label):
    num_qb = len(label)
    size = 2 ** num_qb
    matrix = np.ndarray(shape=(size, size), dtype=np.cdouble)
    for i in range(size):
        for j in range(size):
            prod = 1.0 + 1j * 0.0
            for k in range(num_qb):
                i_k = (i >> k) & 1
                j_k = (j >> k) & 1
                prod *= PAULI_LOOKUP[label[k]][i_k][j_k]
            matrix[i][j] = prod

    return matrix

def hamil_to_matrix(pauli_op):
    num_qb = len(pauli_op[0][0])
    size = 2 ** num_qb
    matrix = np.array([[0.0 + 1j * 0.0 for _ in range(size)] for _ in range(size)])
    for label, coeff in pauli_op:
        matrix += coeff * paulistr_to_matrix(label)
    return matrix


import math

if __name__ == "__main__":
    args = parse_arguments()
    num_atoms = int(args.num_atoms)
    pauli_op = get_trotter_op_Hchain(n_atom=num_atoms, dt=1.0, order=int(args.order))
    t = pauli_op[1]
    #pauli_op[1] = ('YY', t[1])

    h5out = H5Output(pauli_op)
    num_qb = len(pauli_op[0][0])

    if num_qb <= 2:

        print(f">>> num_qubits: {num_qb}")
        matrix = hamil_to_matrix(pauli_op)
        eivals, _ = np.linalg.eig(matrix)
        eival_min = sorted(eivals)[0]
        print(f">>> {eival_min=}")

        pauli_op = pauli_op[1:]

        det = get_hf_qubit_bitstring(norb=num_atoms, nelec=num_atoms)
        basis_state = 0
        for k in range(num_qb):
            basis_state += det[k] << k
        state = np.array([0.0 for _ in range(2 ** num_qb)])
        state[basis_state] = 1.0
        comp_state = state
        for i in range(32):
            for label, coeff in pauli_op:
                c = h5out.norm * coeff
                matrix = paulistr_to_matrix(label)
                comp_state = math.cos(c) * comp_state + 1j * math.sin(c) * matrix.dot(comp_state)
            expect_val = np.vdot(state, comp_state)
            print(f"{i=}, {expect_val=}")

    _ = get_exact_energy_Hchain(num_atoms)

    if args.output:
        h5out.write_h5file(args.output)
