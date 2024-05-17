#!/usr/bin/env python

import argparse
import h5py
import shutil
import random

import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="qdrift_sample",
        description="Sample Hamiltonian for QDRIFT",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="Input file")
    parser.add_argument("--samples", type=int, default=128, help="Number of sampled terms")
    parser.add_argument("--outs", type=int, default=16, help="Number of output files")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    with  h5py.File(args.filename, "r") as f:
        hamil = f['/pauli_hamil']

        coeffs = list(hamil['coeffs'])
        norm = hamil.attrs['normalization']
        paulis = list(hamil['paulis'])

    num_qubits = len(paulis[0])
    weigths = [abs(x) * norm for x in coeffs]
    popul = range(0, len(weigths))
    num_samples = args.samples

    for i in range(1, args.outs + 1):
        filename = args.filename + f".{i:04}"
        shutil.copy(args.filename, filename)
        with h5py.File(filename, "a") as f:
            del f["/pauli_hamil/coeffs"]
            del f["/pauli_hamil/paulis"]
            grp = f["/pauli_hamil"]
            idx = random.choices(popul, weigths, k=num_samples)
            dset_coeffs = grp.create_dataset("coeffs", (num_samples,), dtype="d")
            dset_coeffs[...] = np.array(coeffs)[idx]
            dset_paulis = grp.create_dataset(
                "paulis", (num_samples, num_qubits), dtype=np.dtype("u1")
            )
            dset_paulis[...] = np.array(paulis)[idx]
