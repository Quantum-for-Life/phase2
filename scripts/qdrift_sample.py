#!/usr/bin/env python

import argparse
import h5py
import shutil
import math


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="qdrift_sample",
        description="Sample Hamiltonian for QDRIFT",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="Input file")
    parser.add_argument("--delta", type=float)
    parser.add_argument("--epsilon", type=float)
    parser.add_argument("--num-samples", type=int)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    with h5py.File(args.filename, "r") as f:
        norm = f['pauli_hamil'].attrs['normalization']

    delta = args.delta
    epsilon = args.epsilon
    J = math.ceil(math.log2(delta / (epsilon * norm)))
    x = math.pow(2, -1 * J)
    print(J)
    print(x)

    for i in range(0, J + 1):
        filename = args.filename + f"-{i:04}"
        shutil.copy(args.filename, filename)
        with h5py.File(filename, "a") as f:
            grp = f.create_group("trotter_steps")
            grp.attrs["step_size"] = x
            grp.attrs["num_samples"] = args.num_samples
            grp.attrs["depth"] = int(math.pow(2, i + J))
