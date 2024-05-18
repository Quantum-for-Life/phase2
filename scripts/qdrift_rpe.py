#!/usr/bin/env python

import argparse
import h5py
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
        grp = f['pauli_hamil']
        norm = grp.attrs['normalization']
        offset = grp.attrs['offset']
        pos_shift = grp.attrs['positive_shift']

    delta = args.delta
    epsilon = args.epsilon
    J = math.ceil(math.log2(delta / epsilon))
    x = math.pow(2, -1 * J)

    thetas = [0.0]
    # zets = []
    # phis = []
    for i in range(0, J + 1):
        filename = args.filename + f"-{i:04}"
        with h5py.File(filename, "a") as f:
            grp = f["trotter_steps"]
            depth = grp.attrs["depth"]
            z_reim = grp["values"][0]
            z = z_reim[0] + 1j * z_reim[1]
            phi = math.acos((z / abs(z)).real)
            S = [2 ** (-i) * (2 * math.pi * k + phi) for k in range(0, 2 ** i)]
            # print(S)
            D = [abs(thetas[i] - S[k]) for k in range(0, 2 ** i)]
            # print(D)
            k = 0
            while k < 2 ** i:
                if abs(thetas[i] - S[k]) == min(D):
                    break
                k += 1
            thetas.append(S[k])

    E0 = (math.sqrt(1 - x * x)) / (norm * x) * math.tan(x * thetas[J])
    print(E0)
    print(E0 - pos_shift)
    print(E0 + offset)
