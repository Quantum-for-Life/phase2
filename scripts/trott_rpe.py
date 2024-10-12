#!/usr/bin/env python

import math
import cmath

import argparse

import h5py

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="trott_validate",
        description="",
        epilog="Quantum-for-Life",
    )
    parser.add_argument("filename", type=str, help="input file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    samples = []
    with h5py.File(args.filename, "r") as f:
        dset_values = f["circ_trott/values"]
        J = int(math.log2(len(dset_values)))
        for i in range(0, J + 1):
            val = dset_values[2 ** i - 1]
            samples.append(complex(val[0], val[1]))

        time_factor = f["circ_trott"].attrs["time_factor"]
        ph = f["pauli_hamil"]
        norm = ph.attrs["normalization"]
        offset = ph.attrs["offset"]

    # print(J)
    thetas = [0.0]
    for i in range(0, J + 1):
        # print(i)
        z = samples[i]
        phi = cmath.phase(z)  # phi \in [-\pi, \pi]
        if phi < 0:
            phi = 2 * math.pi + phi
        S = [2 ** (-i) * (2 * math.pi * k + phi) for k in range(0, 2 ** i)]
        D = (abs(thetas[i] - S[k]) for k in range(0, 2 ** i))
        minD = min(D)
        k = 0
        while k < 2 ** i:
            if abs(thetas[i] - S[k]) == minD:
                break
            k += 1
        thetas.append(S[k])
    thJ = thetas[J]
    if thJ > math.pi:
        thJ = - (2 * math.pi - thJ)
    x = 2**(-J)
    E0 = math.sqrt(1-x**2)/x / norm * math.tan(x*thJ)
    #E0 = thJ / (norm * time_factor)
    E0 = E0 / time_factor
    print(E0, E0 + offset)

