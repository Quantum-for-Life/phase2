#!/usr/bin/env python

import argparse
import cmath
import math

import h5py


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="cmpsit-rpe",
        description="Robust phase estimation for Composite",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="Input file")
    parser.add_argument("--step-pow", type=int)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    with h5py.File(args.filename, "r") as f:
        grp = f['pauli_hamil']
        norm = grp.attrs['normalization']
        offset = grp.attrs['offset']

    J = args.step_pow

    thetas = [0.0]
    for i in range(0, J + 1):
        filename = args.filename + f"-{i}"
        with h5py.File(filename, "a") as f:
            grp = f["circ_cmpsit"]
            step_size = grp.attrs["step_size"]
            num_samples = len(grp["values"])
            z_re = sum(x[0] for x in grp["values"]) / num_samples
            z_im = sum(x[1] for x in grp["values"]) / num_samples
        z = z_re + 1j * z_im
        t = step_size * 2**i
        print(f"{step_size=}, 2^i={2**i}, {t=}, {z=}")
        phi = cmath.phase(z)  # phi \in [-\pi, \pi]
        if phi < 0:
            phi = 2 * math.pi + phi
        S = [2 ** (-i) * (2 * math.pi * k + phi) for k in range(0, 2 ** i)]
        D = [abs(thetas[i] - S[k]) for k in range(0, 2 ** i)]
        k = 0
        while k < 2 ** i:
            if abs(thetas[i] - S[k]) == min(D):
                break
            k += 1
        thetas.append(S[k])
    thJ = thetas[J]
    if thJ > math.pi:
        thJ = - (2 * math.pi - thJ)
    x = step_size
    E0 = math.sqrt(1-x**2)/x * math.tan(x*thJ) / norm
    E = E0  + offset
    print(f"{E0},{E}")
