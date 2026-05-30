#!/usr/bin/env python

"""Single-point energy estimate for the composite integrator.

Averages the per-sample overlaps <psi|U|psi> stored in
/circ_cmpsit/values, takes the accumulated phase, and maps it to a
ground-state energy estimate of the simulated Hamiltonian:

    E0 = arg(mean(values)) / (steps * angle_det * normalization)

The estimate is accurate when the randomised tail shares the
deterministic time scale (angle_rand ~ angle_det) and the total
phase stays inside (-pi, pi); a larger evolution time wraps the
phase, and a larger angle_rand over-rotates the randomised part.

Prints `E0,E0+offset` (CSV).  E0 is the energy of the normalised,
identity-removed Hamiltonian the simulator evolves; `offset` is the
FCIDUMP scalar term recorded by parse_fcidump.py.
"""

import argparse
import cmath

import h5py
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="cmpsit_rpe",
        description="Composite-integrator energy estimate",
        epilog="Quantum-for-Life",
    )
    parser.add_argument("filename", type=str, help="input file")
    parser.add_argument("--group", type=str, default="circ_cmpsit",
                        help="results group (default: circ_cmpsit)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    with h5py.File(args.filename, "r") as f:
        values = np.array(
            [complex(z[0], z[1]) for z in f[f"{args.group}/values"]])
        grp = f[args.group]
        angle_det = grp.attrs["angle_det"]
        steps = grp.attrs["steps"]
        ph = f["pauli_hamil"]
        norm = ph.attrs["normalization"]
        offset = ph.attrs["offset"]

    z = values.mean()
    phi = cmath.phase(z)
    E0 = phi / (steps * angle_det * norm)
    print(f"{E0},{E0 + offset}")
