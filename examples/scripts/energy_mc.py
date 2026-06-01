#!/usr/bin/env python

"""Ground-state energy from a Monte-Carlo overlap (qDRIFT / composite).

The randomised drivers (`qdrift`, `cmpsit`) reset the state before
every sample and measure once, so /circ_<algo>/values holds
independent samples of <psi| U |psi> at a single effective evolution
time T -- not a time series.  Each channel reproduces exp(i T H_norm)
in expectation, so the sample mean approximates the overlap and its
phase gives the energy of the largest-overlap eigenstate:

    z  = mean(values)
    E  = arg(z) / (T * normalization) + offset

The effective time depends on the algorithm:

    qdrift:  T = depth * asin(step_size)   (depth gates, angle asin(step_size))
    cmpsit:  T = steps * angle_det         (deterministic-part time scale)

`normalization` and `offset` come from /pauli_hamil.  The estimate is
reliable while the mean stays coherent (|z| not too small): keep T of
order 1 and the per-sample phase inside (-pi, pi).  Large T (small
step_size compensated by large depth, or large angle_rand) decoheres
the average and biases the readout.  Statistical error falls as
1/sqrt(samples); a few hundred samples reach ~0.05 Ha here.

Prints `E,E_ref,dE`: the estimate E = arg(z)/(T*norm) + offset, the
trial-state energy <psi|H|psi>, and their difference dE = E - E_ref
(how far the evolution drove the state below the reference).
"""

import argparse
import cmath
import math

import h5py
import numpy as np

from energy_ref import reference_energy


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="energy_mc",
        description="Monte-Carlo energy from qDRIFT / composite samples",
        epilog="Quantum-for-Life",
    )
    parser.add_argument("filename", type=str, help="input file")
    parser.add_argument("--group", type=str, default="circ_qdrift",
                        help="results group (circ_qdrift or circ_cmpsit)")
    return parser.parse_args()


def effective_time(group, attrs):
    """Effective evolution time T for the algorithm's sample."""
    if group == "circ_qdrift":
        return attrs["depth"] * math.asin(attrs["step_size"])
    if group == "circ_cmpsit":
        return attrs["steps"] * attrs["angle_det"]
    raise ValueError(f"unknown group {group!r}")


if __name__ == "__main__":
    args = parse_arguments()

    with h5py.File(args.filename, "r") as f:
        values = np.array(
            [complex(z[0], z[1]) for z in f[f"{args.group}/values"]])
        T = effective_time(args.group, f[args.group].attrs)
        ph = f["pauli_hamil"]
        norm = ph.attrs["normalization"]
        offset = ph.attrs["offset"]

    z = values.mean()
    e = cmath.phase(z) / (T * norm) + offset
    e_ref = reference_energy(args.filename)
    print(f"{e},{e_ref},{e - e_ref}")
