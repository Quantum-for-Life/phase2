#!/usr/bin/env python

"""Ground-state energy from a Trotter overlap time series (FFT).

The Trotter drivers (`trott`, `trott2`) record a uniform time series

    values[s] = <psi| U^(s+1) |psi>,   U = exp(i * delta * H_norm)

one point per step at spacing dt = delta.  Each eigenstate |n> in the
reference state contributes a tone exp(i * lambda_n * delta * s), so a
discrete Fourier transform resolves the spectrum: a peak at FFT
frequency f maps to an eigenvalue lambda = tau * f / delta of the
normalised Hamiltonian H_norm, and to a physical energy

    E = lambda / normalization + offset = tau * f / (delta * norm) + offset

(tau = 2*pi).  `normalization` and `offset` come from /pauli_hamil.
np.fft.fftfreq returns signed frequencies, so a negative-eigenvalue
ground state lands at negative f and the formula recovers E < offset
directly -- no sign flip.

The dominant peak (largest |FFT|, excluding the DC bin) is the
eigenstate with the largest reference-state overlap.  For a reference
determinant close to the ground state that peak IS the ground-state
energy.  Frequency resolution is one bin = tau / (steps * delta * norm)
Ha; raise `steps` (or `delta`) until the line sits on a bin.

Prints `E,E_ref,dE`: the dominant-peak energy, the trial-state energy
<psi|H|psi>, and their difference dE = E - E_ref (how far the evolution
drove the state below the reference).  With --peaks, also lists every
detected peak as `E,amplitude` sorted by energy.
"""

import argparse
from math import tau

import h5py
import numpy as np
import scipy.signal

from energy_ref import reference_energy


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="energy_fft",
        description="FFT ground-state energy from a Trotter series",
        epilog="Quantum-for-Life",
    )
    parser.add_argument("filename", type=str, help="input file")
    parser.add_argument("--group", type=str, default="circ_trott",
                        help="results group (circ_trott or circ_trott2)")
    parser.add_argument("--peaks", action="store_true",
                        help="list every detected peak, not just the dominant")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    with h5py.File(args.filename, "r") as f:
        values = np.array(
            [complex(z[0], z[1]) for z in f[f"{args.group}/values"]])
        delta = f[args.group].attrs["delta"]
        ph = f["pauli_hamil"]
        norm = ph.attrs["normalization"]
        offset = ph.attrs["offset"]

    n = len(values)
    amp = np.abs(np.fft.fft(values))
    freqs = np.fft.fftfreq(n)
    energy = freqs / (norm * delta / tau) + offset

    # The DC bin carries the off-resonant baseline; exclude it before
    # picking the dominant tone.
    amp[0] = 0.0
    dom = int(np.argmax(amp))
    e = energy[dom]
    e_ref = reference_energy(args.filename)
    print(f"{e},{e_ref},{e - e_ref}")

    if args.peaks:
        peaks = scipy.signal.find_peaks(amp)[0]
        for p in sorted(peaks, key=lambda b: energy[b]):
            print(f"{energy[p]},{amp[p]}")
