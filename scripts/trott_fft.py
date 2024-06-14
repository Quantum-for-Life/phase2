#!/usr/bin/env python

import argparse
from math import tau

import h5py
import numpy as np
import scipy


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

    with h5py.File(args.filename, "r") as f:
        values = np.array(
            [complex(z[0], z[1]) for z in f["circ_trott/values"]])
        time_factor = f["circ_trott"].attrs["time_factor"]
        ph = f["pauli_hamil"]
        norm = ph.attrs["normalization"]
        offset = ph.attrs["offset"]

    fft_size = len(values)
    y_fft = np.fft.fft(values)
    y_fft_abs = [abs(yf) for yf in y_fft]

    peaks = scipy.signal.find_peaks(y_fft_abs)[0]
    fft_freqs = np.fft.fftfreq(fft_size)
    peaks_freqs = [fft_freqs[p] / (norm * time_factor / tau) for p in peaks]
    peaks_freqs.sort()

    print([p + offset for p in peaks_freqs])
