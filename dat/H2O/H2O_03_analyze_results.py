import math
from math import tau

import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", type=str, nargs='+')
    return parser.parse_args()


if __name__ == "__main__":

    data_files = []
    args = parse_arguments()
    for filename in args.files:
        data = {}
        with h5py.File(filename, "r") as f:
            data["times"] = np.array(f["time_series/times"])
            data["values"] = np.array(
                [complex(z[0], z[1]) for z in f["time_series/values"]])
            ph = f["pauli_hamil"]
            data["norm"] = ph.attrs["normalization"]
            data["offset"] = ph.attrs["offset"]
        data_files.append(data)

    norm_table = [data["norm"] for data in data_files]
    offset_table = [data["offset"] for data in data_files]
    assert all(norm == norm_table[0] for norm in norm_table)
    assert all(offset == offset_table[0] for offset in offset_table)
    norm = norm_table[0]
    offset = offset_table[0]

    times = [t for data in data_files for t in data["times"]]
    values = [v for data in data_files for v in data["values"]]

    size = len(times)
    assert len(values) == size

    ts = [(times[i], values[i]) for i in range(size)]
    ts.sort()
    times = [t for t, _ in ts]
    values = [v for _, v in ts]

    x = np.array(range(int(min(times)), int(max(times)) + 1))
    y = np.interp(x, times, values)

    # Upsampling
    milihartree = 0.001
    ups_fac = math.ceil(tau / (norm * milihartree) / len(x))
    ups_len = len(x) * ups_fac

    ups_spacing = x[-1] / (ups_len - 1)
    x_ups = np.array([x[0] + i * ups_spacing for i in range(ups_len)])
    y_ups = np.interp(x_ups, x, y)

    fft_size = len(x_ups)
    y_fft = np.fft.fft(y_ups)
    y_fft_abs = [abs(yf) for yf in y_fft]
    # plt.plot(y_fft_abs)
    # plt.savefig("y_fft_abs.pdf")

    fft_freqs = np.fft.fftfreq(fft_size)
    peaks = scipy.signal.find_peaks(y_fft_abs, threshold=fft_size * 0.1)[0]
    peaks_freqs = [fft_freqs[p] / norm * tau / ups_spacing for p in peaks]
    peaks_freqs.sort()

    print([p + offset for p in peaks_freqs])
