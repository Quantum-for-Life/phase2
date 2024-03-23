import math
from math import tau

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy

with h5py.File("./simul.h5", "r") as f:
    values = np.array([complex(z[0], z[1]) for z in f["trotter_steps/values"]])
    time_factor = f["trotter_steps"].attrs["time_factor"]
    ph = f["pauli_hamil"]
    norm = ph.attrs["normalization"]
    offset = ph.attrs["offset"]

print(f"{norm=}")
print(values)

fft_size = len(values)

y_fft = np.fft.fft(values)
y_fft_abs = [abs(yf) for yf in y_fft]
#
# plt.plot(y_fft_abs)
# plt.show()

fft_freqs = np.fft.fftfreq(fft_size)
peaks = scipy.signal.find_peaks(y_fft_abs)[0]
peaks_freqs = [fft_freqs[p] / norm / time_factor * tau for p in peaks]
peaks_freqs.sort()

print(peaks_freqs)
print([p + offset for p in peaks_freqs])
