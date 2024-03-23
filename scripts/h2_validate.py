from math import tau

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy

with h5py.File("simul.h5", "r") as f:
    times = np.array(f["time_series/times"])
    values = np.array([complex(z[0], z[1]) for z in f["time_series/values"]])
    ph = f["pauli_hamil"]
    norm = ph.attrs["normalization"]
    offset = ph.attrs["offset"]

print(f"{norm=}")
size = len(times)
ts = [(times[i], values[i]) for i in range(size)]
ts.sort()
times = [t for t, _ in ts]
values = [v for _, v in ts]

x = np.array(range(int(min(times)), int(max(times)) + 1))
fft_size = len(x)
y = np.interp(x, times, values)

y_fft = np.fft.fft(y)
y_fft_abs = [abs(yf) for yf in y_fft]

plt.plot(y_fft_abs)
plt.show()

fft_freqs = np.fft.fftfreq(fft_size)
peaks = scipy.signal.find_peaks(y_fft_abs, threshold=fft_size * 0.01)[0]
peaks_freqs = [fft_freqs[p] / norm * tau for p in peaks]
peaks_freqs.sort()

print(peaks_freqs)
print([p + offset for p in peaks_freqs])
