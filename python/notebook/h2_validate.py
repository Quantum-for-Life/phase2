import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy

eigs = [-1.8382052561377937, -1.6712704942150265, -1.6712704942150265,
        -1.422886756088479, -1.422886756088479, -1.422886756088479,
        -1.3370700060532204, -1.2616776969141434, -1.2616776969141434,
        -1.245391395908435, -1.245391395908435, -1.2101011705960758,
        -0.657553722272001, -0.47745562462465924, -0.4774556246246589,
        -3.469446951953614e-16]

with h5py.File("../build/simul.h5", "r") as f:
    times = np.array(f["time_series/times"])
    values = np.array([complex(z[0], z[1]) for z in f["time_series/values"]])

abs_values_fft = abs(np.fft.fft(values, len(values)))
sample_freq = np.fft.fftfreq(len(abs_values_fft))

peaks = scipy.signal.find_peaks(abs_values_fft)[0]
peak_freqs = [sample_freq[i] for i in peaks]
alpha = (max(eigs) - min(eigs)) / (max(peak_freqs) - min(peak_freqs))
peak_freqs_rescaled = [(x - max(peak_freqs)) * alpha + max(eigs) for x in
                       peak_freqs]
plt.vlines(peak_freqs_rescaled, 0, 1, linestyles='dotted', colors="k")
plt.show()
plt.vlines(eigs, 0, 1, linestyles="dotted", )
plt.show()
