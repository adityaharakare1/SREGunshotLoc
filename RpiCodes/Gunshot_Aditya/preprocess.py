import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import butter, lfilter, freqz

import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def noise_remove(raw_data, DEBUG, Fs):
    N = 700
    op = np.zeros((4,700))
    for i in range(4):
        yf = fft(raw_data[i, :])
        if DEBUG:
            plt.plot(raw_data[0,:])
            plt.title('Original Signal')
            plt.ylabel('A(t)')
            plt.xlabel('t')
            plt.show()
            xf = fftfreq(N, 1 / Fs)
            plt.plot(xf, np.abs(yf))
            plt.title('Spectrum of Signal with Noise (-15dB)')
            plt.ylabel('A(f)')
            plt.xlabel('f')
            plt.show()
        y2 = butter_lowpass_filter(raw_data[i,:], 5000, Fs, 5)
        y2f = fft(y2)
        if DEBUG:
            plt.plot(y2)
            plt.title('Filtered Signal')
            plt.ylabel('A(t)')
            plt.xlabel('t')
            plt.show()
            xf = fftfreq(N, 1 / Fs)
            plt.plot(xf, np.abs(y2f))
            plt.title('Spectrum of Filtered Signal')
            plt.ylabel('A(f)')
            plt.xlabel('f')
            plt.show()
        op[i,:] = y2
    return op
    # return raw_data