import numpy as np
import copy, math
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

MaxCorrIter = 4

def detect_reflection(data, DEBUG):
    corrArray = np.correlate(data, data, 'full')  # auto-correlate 'x'
    corrArray = corrArray[(len(data) - 1):]  # keep zero & +ve delays only
    corrArray = corrArray / corrArray[0]  # normalize by 0-delay value
    corrArray.clip(min=0.)  # clip negative values


    iPks, _ = find_peaks(corrArray, height=0.15)

    if DEBUG:
        print(iPks)
        plt.plot(corrArray)
        plt.xlabel("n")
        plt.ylabel("R[n]")
        plt.title("Autocorrelation of Signal")
        plt.show()

    if len(iPks) < 1:
        return None, None

    jPks = np.argmax(corrArray[iPks])  # index of max peak; 0-delay peak absent
    iDelay = iPks[jPks]  # note no. of samples to peak from 0-delay
    ratioCorr = corrArray[iPks[jPks]]  # ratio of correlations for this
    return iDelay, ratioCorr

def correct_reflection(data, iDelay, ratioCorr, DEBUG):

    if ratioCorr <= 0. or ratioCorr > 0.5:
        # nothing to correct, since ratio of correlations is out of range
        return data
    ampRatio = (1 - math.sqrt(1 - 4 * ratioCorr ** 2)) / (2. * ratioCorr)
    y = copy.deepcopy(data)  # start by copying original signal
    if DEBUG:
        plt.plot(data)
    for iter in range(1, MaxCorrIter + 1):
        if iter * iDelay >= len(data):  # no more shifts left in data
            break
        xsub = np.zeros_like(data)  # initialize shifted & scaled data
        # shift 'x' by 'iter' times no. of delayed samples, and scale it by
        # the negative of the amplitude ratio to the power of 'iter'
        xsub[iter * iDelay:] = (-ampRatio) ** iter \
                                    * data[:-(iter * iDelay)]
        y = y + xsub  # add to running 'y' (modified 'x')
        if DEBUG:
            plt.plot(y)
    if DEBUG:
        plt.legend(["Original Signal","Correction 1", "Correction 2", "Correction 3"], loc="upper right")
        plt.xlabel("t")
        plt.ylabel("A(t)")
        plt.title("Reflection Removal")
        plt.show()
    return y