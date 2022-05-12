import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from event_categorize import categorize
from WCM import calSrcLocAvg
from preprocess import noise_remove
from reflection import detect_reflection, correct_reflection

TOT_CHN = 4
BUFFER_SIZE = 1400
Fs = 80000
DEBUG = False

# Get data from .mat file
mat_contents = sio.loadmat('d168.mat')
raw_data = mat_contents['data_muzzle']
data = np.zeros_like(raw_data)
# preprocessed_data = noise_remove(raw_data, DEBUG, Fs)
preprocessed_data = raw_data

for i in range(0,4):
    iDelay, ratioCorr = detect_reflection(preprocessed_data[i,:], DEBUG)
    if iDelay is not None:
        data[i,:] = correct_reflection(preprocessed_data[i,:], iDelay, ratioCorr, DEBUG)
    else:
        data[i,:] = preprocessed_data[i, :]
# print(data.shape)

if DEBUG:
    # Visualize the raw data
    plt.plot(data[0,:], label='1')
    plt.plot(data[1,:])
    plt.plot(data[2,:])
    plt.plot(data[3,:])
    plt.title('Reflection Removed Data Data')
    plt.ylabel('A(t)')
    plt.xlabel('t')
    plt.legend(["Microphone 1", "Microphone 2", "Microphone 3", "Microphone 4"], loc="upper right")
    plt.show()

category = categorize(data, Fs, DEBUG)

if category == 'MuzzleBlast':
    loc = calSrcLocAvg(data)