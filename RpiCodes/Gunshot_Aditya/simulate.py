import scipy.io as sio
import matplotlib.pyplot as plt
from event_categorize import categorize
from WCM import calSrcLocAvg

TOT_CHN = 4
BUFFER_SIZE = 1400
Fs = 80000
DEBUG = False

# Get data from .mat file
mat_contents = sio.loadmat('muzzle_test_14.mat')
data = mat_contents['data_muzzle']
# print(data.shape)

if DEBUG:
    # Visualize the raw data
    plt.plot(data[0,:], label='1')
    plt.plot(data[1,:])
    plt.plot(data[2,:])
    plt.plot(data[3,:])
    plt.show()

category = categorize(data, Fs, DEBUG)

if category == 'MuzzleBlast':
    loc = calSrcLocAvg(data)