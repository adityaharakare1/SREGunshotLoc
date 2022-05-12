import numpy as np
from event_categorize import categorize
import scipy.io as sio

mat_contents = sio.loadmat('data_muzzle_6.mat')
data = mat_contents['data_muzzle']
# mat_contents = sio.loadmat('Noise1.mat')
# data = mat_contents['xt_n']
# mat_contents = sio.loadmat('NW1.mat')
# data = mat_contents['shockwave']


TOT_CHN = 4
BUFFER_SIZE = 1400
Fs = 80000
DEBUG = True
category = categorize(data, Fs, DEBUG)
