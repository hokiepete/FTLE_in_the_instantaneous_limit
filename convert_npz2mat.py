import numpy as np
import scipy.io as sio

f = np.load('dg_ftle_data_long.npz')
time = f['time']
ftle = f['ftle']
f.close()

sio.savemat('dg_ftle_data_long.mat',{'ftle':ftle,'time':time})
