import numpy as np
import scipy.io as sio

#f = np.load('dg_ftle_data_long.npz')
f = np.load('wrf_ftle_data.npz')
time = f['time']
ftle = f['ftle']
f.close()

#sio.savemat('dg_ftle_data_long.mat',{'ftle':ftle,'time':time})
sio.savemat('wrf_ftle_data.mat',{'ftle':ftle,'time':time})
