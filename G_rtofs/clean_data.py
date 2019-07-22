import numpy as np
from hdf5storage import loadmat, savemat
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
f = loadmat('rtofs_vel.mat')
u = f['u']
v = f['v']
x = f['x']
y = f['y']
t = f['t']
ydim,xdim,tdim = u.shape

plt.subplot(211)
plt.pcolormesh(u[:,:,-1])
uu = np.empty(u.shape)
uu[:]=np.nan
vv = np.empty(v.shape)
vv[:]=np.nan
for i in range(1,ydim-1):
    print(i)
    for j in range(1,xdim-1):
        for k in range(tdim):
            
            if ~np.isnan(u[i,j,k]):
                if (~np.isnan(u[i-1,j,k]) or ~np.isnan(u[i+1,j,k])) and (~np.isnan(u[i,j-1,k]) or ~np.isnan(u[i,j+1,k])):
                    uu[i,j,k] = u[i,j,k]
            '''
            if (~np.isnan(u[i,j,k]) and (np.isnan(u[i-1,j,k]) and np.isnan(u[i+1,j,k])) or (np.isnan(u[i,j-1,k]) and np.isnan(u[i,j+1,k]))):
                u[i,j,k] = np.nan
            if (~np.isnan(v[i,j,k]) and (np.isnan(v[i-1,j,k]) and np.isnan(v[i+1,j,k])) or (np.isnan(v[i,j-1,k]) and np.isnan(v[i,j+1,k]))):
                v[i,j,k] = np.nan
            '''
#savemat('rtofs_clean_vel.mat',{'x':x,'y':y,'t':t,'u':u,'v':v})
plt.subplot(212)
plt.pcolormesh(uu[:,:,-1])