from netCDF4 import Dataset
import numpy as np
import pandas as pd

height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m
tdim = 25
xdim = 405
ydim = 325


root = Dataset('hosiendata_wind_velocity.nc','r')
vars = root.variables
#Wind Velocity
u = vars['eastward_vel'][:]
v = vars['northward_vel'][:]
root.close()

u=u[-1,:,:]
v=v[-1,:,:]


x = np.linspace(0,grid_spacing*(xdim-1),xdim)
dx = x[1]-x[0]
y = np.linspace(0,grid_spacing*(ydim-1),ydim)
dy = y[1]-y[0]

dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)


s1 = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = np.min(np.linalg.eig(S)[0])

        else:
            #s1[i,j] = np.nan #np.ma.masked
            s1[i,j] =  999999
s1 = np.ma.masked_where(s1==999999,s1)            
s1 = s1 - s1.max(axis=None)
s1 = s1/s1.min(axis=None)
s1 = s1.filled(np.nan)

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()
data = [s1.ravel()]
name = ['s1']

timelen = ftle.shape[0]
for tt in range(timelen):
    if tt == 0:
        continue
    
    #f = ftle[timelen-1-tt,:,:] - ftle[timelen-1-tt,:,:].min(axis=None)
    f = ftle[-tt-1,:,:] - ftle[-tt-1,:,:].min(axis=None)
    f = f/f.max(axis=None)
    f =  np.ma.filled(f,np.nan)
    data.append(f.ravel())
    name.append('{0:2.3f}hr'.format(tt/6))


#Alldata = pd.DataFrame(np.transpose([s1.ravel(),ftle_2hr.ravel(),ftle_4hr.ravel(),ftle_6hr.ravel()]),columns=['s1','2hr','4hr','6hr'])
Alldata = pd.DataFrame(np.transpose(data),columns=name)
Alldata.corr().to_csv('Correlation_and_stats_v2.csv',mode='w')
B=Alldata.describe()
B.to_csv('Correlation_and_stats_v2.csv',mode='a')

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.pcolormesh(s1)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(f)
plt.colorbar()





