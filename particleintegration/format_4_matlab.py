from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
import mapping_functions as mf
#from scipy.io import savemat    
from hdf5storage import savemat
height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km

tdim = 25
xdim = 102
ydim = 82

root = Dataset('wrf_2011_07_01','r')
cen_lat = getattr(root,'CEN_LAT')
cen_lon = getattr(root,'CEN_LON')
true_lat1 = getattr(root,'TRUELAT1')
true_lat2 = getattr(root,'TRUELAT2')
ref_lat = getattr(root,'MOAD_CEN_LAT')
ref_lon = getattr(root,'STAND_LON')
vars = root.variables
#Wind Velocity
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
#Water Vapor Flux, Vertically Integrated
#u = vars['UQ_Q'][:,:,:]
#v = vars['VQ_Q'][:,:,:]
lat_in = vars['XLAT'][0,:,:]
lon_in = vars['XLONG'][0,:,:]
root.close()

root = Dataset('wrf_2011_07_02','r')
vars = root.variables
#Wind Velocity
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
#Water Vapor Flux, Vertically Integrated
#u = np.concatenate((u,vars['UQ'][:,:,:]))
#v = np.concatenate((v,vars['VQ'][:,:,:]))
root.close()
#FLUX DOES NOT NEED UNSTAGGERING
u = mf.unstagger(np.double(u[:25,:,:]),2)
v = mf.unstagger(np.double(v[:25,:,:]),1)
#u = mf.unstagger(u[:25,:,:],2)
#v = mf.unstagger(v[:25,:,:],1)
#u = u[:25,:,:]
#v = v[:25,:,:]
u = np.moveaxis(u,0,-1)
v = np.moveaxis(v,0,-1)
x = np.linspace(0,(xdim-1)*grid_spacing,xdim)
y = np.linspace(0,(ydim-1)*grid_spacing,ydim)
t = np.linspace(0,24,25)
t=t*3600
[x,y,t] = np.meshgrid(x,y,t)
u = u.data
v = v.data
savemat('wrf_vel_data.mat',{'u':u,'v':v,'x':x,'y':y,'time':t},format='7.3')
from matplotlib.pyplot import pcolormesh, colorbar
pcolormesh(x[:,:,1])
colorbar()
