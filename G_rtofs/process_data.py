import mapping_functions as mf
from netCDF4 import Dataset
import numpy as np
from hdf5storage import savemat
from scipy.interpolate import griddata

#"""
#For plots
#25N
#-90E
xdim=101
ydim=xdim#int(np.ceil(xdim*4/5))
tdim = 73
lon_0 = -90
lat_0 = 27.5
std_lat1=25
std_lat2=30


u_out = np.empty((ydim,xdim,tdim))
v_out = np.empty((ydim,xdim,tdim))
t = np.empty((tdim))

x_out = np.linspace(-500*1000,500*1000,xdim)
y_out = np.linspace(-600*1000,400*1000,ydim)
y_low_thresh = y_out.min()-20*1000
x_low_thresh = x_out.min()-20*1000
y_hi_thresh = y_out.max()+20*1000
x_hi_thresh = x_out.max()+20*1000

xx,yy = np.meshgrid(x_out,y_out)
Xi= list(zip(xx.ravel(),yy.ravel()))
del xx, yy
"""
#For integration
#25N
#-90E
xdim=201
ydim=xdim#int(np.ceil(xdim*4/5))
tdim = 73
lon_0 = -90
lat_0 = 27.5
std_lat1=25
std_lat2=30


u_out = np.empty((ydim,xdim,tdim))
v_out = np.empty((ydim,xdim,tdim))
t = np.empty((tdim))

x_out = np.linspace(-1000*1000,1000*1000,xdim)
y_out = np.linspace(-1100*1000,900*1000,ydim)
y_low_thresh = y_out.min()-20*1000
x_low_thresh = x_out.min()-20*1000
y_hi_thresh = y_out.max()+20*1000
x_hi_thresh = x_out.max()+20*1000

xx,yy = np.meshgrid(x_out,y_out)
Xi= list(zip(xx.ravel(),yy.ravel()))
del xx, yy
#"""
for i in range(tdim):
    print(i)
    root = Dataset('data/07_25_2019/rtofs_glo_2ds_f{0:03d}_1hrly_prog.nc'.format(i),'r')
    vars = root.variables #dictionary, all variables in dataset\
    t[i] = vars['MT'][0]
    lon = vars['Longitude'][:-1,:]-360
    lat = vars['Latitude'][:-1,:]
    u = vars['u_velocity'][:,:,:-1,:].squeeze()
    v = vars['v_velocity'][:,:,:-1,:].squeeze()
    root.close()
    
    
    x_in,y_in = mf.lonlat2m(lon_0,lat_0,lon,lat,std_lat1=std_lat1,std_lat2=std_lat2)
    del lon, lat
    
    u = u[x_in>=x_low_thresh]
    v = v[x_in>=x_low_thresh]
    y_in = y_in[x_in>=x_low_thresh]
    x_in = x_in[x_in>=x_low_thresh]
    
    u = u[x_in<=x_hi_thresh]
    v = v[x_in<=x_hi_thresh]
    y_in = y_in[x_in<=x_hi_thresh]
    x_in = x_in[x_in<=x_hi_thresh]
    
    u = u[y_in>=y_low_thresh]
    v = v[y_in>=y_low_thresh]
    x_in = x_in[y_in>=y_low_thresh]
    y_in = y_in[y_in>=y_low_thresh]
    
    u = u[y_in<=y_hi_thresh]
    v = v[y_in<=y_hi_thresh]
    x_in = x_in[y_in<=y_hi_thresh]
    y_in = y_in[y_in<=y_hi_thresh]
    
    v[u.mask] = np.nan
    u[u.mask] = np.nan
    
    u[v.mask] = np.nan
    v[v.mask] = np.nan
    
    points = list(zip(x_in,y_in))
    
    uu = griddata(points,u,Xi,method='linear')
    vv = griddata(points,v,Xi,method='linear')
    del points

    u_out[:,:,i] = np.reshape(uu,(ydim,xdim))
    v_out[:,:,i] = np.reshape(vv,(ydim,xdim))
    del uu,vv


xx,yy = np.meshgrid(x_out,y_out)
lon, lat = mf.m2lonlat(lon_0,lat_0,xx,yy,std_lat1=std_lat1,std_lat2=std_lat2)
x_out -= x_out.min()
y_out -= y_out.min()
savemat('rtofs_vel.mat',
        {'x':x_out,
         'y':y_out,
         't':t,
         'u':u_out,
         'v':v_out,
         'lon':lon,
         'lat':lat,
         'lon_0':lon_0,
         'lat_0':lat_0,
         'std_lat1':std_lat1,
         'std_lat2':std_lat2
         })
    
