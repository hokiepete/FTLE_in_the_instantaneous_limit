import mapping_functions as mf
#from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
plt.close('all')
#25N
#-90E


root = Dataset('data/07_25_2019/rtofs_glo_2ds_f000_1hrly_prog.nc','r')
vars = root.variables #dictionary, all variables in dataset\
lon = vars['Longitude'][:-1,:]-360
lat = vars['Latitude'][:-1,:]
u = vars['u_velocity'][:,:,:-1,:].squeeze()
v = vars['v_velocity'][:,:,:-1,:].squeeze()
del vars
root.close()

x,y = mf.lonlat2km(-90,25,lon,lat,std_lat1=20,std_lat2=30)
del lon, lat


plt.figure()
plt.pcolormesh(x,y,np.sqrt(u**2+v**2),vmin=0,vmax=3)#np.sqrt(u**2+v**2))
plt.xlim([-900,1100])
plt.ylim([-800,800])
plt.colorbar()

u = u[x>=-1000]
v = v[x>=-1000]
y = y[x>=-1000]
x = x[x>=-1000]

u = u[x<=1200]
v = v[x<=1200]
y = y[x<=1200]
x = x[x<=1200]

u = u[y>=-900]
v = v[y>=-900]
x = x[y>=-900]
y = y[y>=-900]

u = u[y<=900]
v = v[y<=900]
x = x[y<=900]
y = y[y<=900]

v[u.mask] = np.nan
u[u.mask] = np.nan

u[v.mask] = np.nan
v[v.mask] = np.nan

points = list(zip(x,y))
xdim=400
ydim=int(np.ceil(xdim*4/5))
xx = np.linspace(-900,1100,xdim)
dx = xx[1]-xx[0]
yy = np.linspace(-800,800,ydim)
dy = yy[1]-yy[0]
xx,yy = np.meshgrid(xx,yy)

Xi= list(zip(xx.ravel(),yy.ravel()))

uu = griddata(points,u,Xi)
vv = griddata(points,v,Xi)
del points, Xi
x = np.reshape(xx,(ydim,xdim))
y = np.reshape(yy,(ydim,xdim))
u = np.reshape(uu,(ydim,xdim))
v = np.reshape(vv,(ydim,xdim))
del uu,vv,xx,yy

plt.figure()
plt.pcolormesh(x,y,np.sqrt(u**2+v**2),vmin=0,vmax=3)#np.sqrt(u**2+v**2))
plt.xlim([-900,1100])
plt.ylim([-800,800])
plt.colorbar()


dx *=1000
dy *=1000
dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)
s1 = np.empty((ydim,xdim))
for i in range(ydim):
    for j in range(xdim):
        if ~np.isnan(dudy[i,j]) and ~np.isnan(dudx[i,j]) and ~np.isnan(dvdy[i,j]) and ~np.isnan(dvdx[i,j]):
            Grad_v = np.array([[dudx[i,j],dudy[i,j]],[dvdx[i,j],dvdy[i,j]]])
            S = 0.5*(Grad_v + Grad_v.T)
            s1[i,j] = np.min(np.linalg.eig(S)[0])
            
        else:
            s1[i,j]=np.nan
            
            
        
plt.figure()
plt.pcolormesh(x,y,-s1)#np.sqrt(u**2+v**2))
plt.xlim([-900,1100])
plt.ylim([-800,800])
plt.colorbar()
#plt.scatter(0,0)
