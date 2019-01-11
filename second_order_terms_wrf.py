# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 11:28:01 2018

@author: pnola
def kevrekidis(t,Y, eps=0.01):
  return np.array([-Y[0]-Y[1]+2, 1/eps*(Y[0]**3-Y[1])])

def ex11(t,Y):
  return -np.array([np.tanh(Y[0]**2/4)+Y[0],Y[0]+2*Y[1]])

def rotHoop(t,Y, eps=0.1, gamma=2.3):
  return np.array([Y[1], 1/eps*(np.sin(Y[0])*(gamma*np.cos(Y[0])-1)-Y[1])])


def verhulst(t,Y, eps=0.01):
  return [1, 1/eps*(Y[0]*Y[1]-Y[1]**2)]


"""

import numpy as np
import matplotlib.pyplot as plt
#from velocities import van_der_pol_oscillator as vel_func
from velocities import double_gyre as vel_func
#from velocities import auto_double_gyre as vel_func
from velocities import co

from netCDF4 import Dataset
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time
import calendar
plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}


height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

#xdim = 405
#ydim = 325
time_step = 21 #2100hrs

tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ncfile="wrf_2011_07_01"#"ftle_80m.nc"
root = Dataset(ncfile,'r')
vars = root.variables
u = mf.unstagger(vars['U'][:,height_level,:,:],axis=2)
v = mf.unstagger(vars['V'][:,height_level,:,:],axis=1)
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()
[tdim,ydim,xdim]=u.shape



#u=u[time_step,:,:]
#v=v[time_step,:,:]


#lon = lon[-1,:,:]
#lat = lat[-1,:,:]
#latin = lat[:25,:,:]
#longin = lon[:25,:,:]
dt=3600 #sec
#x = np.linspace(-grid_spacing*(xdim-1)/2,grid_spacing*(xdim-1)/2,xdim)
x = np.linspace(0,grid_spacing*(xdim-1),xdim)
dx = x[1]-x[0]
#y = np.linspace(-grid_spacing*(ydim-1)/2,grid_spacing*(ydim-1)/2,ydim)
y = np.linspace(0,grid_spacing*(ydim-1),ydim)
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)

dudt = np.gradient(u,dt,axis=0)
dvdt = np.gradient(v,dt,axis=0)

dudy,dudx = np.gradient(u,dy,dx,axis=(1,2))
dvdy,dvdx = np.gradient(v,dy,dx,axis=(1,2))

Du = dudt+u*dudx+v*dudy
Dv = dvdy+u*dvdx+v*dvdy

dDudy,dDudx = np.gradient(Du,dy,dx,axis=(1,2))
dDvdy,dDvdx = np.gradient(Dv,dy,dx,axis=(1,2))


s1 = np.ma.empty([tdim,ydim,xdim])
b1 = np.ma.empty([tdim,ydim,xdim])
s2 = np.ma.empty([tdim,ydim,xdim])
b2 = np.ma.empty([tdim,ydim,xdim])
for t in range(tdim):
    for i in range(ydim):
        for j in range(xdim):
            if (dDudx[t,i,j] and dDudy[t,i,j] and dDvdx[t,i,j] and dDvdy[t,i,j]) is not np.ma.masked:    
                Grad_v = np.array([[dudx[t,i,j], dudy[t,i,j]], [dvdx[t,i,j], dvdy[t,i,j]]])
                Grad_D = np.array([[dDudx[t,i,j], dDudy[t,i,j]], [dDvdx[t,i,j], dDvdy[t,i,j]]])
                S = 0.5*(Grad_v + np.transpose(Grad_v))
                B = (0.5*(Grad_D + np.transpose(Grad_D)+np.matmul(np.transpose(Grad_v),Grad_v)))
                s1[t,i,j] = np.min(np.linalg.eig(S)[0])
                b1[t,i,j] = np.min(np.linalg.eig(B)[0])
                s2[t,i,j] = np.max(np.linalg.eig(S)[0])
                b2[t,i,j] = np.max(np.linalg.eig(B)[0])
            else:
                s1[t,i,j] = np.ma.masked
                b1[t,i,j] = np.ma.masked
                s2[t,i,j] = np.ma.masked
                b2[t,i,j] = np.ma.masked

s1=-3600*s1[0,:,:]
b1=-3600*b1[0,:,:]
s2=3600*s2[0,:,:]
b2=3600*b2[0,:,:]


lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'h',
            area_thresh=1000.,
            )

#lon,lat = np.meshgrid(lon,lat)
plt.subplot(221)
cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title('$-s_{1}$')
plt.colorbar()

plt.subplot(222)
cs=m.contourf(lon,lat,b1,levels=np.linspace(b1.min(axis=None),b1.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title('$-b_{1}$')
plt.colorbar()

plt.subplot(223)
cs=m.contourf(lon,lat,s2,levels=np.linspace(s2.min(axis=None),s2.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title('$s_{2}$')
plt.colorbar()

plt.subplot(224)
cs=m.contourf(lon,lat,b2,levels=np.linspace(b2.min(axis=None),b2.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title('$b_{2}$')
plt.colorbar()







#t=1
#hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
#plt.savefig('t.png', transparent=False, bbox_inches='tight')

'''

plt.close('all')
dimx = 301
#dimy = int(np.ceil(dimx/2.5))
dimy = int(np.ceil(dimx/2))
#dimy = 0.5*dimx
t0 = 0
tf =-0.1#days
dimt =2
t = np.linspace(t0,tf,dimt)
dt = t[1]-t[0]
x = np.linspace(0,2,dimx)
y = np.linspace(0,1,dimy)
dx = x[1]-x[0]
dy = y[1]-y[0]

#[x,y,t]=np.meshgrid(x,y,t)
[x,y]=np.meshgrid(x,y)

dudy = np.empty([dimy,dimx,dimt])
dudx = np.empty([dimy,dimx,dimt])
dvdy = np.empty([dimy,dimx,dimt])
dvdx = np.empty([dimy,dimx,dimt])
for tt in range(dimt):
    u,v = vel_func(t[tt],[x,y])
    dudy[:,:,tt],dudx[:,:,tt] = np.gradient(u,dy,dx)#,axis=(0,1))
    dvdy[:,:,tt],dvdx[:,:,tt] = np.gradient(v,dy,dx)#,axis=(0,1))


dudydt = np.gradient(dudy,dt,axis=2)
dudxdt = np.gradient(dudx,dt,axis=2)
dvdydt = np.gradient(dvdy,dt,axis=2)
dvdxdt = np.gradient(dvdx,dt,axis=2)

s1 = np.ma.empty([dimy,dimx,dimt])
b1 = np.ma.empty([dimy,dimx,dimt])
s2 = np.ma.empty([dimy,dimx,dimt])
b2 = np.ma.empty([dimy,dimx,dimt])
for t in range(dimt):
    for i in range(dimy):
        for j in range(dimx):
            if (dudxdt[i,j,t] and dudydt[i,j,t] and dvdxdt[i,j,t] and dvdydt[i,j,t]) is not np.ma.masked:    
                Grad = np.array([[dudx[i,j,t], dudy[i,j,t]], [dvdx[i,j,t], dvdy[i,j,t]]])
                Grad_dt = np.array([[dudxdt[i,j,t], dudydt[i,j,t]], [dvdxdt[i,j,t], dvdydt[i,j,t]]])
                S = 0.5*(Grad + np.transpose(Grad))
                D = (0.5*(Grad_dt + np.transpose(Grad_dt)+np.matmul(np.transpose(Grad),Grad)))
                s1[i,j,t] = np.min(np.linalg.eig(S)[0])
                b1[i,j,t] = np.min(np.linalg.eig(D)[0])
                s2[i,j,t] = np.max(np.linalg.eig(S)[0])
                b2[i,j,t] = np.max(np.linalg.eig(D)[0])
            else:
                s1[i,j,t] = np.ma.masked
                b1[i,j,t] = np.ma.masked
                s2[i,j,t] = np.ma.masked
                b2[i,j,t] = np.ma.masked
    
plt.register_cmap(name='co', data=co())


s1=s1[:,:,0]
b1=b1[:,:,0]
s2=s2[:,:,0]
b2=b2[:,:,0]
y=y[:,:]
x=x[:,:]


plt.figure()

plt.subplot(221)
#plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
plt.pcolormesh(x,y,-3600*s1)
plt.colorbar()
plt.title('s$_{1}$')

plt.subplot(222)
#plt.contourf(x,y,sigma,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.pcolormesh(x,y,-3600*b1)#,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('b$_{1}$')

plt.subplot(223)
#plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
plt.pcolormesh(x,y,3600*s2)
plt.colorbar()
plt.title('s$_{2}$')

plt.subplot(224)
#plt.contourf(x,y,sigma,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.pcolormesh(x,y,3600*b2)#,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('b$_{2}$')

#plt.gca().set_aspect('equal', adjustable='box', anchor='C')



for i in range(ydim):
    for j in range(xdim):
        if (dudx[t,i,j] and dudy[t,i,j] and dvdx[t,i,j] and dvdy[t,i,j]) is not np.ma.masked:    
            Grad = np.array([[dudx[t,i, j], dudy[t,i, j]], [dvdx[t,i, j], dvdy[t,i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = -3600*np.min(np.linalg.eig(S)[0])

        else:
            s1[i,j] = np.ma.masked
            

'''
#"""