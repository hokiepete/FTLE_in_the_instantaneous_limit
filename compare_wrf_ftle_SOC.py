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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
from mpl_toolkits.basemap import Basemap
import time
import calendar
import matplotlib
matplotlib.use('Agg')
plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}


height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
#grid_spacing = 12*1000 #km 2 m
grid_spacing = 12 #km
#xdim = 405
#ydim = 325
time_step = 24 #2100hrs

figwidth = 6

tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ncfile="wrf_2011_07_01"#"ftle_80m.nc"
root = Dataset(ncfile,'r')
vars = root.variables
#u = mf.unstagger(vars['U'][:,height_level,:,:],axis=2)
#v = mf.unstagger(vars['V'][:,height_level,:,:],axis=1)
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()

ncfile="wrf_2011_07_02"#"ftle_80m.nc"
root = Dataset(ncfile,'r')
vars = root.variables
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
root.close()


u = 24*3.6*mf.unstagger(u,axis=2)
v = 24*3.6*mf.unstagger(v,axis=1)

[tdim,ydim,xdim]=u.shape

dt=1/24#3600 #sec
x = np.linspace(0,grid_spacing*(xdim-1),xdim)
dx = x[1]-x[0]
y = np.linspace(0,grid_spacing*(ydim-1),ydim)
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)


dudy,dudx = np.gradient(u,dy,dx,axis=(1,2))
dvdy,dvdx = np.gradient(v,dy,dx,axis=(1,2))

dudydt = np.gradient(dudy,dt,axis=0)
dudxdt = np.gradient(dudx,dt,axis=0)
dvdydt = np.gradient(dvdy,dt,axis=0)
dvdxdt = np.gradient(dvdx,dt,axis=0)

u=u[time_step,:,:]
v=v[time_step,:,:]
dudy=dudy[time_step,:,:]
dudx=dudx[time_step,:,:]
dvdy=dvdy[time_step,:,:]
dvdx=dvdx[time_step,:,:]
dudydt=dudydt[time_step,:,:]
dudxdt=dudxdt[time_step,:,:]
dvdydt=dvdydt[time_step,:,:]
dvdxdt=dvdxdt[time_step,:,:]


[ydim,xdim]=u.shape
#FigSize=(2*figwidth, int(np.ceil(ydim/xdim*figwidth)))


s1 = np.ma.empty([ydim,xdim])
s2 = np.ma.empty([ydim,xdim])
corr1 = np.ma.empty([ydim,xdim])
corr2 = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        if (dudxdt[i,j] and dudydt[i,j] and dvdxdt[i,j] and dvdydt[i,j]) is not np.ma.masked:    
            Grad = np.array([[dudx[i,j], dudy[i,j]], [dvdx[i,j], dvdy[i,j]]])
            Grad_dt = np.array([[dudxdt[i,j], dudydt[i,j]], [dvdxdt[i,j], dvdydt[i,j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            B = (0.5*(Grad_dt + np.transpose(Grad_dt)+np.dot(np.transpose(Grad),Grad)))
            eigenValues, eigenVectors_temp = np.linalg.eig(S)
            idx = eigenValues.argsort()
            eigenMin = eigenVectors_temp[:,idx[0]]
            eigenMax = eigenVectors_temp[:,idx[1]]
            s1[i,j] = eigenValues[idx[0]]
            s2[i,j] = eigenValues[idx[1]]
            corr1[i,j] = -s1[i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin)
            corr2[i,j] = -s2[i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax)
            #sigb[i,j] = -s1[i,j] - (tf-t0)*(-s1[t,i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin))
            #sigf[i,j] = s2[i,j] - (tf-t0)*(-s2[t,i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax))
        else:
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked
            corr1[i,j] = np.ma.masked
            corr2[i,j] = np.ma.masked



lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'c',
            area_thresh=1000.,
            )

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
flat = vars["lat"][:]#.reshape([ydim,xdim])
flon = vars["lon"][:]#.reshape([ydim,xdim])
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
time = vars['time'][:]
flon, flat = np.meshgrid(flon,flat)
root.close()

for i,k in enumerate(time):
    print(i,k)
    fig=plt.figure(figsize=(12,5))
    #fig, ax = plt.subplots()
    plt.subplot(121)
    #k=-k
    sig=-s1-k*corr1
    #cs=m.pcolormesh(lon,lat,sig,latlon=True)
    cs=m.contourf(lon,lat,sig,levels=np.linspace(sig.min(axis=None),sig.max(axis=None),301),latlon=True)
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    #plt.colorbar()
    plt.colorbar( fraction=0.04, pad=0.04)
    plt.title('$s_{1}-T*\\left(-s_{1}^{2}+0.5 \\langle \\langle \\xi_{s_{1}}, B \\rangle, \\xi_{s_{1}}\\rangle\\right)$'+' T = {0:1.3f} Days'.format(k))
    
    f=ftle[-1-i,:,:]
    plt.subplot(122)
    if i == 0: 
        cs=m.pcolormesh(flon,flat,f,latlon=True)
    else:
        cs=m.contourf(flon,flat,f,levels=np.linspace(f.min(axis=None),f.max(axis=None),301),latlon=True)
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    #plt.colorbar()
    plt.colorbar( fraction=0.04, pad=0.04)
    plt.title('$\\sigma$')
    
    plt.savefig('ftle_s1+HOT_comp_{:}.png'.format(i))
    plt.close('all')
#plt.axis('tight')