# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 22:16:30 2018

@author: pnola
"""

import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
from mpl_toolkits.basemap import Basemap
import time
import calendar

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
ftle = vars['FTLE'][::6,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
time = vars['time'][:]
root.close()

dims = ftle.shape

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
lon, lat = np.meshgrid(lon,lat)
for t in range(dims[0]):
    plt.figure()
    #lon,lat = np.meshgrid(lon,lat)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True)
    cs=m.pcolormesh(lon,lat,ftle[t,:,:],latlon=True)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    plt.colorbar()
    plt.savefig('ftle_{:}'.format(t))
    plt.close('all')


import winsound
duration = 500  # millisecond
freq = 700  # Hz
winsound.Beep(freq, duration)
