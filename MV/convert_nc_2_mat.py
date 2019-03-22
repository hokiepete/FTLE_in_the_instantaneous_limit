from netCDF4 import Dataset
import numpy as np
from hdf5storage import savemat
import time
import calendar
siavashtime = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

def cot(th):
    return 1.0/np.tan(th)

def sec(th):
    return 1.0/np.cos(th)

def deg2rad(deg):
    return deg*np.pi/180.0

def rad2deg(rad):
    return rad*180.0/np.pi

def lonlat2km(reflon,reflat,lon,lat ):
    #LONLAT2KM Summary of this function goes here
    #   Uses Lambert Conformal Projection
    #stdlat1  =  deg2rad(30)
    #stdlat2  =  deg2rad(60)
    stdlat1  =  deg2rad(40)
    stdlat2  =  deg2rad(42)
    R=6371
    reflon = deg2rad(reflon)
    reflat = deg2rad(reflat)
    lon = deg2rad(lon)
    lat = deg2rad(lat)
    n = np.log(np.cos(stdlat1)*sec(stdlat2)) / np.log(np.tan(0.25*np.pi+0.5*stdlat2)*cot(0.25*np.pi+0.5*stdlat1))
    F=(np.cos(stdlat1)*(np.tan(0.25*np.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*np.pi+0.5*reflat)**n)
    p = R*F*(cot(0.25*np.pi+0.5*lat)**n)
    th = n*(lon-reflon)
    x=p*np.sin(th)
    y=p0-p*np.cos(th)
    return x,y
ncfile = "windagedata.nc"

root = Dataset(ncfile,'r')
vars = root.variables
lat = vars["lat"][:]
lon = vars["lon"][:]
u = vars["eastward_vel"][:].transpose([1,2,0])
v = vars["northward_vel"][:].transpose([1,2,0]) 
latorg = 0.5*(lat[0]+lat[-1])
lonorg = 0.5*(lon[0]+lon[-1])
time = 3600*24*(vars['time'][:]-74)
root.close()

yy = np.linspace(0,2000*23.9,240)
xx = np.linspace(0,2000*26.2,263)
#xx,yy = np.meshgrid(xx,yy)

savemat('wrf_vel_data.mat',{'u':u.data,'v':v.data,'x':xx,'y':yy,'time':time.data},format='7.3')

