# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 19:34:10 2016

@author: pnolan86

source: https://www.youtube.com/watch?v=mlAuOKD1ff8

"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from time import gmtime, strftime, strptime
from mpl_toolkits.basemap import Basemap
import calendar
#import scipy.ndimage.filters as filter
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

def km2lonlat(reflon,reflat,x,y ):
    #KM2LONLAT Summary of this function goes here
    #   Inverse Lambert Conformal Projection
    #stdlat1  =  deg2rad(30)
    #stdlat2  =  deg2rad(60)
    stdlat1  =  deg2rad(40)
    stdlat2  =  deg2rad(42)
    R=6371
    reflon = deg2rad(reflon)
    reflat = deg2rad(reflat)
    n = np.log(np.cos(stdlat1)*sec(stdlat2)) / np.log(np.tan(0.25*np.pi+0.5*stdlat2)*cot(0.25*np.pi+0.5*stdlat1))
    F=(np.cos(stdlat1)*(np.tan(0.25*np.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*np.pi+0.5*reflat)**n)
    p = np.sign(n)*np.sqrt(x**2+(p0-y)**2)
    th = np.arctan(x/(p0-y))
    lon = th/n + reflon
    lat = 2*np.arctan((R*F/p)**(1/n))-np.pi/2
    lon = rad2deg(lon)
    lat = rad2deg(lat)
    return lon,lat

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \cdot 10^{{{}}}$'.format(a, b)
    
    #return x

tstart = calendar.timegm(strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

### GENERATE RHODOT
ncfile="windagedata.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t=14
dx=200 #m
dy=200 #m
lat = vars["lat"][:]
lon = vars["lon"][:]
time = 86400*vars["time"][:]+tstart
u = vars["eastward_vel"][t,:,:]
v = vars["northward_vel"][t,:,:]
reflat = 0.5*(max(lat)+min(lat))#midpoint lat
reflon = 0.5*(max(lon)+min(lon))#midpoint lon
ydim = lat.shape[0]
xdim = lon.shape[0]
lon, lat = np.meshgrid(lon,lat)

###MIT SIMULATION
dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)
lon_min = lon.min()
lon_max = lon.max()
lat_min = lat.min()
lat_max = lat.max()

s1 = np.ma.empty([ydim,xdim])
s1Eig = np.ma.empty([ydim,xdim,2])
sn = np.ma.empty([ydim,xdim])
snEig = np.ma.empty([ydim,xdim,2])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors = np.linalg.eig(S)
            idx = eigenValues.argsort()   
            eigenVectors = eigenVectors[:,idx]            
            eigenValues = eigenValues[idx]
            s1[i,j] = 3600*eigenValues[0]
            s1Eig[i,j,:] = eigenVectors[:,0]
            sn[i,j] = 3600*eigenValues[-1]
            snEig[i,j,:] = eigenVectors[:,-1]

        else:
            s1Eig[i,j,0] = np.ma.masked
            s1Eig[i,j,1] = np.ma.masked            
            s1[i,j] = np.ma.masked
            snEig[i,j,0] = np.ma.masked
            snEig[i,j,1] = np.ma.masked            
            sn[i,j] = np.ma.masked

ds1dy,ds1dx = np.gradient(s1,dy,dx,edge_order=2)
ds1dydy,ds1dydx = np.gradient(ds1dy,dy,dx,edge_order=2)
ds1dxdy,ds1dxdx = np.gradient(ds1dx,dy,dx,edge_order=2)

dsndy,dsndx = np.gradient(sn,dy,dx,edge_order=2)
dsndydy,dsndydx = np.gradient(dsndy,dy,dx,edge_order=2)
dsndxdy,dsndxdx = np.gradient(dsndx,dy,dx,edge_order=2)




s1dirdiv = np.ma.empty([ydim,xdim])
s1concav = np.ma.empty([ydim,xdim])
sndirdiv = np.ma.empty([ydim,xdim])
snconcav = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        if (ds1dx[i,j] and ds1dy[i,j] and ds1dxdy[i,j] and ds1dydy[i,j] and ds1dxdx[i,j] and ds1dydx[i,j]) is not np.ma.masked:    
            s1dirdiv[i,j] = np.dot([ds1dx[i,j],ds1dy[i,j]],s1Eig[i,j,:])
            s1concav[i,j] = np.dot(np.dot([[ds1dxdx[i,j],ds1dxdy[i,j]],[ds1dydx[i,j],ds1dydy[i,j]]],s1Eig[i,j,:]),s1Eig[i,j,:])
        #print aconcav[i,j]
        else:
            s1dirdiv[i,j] = np.ma.masked
            s1concav[i,j] = np.ma.masked

        if (dsndx[i,j] and dsndy[i,j] and dsndxdy[i,j] and dsndydy[i,j] and dsndxdx[i,j] and dsndydx[i,j]) is not np.ma.masked:    
            sndirdiv[i,j] = np.dot([dsndx[i,j],dsndy[i,j]],snEig[i,j,:])
            snconcav[i,j] = np.dot(np.dot([[dsndxdx[i,j],dsndxdy[i,j]],[dsndydx[i,j],dsndydy[i,j]]],snEig[i,j,:]),snEig[i,j,:])
        #print aconcav[i,j]
        else:
            sndirdiv[i,j] = np.ma.masked
            snconcav[i,j] = np.ma.masked
thresh=0.9
s1_percentile = -s1.compressed()
percentile = int(np.ceil(thresh*s1_percentile.size))
s1thresh = -np.sort(s1_percentile)[percentile]
s1dirdiv = np.ma.masked_where(s1concav<0,s1dirdiv)
s1dirdiv = np.ma.masked_where(s1>=s1thresh,s1dirdiv)

sn_percentile = sn.compressed()
percentile = int(np.ceil(thresh*sn_percentile.size))
snthresh = np.sort(sn_percentile)[percentile]
sndirdiv = np.ma.masked_where(snconcav>0,sndirdiv)
sndirdiv = np.ma.masked_where(sn<=snthresh,sndirdiv)

###GENERATE MAP
plt.close('all')
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )
parallels = np.arange(41.1,lat_max+0.1,0.1)
meridians = np.arange(-70.2,lon_min-0.1,-0.1)

#s1=3600*s1

"""
#PLOT iLCS segments to see which ones are worth investigating

aridge = m.contour(lon,lat,s1dirdiv,levels=0,latlon=True,colors='blue')
#rridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
plt.close('all')
pp = aridge.collections[0].get_paths()
index = 0
for p in pp:#[13,17,40,46,47,55]:
    v = p.vertices
    if np.size(v[:,0])>=5:
        x = v[:,0]
        y = v[:,1]
        m.contour(lon,lat,s1dirdiv,levels=0,latlon=True,colors='blue')
        m.plot(x,y,'red')
        m.drawcoastlines()
        
        plt.savefig('a{0:04d}'.format(index), transparent=True, bbox_inches='tight')
        plt.close('all')
    index += 1

plt.close('all')
rridge = m.contour(lon,lat,sndirdiv,levels=0,latlon=True,colors='red')
pp = rridge.collections[0].get_paths()
index = 0
for p in pp:#[13,17,40,46,47,55]:
    v = p.vertices
    if np.size(v[:,0])>=5:
        x = v[:,0]
        y = v[:,1]
        m.contour(lon,lat,sndirdiv,levels=0,latlon=True,colors='blue')
        m.plot(x,y,'red')
        m.drawcoastlines()
        
        plt.savefig('r{0:04d}'.format(index), transparent=True, bbox_inches='tight')
        plt.close('all')
    index += 1
#"""



"""
#PLOT iLCS segments to see which ones are worth investigating
plt.figure()
aridge = m.contour(lon,lat,s1dirdiv,levels=0,latlon=True,colors='blue')
#rridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
pp = aridge.collections[0].get_paths()
index = 0
m.contour(lon,lat,s1dirdiv,levels=0,latlon=True,colors='blue')
for p in [13,15,17,26,42,45,49,50]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    m.plot(x,y,'red')    
    index += 1
m.drawcoastlines()

plt.figure()
rridge = m.contour(lon,lat,sndirdiv,levels=0,latlon=True,colors='red')
pp = rridge.collections[0].get_paths()
index = 0
m.contour(lon,lat,sndirdiv,levels=0,latlon=True,colors='blue')
for p in [35,37,40,45,48,56,81,99,104,113]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    m.plot(x,y,'red')
    index += 1
m.drawcoastlines()

#"""

"""
#Output for advection
import csv
aridge = m.contour(lon,lat,s1dirdiv,levels=0,latlon=True,colors='blue')
pp = aridge.collections[0].get_paths()
index = 0
for p in [13,15,17,26,42,45,49,50]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    x,y = m(x,y,inverse=True)
    with open('attracting_ridge_{:02d}.txt'.format(p), 'w', newline='') as f:
        writer = csv.writer(f)
        for row in zip(x,y):
            writer.writerow(row)
        f.close()

rridge = m.contour(lon,lat,sndirdiv,levels=0,latlon=True,colors='red')
pp = rridge.collections[0].get_paths()
for p in [35,37,40,45,48,56,81,99,104,113]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    x,y = m(x,y,inverse=True)
    with open('repelling_ridge_{:02d}.txt'.format(p), 'w', newline='') as f:
        writer = csv.writer(f)
        for row in zip(x,y):
            writer.writerow(row)
        f.close()
#"""

""
#TRACERS
origin_x = [46787.5,55175.7,4028.29,7813.2,14462.4,16661.7,42491.1]
origin_y = [35571,45493.5,36133.6,40327.7,37719.1,36696.2,32706.7]
theta = np.linspace(0,2*np.pi,51)
radius = 750#m
x=[]
y=[]
for i in range(len(origin_x)):
    x.extend(radius*np.cos(theta)+origin_x[i])
    y.extend(radius*np.sin(theta)+origin_y[i])
    
lon_t, lat_t = m(x,y,inverse=True)
import csv
with open('tracers.txt', 'w', newline='') as f:
    writer = csv.writer(f)
    for row in zip(lon_t, lat_t):
        writer.writerow(row)
    f.close()


#"""
 
""

###GENERATE PLOTS
pltsize = [15,12]
#pltsize = [5,4]

fig = plt.figure(figsize=pltsize, dpi=150)
ax = plt.subplot(111)
#geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.25*colorlevel,vmax=0.25*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
geatmap = m.contourf(lon,lat,s1,levels=np.linspace(s1.min(),s1.max(),301),latlon=True)#,shading='gourand')
plt.colorbar()
aridge = m.contour(lon,lat,s1dirdiv,levels =[0],colors='blue',latlon=True,alpha=0.6)
#rridge = m.contour(flon,flat,rdirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
#plt.savefig('output1_v2.3.eps', transparent=True, bbox_inches='tight')

fig = plt.figure(figsize=pltsize, dpi=150)
ax = plt.subplot(111)
#geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.25*colorlevel,vmax=0.25*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
geatmap = m.contourf(lon,lat,sn,levels=np.linspace(sn.min(),sn.max(),301),latlon=True)#,shading='gourand')
plt.colorbar()
aridge = m.contour(lon,lat,sndirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
#rridge = m.contour(flon,flat,rdirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
#plt.savefig('output1_v2.3.eps', transparent=True, bbox_inches='tight')

fig = plt.figure(figsize=pltsize, dpi=150)
ax = plt.subplot(111)
#geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.25*colorlevel,vmax=0.25*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
aridge = m.contour(lon,lat,s1dirdiv,levels =[0],colors='blue',latlon=True)
rridge = m.contour(lon,lat,sndirdiv,levels =[0],colors='red',latlon=True)
m.scatter(x,y)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
#plt.savefig('output1_v2.3.eps', transparent=True, bbox_inches='tight')

#"""