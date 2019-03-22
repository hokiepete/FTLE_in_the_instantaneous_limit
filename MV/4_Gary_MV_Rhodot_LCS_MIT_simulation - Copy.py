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
import matplotlib.ticker as ticker
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

cdict = {'red':  [(0.0, 0.0000, 0.0000),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 1.0000, 1.0000)],
        'green': [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.5450, 0.5450)],
        'blue':  [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.0000, 0.0000)]}
plt.register_cmap(name='CO', data=cdict)

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

#s1 = np.ma.empty(dim)
rhodot = np.ma.empty([ydim,xdim])
nudot = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            #s1[i,j] = np.min(np.linalg.eig(S)[0])
            rhodot[i, j] = np.dot(Utemp, np.dot(np.dot(np.transpose(J), np.dot(S, J)), Utemp))/np.dot(Utemp, Utemp)
            nudot[i, j] = np.dot(Utemp,np.dot(np.trace(S)*np.identity(2) - 2*S, Utemp))/np.dot(Utemp, Utemp)

        else:
            rhodot[i,j] = np.ma.masked
            nudot[i,j] = np.ma.masked
            #s1[i,j] = np.ma.masked

rhodot = 3600*rhodot
#s1 = 3600*s1
colorlevel = np.max(np.fabs([rhodot.min(),rhodot.max()]))

### CALCULATE LCS
athresh=0.38
rthresh=0.38#40
#rthresh=0.45
aftlesnap=-3 #-2hr
rftlesnap= 2 #+2hr
ncfile = 'MV_FTLE_-6hrs_NoWindage.nc'
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
aftle = vars['FTLE'][:,aftlesnap,0]
flon = vars['initial_lon'][:]
flat = vars['initial_lat'][:]
ydim = 479
xdim = 525
acg = vars['CauchyGreen'][:,aftlesnap,:,:]
aftle = np.reshape(aftle,[ydim,xdim])/24
flon = np.reshape(flon,[ydim,xdim])
flat = np.reshape(flat,[ydim,xdim])
acg = np.reshape(acg,[ydim,xdim,2,2])
root.close()
ncfile = 'MV_FTLE_+6hrs_NoWindage.nc'

root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
rftle = vars['FTLE'][:,rftlesnap,0]
rcg = vars['CauchyGreen'][:,rftlesnap,:,:]
rftle = np.reshape(rftle,[ydim,xdim])/24
rcg = np.reshape(rcg,[ydim,xdim,2,2])
root.close()

adfdy,adfdx = np.gradient(aftle,dy/2.0,dx/2.0,edge_order=2)
adfdydy,adfdydx = np.gradient(adfdy,dy/2.0,dx/2.0,edge_order=2)
adfdxdy,adfdxdx = np.gradient(adfdx,dy/2.0,dx/2.0,edge_order=2)
rdfdy,rdfdx = np.gradient(rftle,dy/2.0,dx/2.0,edge_order=2)
rdfdydy,rdfdydx = np.gradient(rdfdy,dy/2.0,dx/2.0,edge_order=2)
rdfdxdy,rdfdxdx = np.gradient(rdfdx,dy/2.0,dx/2.0,edge_order=2)

adirdiv = np.ma.empty([ydim,xdim])
rdirdiv = np.ma.empty([ydim,xdim])
aconcav = np.ma.empty([ydim,xdim])
rconcav = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        if (adfdx[i,j] and adfdy[i,j] and adfdxdy[i,j] and adfdydy[i,j] and adfdxdx[i,j] and adfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(acg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            adirdiv[i,j] = np.dot([adfdx[i,j],adfdy[i,j]],eigenVectors[:,0])
            aconcav[i,j] = np.dot(np.dot([[adfdxdx[i,j],adfdxdy[i,j]],[adfdydx[i,j],adfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        #print aconcav[i,j]
        else:
            adirdiv[i,j] = np.ma.masked
            aconcav[i,j] = np.ma.masked

        if (rdfdx[i,j] and rdfdy[i,j] and rdfdxdy[i,j] and rdfdydy[i,j] and rdfdxdx[i,j] and rdfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(rcg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            rdirdiv[i,j] = np.dot([rdfdx[i,j],rdfdy[i,j]],eigenVectors[:,0])
            rconcav[i,j] = np.dot(np.dot([[rdfdxdx[i,j],rdfdxdy[i,j]],[rdfdydx[i,j],rdfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        else:
           rdirdiv[i,j] = np.ma.masked
           rconcav[i,j] = np.ma.masked



adirdiv = np.ma.masked_where(aconcav>0,adirdiv)
adirdiv = np.ma.masked_where(aftle<=athresh,adirdiv)
rdirdiv = np.ma.masked_where(rconcav>0,rdirdiv)
rdirdiv = np.ma.masked_where(rftle<=rthresh,rdirdiv)

#rftle = np.ma.masked_where(rftle<=rthresh,rftle)
#aftle = np.ma.masked_where(aftle<=athresh,aftle)
rftle = np.ma.masked_where(np.ma.getmask(rdirdiv),rftle)
aftle = np.ma.masked_where(np.ma.getmask(adirdiv),aftle)

#rhodot = np.ma.masked_where(nudot < 0,rhodot)
rdim = rhodot.shape
for i in range(rdim[0]):
    for j in range(rdim[1]):
        if (rhodot[i,j]>0 and nudot[i,j]>0) or (rhodot[i,j]<0 and nudot[i,j]<0):
            'nothing'
            #rhodot = rhodot[(rhodot>0 and nudot>0) or (rhodot<0 and nudot<0)]
        else:
            rhodot[i,j] = np.ma.masked

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


###GENERATE PLOTS
pltsize = [15,12]
#pltsize = [5,4]
'''
fig = plt.figure(1,figsize=pltsize, dpi=150)
ax = plt.subplot(221)
heatmap = m.contourf(flon,flat,rftle,levels=np.linspace(rftle.min(),rftle.max(),301),latlon=True)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax.set_aspect('equal', adjustable='box', anchor='C')
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
cbar = plt.colorbar(heatmap,format=ticker.FuncFormatter(fmt)) 
plt.title(strftime("+2hr FTLE @ %a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))

ax = plt.subplot(222)
heatmap = m.contourf(flon,flat,aftle,levels=np.linspace(aftle.min(),aftle.max(),301),latlon=True)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax.set_aspect('equal', adjustable='box', anchor='C')
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
cbar = plt.colorbar(heatmap,format=ticker.FuncFormatter(fmt)) 
plt.title(strftime("-2hr FTLE @ %a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
'''

fig = plt.figure(2,figsize=pltsize, dpi=150)
ax = plt.subplot(111)
#geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.25*colorlevel,vmax=0.25*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.5,vmax=0.5,cmap = 'CO',latlon=True)#,shading='gourand')
plt.colorbar()
aridge = m.contour(flon,flat,adirdiv,levels =[0],colors='blue',latlon=True,alpha=0.6)
rridge = m.contour(flon,flat,rdirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
plt.savefig('output1_v2.3.eps', transparent=True, bbox_inches='tight')

'''
fig = plt.figure(3,figsize=pltsize, dpi=150)
#m.scatter(lon[153,:],lat[153,:],latlon=True)
plt.plot(lon[153,:],rhodot[153,:],color='orange')
plt.plot(flon[306,:],aftle[306,:],color='blue')
plt.plot(flon[306,:],rftle[306,:],color='red')
#plt.plot(flon[306,:],adirdiv[306,:],color='cyan')
#plt.plot(flon[306,:],rdirdiv[306,:],color='magenta')

#plt.axhline(athresh,color='cyan',linestyle='dashed')
#plt.axhline(rthresh,color='magenta',linestyle='dashed')
'''

fig = plt.figure(4,figsize=pltsize, dpi=150)
bx = plt.subplot(111)
#geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.25*colorlevel,vmax=0.25*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
geatmap = m.contourf(lon,lat,rhodot,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.5,vmax=0.5,cmap = 'CO',latlon=True)#,shading='gourand')
plt.colorbar()
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#uiv = bx.quiver(xx[::2*iterplevl,::2*iterplevl],yy[::2*iterplevl,::2*iterplevl],uu[::2*iterplevl,::2*iterplevl],vv[::2*iterplevl,::2*iterplevl])
m.streamplot(lon,lat,u,v,linewidth=0.4,density=1.75,latlon=True,color='dimgrey')#(0.74, 0.75, 0.75))
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
plt.savefig('output2_v2.3.eps', transparent=True, bbox_inches='tight')
#uiv = bx.quiver(xx,yy,uu,vv)
#dbar = plt.colorbar(geatmap,format=ticker.FuncFormatter(fmt))
#'''
'''
##PLOT DIFFERENT FTLEs
#from mpl_toolkits.basemap import Basemap
ncfile="windagedata.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print vars.keys()
t=14
iterplevl = 24
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
import datetime
print(
    datetime.datetime.fromtimestamp(
        time[t]
    ).strftime('%Y-%m-%d %H:%M:%S')
)
###MIT SIMULATION
speed = np.sqrt(u**2+v**2)
ymin = -dy*(ydim-1.0)/2.0
ymax = dy*(ydim-1)/2.0
xmin = -dx*(xdim-1)/2.0
xmax = dx*(xdim-1)/2.0
x = np.linspace(xmin,xmax,xdim)
y = np.linspace(ymin,ymax,ydim)
x,y = np.meshgrid(x,y)
uu=u
vv=v
dudy,dudx = np.gradient(uu,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(vv,dy,dx,edge_order=2)
dim = uu.shape
lon_min = lon.min()
lon_max = lon.max()
lat_min = lat.min()
lat_max = lat.max()

s1 = np.ma.empty(dim)
A = np.ma.empty(dim)
J = np.array([[0, 1], [-1, 0]])
for i in range(dim[0]):
    for j in range(dim[1]):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and uu[i,j] and vv[i,j]) is not np.ma.masked:    
            Utemp = np.array([uu[i, j], vv[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = np.min(np.linalg.eig(S)[0])
            A[i, j] = np.dot(Utemp, np.dot(np.dot(np.transpose(J), np.dot(S, J)), Utemp))/np.dot(Utemp, Utemp)
        else:
            A[i,j] = np.ma.masked
            s1[i,j] = np.ma.masked

A = 3600*A
s1 = 3600*s1
colormin = A.min()
colormax = A.max()
colorlevel = np.max(np.fabs([colormin,colormax]))

plt.close('all')
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'c',
            area_thresh=0.,
            )
parallels = np.arange(41.1,lat_max+0.1,0.1)
meridians = np.arange(-70.2,lon_min-0.1,-0.1)
pltsize = [15,12]
pltsize = [5,4]
fig = plt.figure(1,figsize=pltsize, dpi=150)
ax = plt.subplot(221)
heatmap = m.contourf(lon,lat,speed,levels=np.linspace(0,1,301),cmap = 'plasma',latlon=True)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax.set_aspect('equal', adjustable='box', anchor='C')

cbar = plt.colorbar(heatmap,format=ticker.FuncFormatter(fmt))
plt.title(strftime("%a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))

bx = plt.subplot(222)
geatmap = m.contourf(lon,lat,A,levels=np.linspace(-colorlevel,colorlevel,301),vmin=-0.5*colorlevel,vmax=0.5*colorlevel,cmap = 'CO',latlon=True)#,shading='gourand')
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#uiv = bx.quiver(xx[::2*iterplevl,::2*iterplevl],yy[::2*iterplevl,::2*iterplevl],uu[::2*iterplevl,::2*iterplevl],vv[::2*iterplevl,::2*iterplevl])
m.streamplot(lon,lat,uu,vv,linewidth=0.4,density=3,latlon=True)
bx.set_aspect('equal', adjustable='box', anchor='C')
#uiv = bx.quiver(xx,yy,uu,vv)
dbar = plt.colorbar(geatmap,format=ticker.FuncFormatter(fmt))

cx = plt.subplot(223)
featmap = m.contourf(lon,lat,s1,levels=np.linspace(s1.min(),s1.max(),301),latlon=True)#,shading='gourand')
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
cx.set_aspect('equal', adjustable='box', anchor='C')
#wiv = cx.quiver(lon[::2*iterplevl,::2*iterplevl],lat[::2*iterplevl,::2*iterplevl],uu[::2*iterplevl,::2*iterplevl],vv[::2*iterplevl,::2*iterplevl])
#uiv = bx.quiver(xx,yy,uu,vv)
ebar = plt.colorbar(featmap,format=ticker.FuncFormatter(fmt))
ncfile = 'MV_FTLE_-6hrs_NoWindage.nc'
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print vars.keys()
ftle = vars['FTLE'][:]
flon = vars['initial_lon'][:]
flat = vars['initial_lat'][:]

ftle = np.reshape(ftle[:,:,0],[ydim,xdim,7])/24
flon = np.reshape(flon,[ydim,xdim])
flat = np.reshape(flat,[ydim,xdim])

fig = plt.figure(2,figsize=[15,12], dpi=150)
ax = plt.subplot(221)
ftlesnap=-3 #-2hr
heatmap = m.contourf(flon,flat,ftle[:,:,ftlesnap],levels=np.linspace(ftle[:,:,ftlesnap].min(),ftle[:,:,ftlesnap].max(),301),latlon=True)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax.set_aspect('equal', adjustable='box', anchor='C')
#qiv = ax.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
cbar = plt.colorbar(heatmap,format=ticker.FuncFormatter(fmt))
plt.title(strftime("-2hr FTLE @ %a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))

bx = plt.subplot(222)
ftlesnap=-5
geatmap = m.contourf(flon,flat,ftle[:,:,ftlesnap],levels=np.linspace(ftle[:,:,ftlesnap].min(),ftle[:,:,ftlesnap].max(),301),latlon=True)#,shading='gourand')
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
bx.set_aspect('equal', adjustable='box', anchor='C')
#uiv = bx.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
#uiv = bx.quiver(xx,yy,uu,vv)
dbar = plt.colorbar(geatmap,format=ticker.FuncFormatter(fmt))
plt.title(strftime("-4hr FTLE @ %a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))

cx = plt.subplot(223)
ftlesnap=-7
featmap = m.contourf(flon,flat,ftle[:,:,ftlesnap],levels=np.linspace(ftle[:,:,ftlesnap].min(),ftle[:,:,ftlesnap].max(),301),latlon=True)#,shading='gourand')
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
cx.set_aspect('equal', adjustable='box', anchor='C')
#wiv = cx.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2])
#uiv = bx.quiver(xx,yy,uu,vv)
ebar = plt.colorbar(featmap,format=ticker.FuncFormatter(fmt))
plt.title(strftime("-6hr FTLE @ %a, %d %b %Y %H:%M:%S Zulu Time", gmtime(time[t])))
'''