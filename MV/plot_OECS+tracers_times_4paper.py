
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import calendar
import time
from mpl_toolkits.basemap import Basemap
#from scipy.interpolate import spline
from scipy import interpolate
plt.close('all')

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
marker_size = 4
marker_color = 'c'
xdim = 263
ydim = 240
smoothing_coef = 0.1
interp_number = 20
figwidth = 5+3/8
FigSize=(figwidth, ydim/xdim*figwidth)

tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ncfile="mv_ftle.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
#ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )
'''
m = Basemap(llcrnrlon=-179,
            llcrnrlat=-85,
            urcrnrlon=180,
            urcrnrlat=85,
            projection='merc',
            resolution = 'c',
            area_thresh=1000.,
            )
#'''
#lon,lat = np.meshgrid(lon,lat)
#cs=m.contourf(lon,lat,ftle[-1,:,:],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)
#ncfile="SE_tracers.nc"
#ncfile="tracers_750m.nc"
ncfile="tracers_1km.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t_lat = vars["lat"][:]#.reshape([ydim,xdim])
t_lon = vars["lon"][:]#.reshape([ydim,xdim])
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
root.close()
attracting_ridges_lat=[]
attracting_ridges_lon=[]
for a in [13,15,17,26,42,45,49,50]:
    ncfile="a{:d}.nc".format(a)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    attracting_ridges_lat.append(vars["lat"][:])
    attracting_ridges_lon.append(vars["lon"][:])
    root.close()

repelling_ridges_lat=[]
repelling_ridges_lon=[]
for r in [35,37,40,45,48,56,81,99,104,113]:
    ncfile="r{:d}.nc".format(r)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    repelling_ridges_lat.append(vars["lat"][:])
    repelling_ridges_lon.append(vars["lon"][:])
    root.close()
    
#for t in range(time.shape[0]):
#for c in cs.collections:
#    c.remove()
#cs.set_array(np.ravel(ftle[:,:,t]))
#cs=m.contourf(lon,lat,ftle[:,:,-t],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)

#cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True)
#plt.title("{0}".format(time[-t]),fontsize=18)
parallels = np.arange(round(lat_min,0),lat_max+1,0.2)
meridians = np.arange(round(lon_max,0),lon_min-1,-0.2)

plt.figure(figsize=FigSize)
plt.subplot(221)
t=0
m.drawcoastlines()
m.drawstates()
#m.drawparallels(parallels,labels=[1,0,0,0],**tickfont,linewidth=0.0)
#m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont,linewidth=0.0)
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
for i in range(len(repelling_ridges_lat)):
    x = [elem for elem in repelling_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in repelling_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='r',latlon=True)    
    #m.plot(x,y,c='r',latlon=True)
for i in range(len(attracting_ridges_lat)):
    x = [elem for elem in attracting_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in attracting_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='b',latlon=True)
    #m.plot(x,y,c='b',latlon=True)
m.scatter(t_lon[:,t],t_lat[:,t],color=marker_color,s=marker_size,latlon=True)
plt.annotate('A', xy=(0.91, 0.02), xycoords='axes fraction')


plt.subplot(222)
t=5
m.drawcoastlines()
m.drawstates()
#m.drawparallels(parallels,labels=[1,0,0,0],**tickfont,linewidth=0.0)
#m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont,linewidth=0.0)
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
for i in range(len(repelling_ridges_lat)):
    x = [elem for elem in repelling_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in repelling_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='r',latlon=True)    
    #m.plot(x,y,c='r',latlon=True)
for i in range(len(attracting_ridges_lat)):
    x = [elem for elem in attracting_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in attracting_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='b',latlon=True)
    #m.plot(x,y,c='b',latlon=True)
m.scatter(t_lon[:,t],t_lat[:,t],color=marker_color,s=marker_size,latlon=True)
plt.annotate('B', xy=(0.91, 0.02), xycoords='axes fraction')

plt.subplot(223)
t=10
m.drawcoastlines()
m.drawstates()
#m.drawparallels(parallels,labels=[1,0,0,0],**tickfont,linewidth=0.0)
#m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont,linewidth=0.0)
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
for i in range(len(repelling_ridges_lat)):
    x = [elem for elem in repelling_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in repelling_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='r',latlon=True)    
    #m.plot(x,y,c='r',latlon=True)
for i in range(len(attracting_ridges_lat)):
    x = [elem for elem in attracting_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in attracting_ridges_lat[i][:,t].data if elem < 99999]
    tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
    u_new = np.linspace(u.min(), u.max(), interp_number)
    xnew, ynew = interpolate.splev(u_new, tck, der=0)
    m.plot(xnew,ynew,c='b',latlon=True)
    #m.plot(x,y,c='b',latlon=True)
m.scatter(t_lon[:,t],t_lat[:,t],color=marker_color,s=marker_size,latlon=True)
plt.annotate('C', xy=(0.91, 0.02), xycoords='axes fraction')

plt.subplot(224)
t=20
m.drawcoastlines()
m.drawstates()
#m.drawparallels(parallels,labels=[1,0,0,0],**tickfont,linewidth=0.0)
#m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont,linewidth=0.0)
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
for i in range(len(repelling_ridges_lat)):
    x = [elem for elem in repelling_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in repelling_ridges_lat[i][:,t].data if elem < 99999]
    if len(x)>2:
        tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
        u_new = np.linspace(u.min(), u.max(), interp_number)
        xnew, ynew = interpolate.splev(u_new, tck, der=0)
        m.plot(xnew,ynew,c='r',latlon=True)
    else:
        m.plot(x,y,c='r',latlon=True)
for i in range(len(attracting_ridges_lat)):
    x = [elem for elem in attracting_ridges_lon[i][:,t].data if elem < 99999]
    y = [elem for elem in attracting_ridges_lat[i][:,t].data if elem < 99999]
    if len(x)>2:
        tck, u = interpolate.splprep([x, y], u=None, s=smoothing_coef) 
        u_new = np.linspace(u.min(), u.max(), interp_number)
        xnew, ynew = interpolate.splev(u_new, tck, der=0)
        m.plot(xnew,ynew,c='b',latlon=True)
    else:
        m.plot(x,y,c='b',latlon=True)
m.scatter(t_lon[:,t],t_lat[:,t],color=marker_color,s=marker_size,latlon=True)
plt.annotate('D', xy=(0.91, 0.02), xycoords='axes fraction')


plt.savefig('iLCS+tracers_mv.eps'.format(t), transparent=False, bbox_inches='tight')

#'''
    
    
