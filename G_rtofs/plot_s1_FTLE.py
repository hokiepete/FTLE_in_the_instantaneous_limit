import numpy as np
from hdf5storage import loadmat
import matplotlib
import matplotlib.pyplot as plt
import os
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap
plt.close('all')

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
marker_size = 4
marker_color = 'c'

f = loadmat('rtofs_vel.mat')
lon_0 = f['lon_0']
lat_0 = f['lat_0']
std_lat1 = f['std_lat1']
std_lat2 = f['std_lat2']
lat = f['lat']
lon = f['lon']
x = f['x']
y = f['y']
height=y.max()-y.min()
width=x.max()-x.min()
m = Basemap(
        width=width,
        height=height,
        rsphere=(6378137.00,6356752.3142),
        resolution='i',
        area_thresh=5000.,
        projection='lcc',
        lat_1=std_lat1,
        lat_2=std_lat2,
        lat_0=np.mean(lat),#lat_0,
        lon_0=np.mean(lon)#lon_0
        )
parallels = np.arange(round(lat.min(),0),lat.max()+2,3)
meridians = np.arange(round(lon.max(),0),lon.min()-2,-3)


f = loadmat('error_comparison_terms.mat')
s1 = -np.ma.masked_where(np.isnan(f['s1']),f['s1'])

f = loadmat('FTLE_rtofs.mat')
ftle = np.ma.masked_where(np.isnan(f['sigma']),f['sigma'])
T = f['T']
lev = 301
plt.figure(figsize=[5+3/8,5+3/8])    
plt.subplot(221)
t=0
#m.pcolormesh(lon,lat,-s1[:,:,24],latlon=True)
m.contourf(lon,lat,s1[:,:,24],levels=np.linspace(s1[:,:,24].min(axis=None),s1[:,:,24].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('A', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(222)
t=12
m.contourf(lon,lat,ftle[:,:,4],levels=np.linspace(ftle[:,:,4].min(axis=None),ftle[:,:,4].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('B', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(223)
t=24
m.contourf(lon,lat,ftle[:,:,16],levels=np.linspace(ftle[:,:,16].min(axis=None),ftle[:,:,16].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('C', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(224)
t=48
m.contourf(lon,lat,ftle[:,:,96],levels=np.linspace(ftle[:,:,96].min(axis=None),ftle[:,:,96].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('D', xy=(0.91, 0.02), xycoords='axes fraction')
plt.savefig('ocean_s1_Backward-Time_FTLE_Comparison_v2.png', transparent=False, bbox_inches='tight',dpi=300)
plt.savefig('ocean_s1_Backward-Time_FTLE_Comparison_v2.eps', transparent=False, bbox_inches='tight')
#plt.close('all')
    
    
    
    
    
    
    
    
    
    
    
    
    