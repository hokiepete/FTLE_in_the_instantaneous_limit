from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import calendar
plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

with np.load('hossien_wrf_eulerian_data.npz') as F:
    s1 = np.ma.masked_where(F['s1'][0,:,:]==0.,F['s1'][0,:,:])
    s2 = np.ma.masked_where(F['s2'][0,:,:]==0.,F['s2'][0,:,:])

s1=-3600*s1
s2= 3600*s2

root = Dataset('hosiendata_wind_velocity.nc','r')
vars = root.variables
lon = vars['lon'][:]
lat = vars['lat'][:]

figwidth = 5+3/8#6.5
FigSize=(figwidth, 0.32*figwidth)
plt.figure(figsize=FigSize)
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

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
parallels_spacing = 3
meridian_spacing = -4
parallels = np.arange(round(lat_min,0),lat_max+2,parallels_spacing)
meridians = np.arange(round(lon_max,0),lon_min-2,meridian_spacing)
lon,lat = np.meshgrid(lon,lat)
ax1=plt.subplot(121)
cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=0.0)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=0.0)
#plt.title('-s$_{1}$',**titlefont)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=plt.colorbar(cs, cax=cax, orientation='vertical',format="%.2f");
#cbar.ax.tick_params(labelsize=8)

ax2=plt.subplot(122)
cs=m.contourf(lon,lat,s2,levels=np.linspace(s2.min(axis=None),s2.max(axis=None),301),latlon=True)
#cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=0.0)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=0.0)
#plt.title('s$_{2}$',**titlefont)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=plt.colorbar(cs, cax=cax, orientation='vertical',format="%.2f");

plt.savefig('s1_vs_s2_v2.png', transparent=False, bbox_inches='tight',pad_inches=0.03,dpi=300)
plt.savefig('s1_vs_s2_v2.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)
