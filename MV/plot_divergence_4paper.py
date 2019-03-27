
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import calendar
import time
from mpl_toolkits.basemap import Basemap
plt.close('all')
""
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
tdim = 25
smoothing_coef = 0.1
interp_number = 20

tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
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

dudy ,dudx, qt = np.gradient(u,200,200,3600)
dvdy, dvdx, qt = np.gradient(v,200,200,3600)
del qt
div = dudx+dvdy

s1 = np.ma.empty([ydim,xdim,tdim])

J = np.array([[0, 1], [-1, 0]])
for t in range(tdim):
    for i in range(ydim):
        for j in range(xdim):
            if (dudx[i,j,t] and dudy[i,j,t] and dvdx[i,j,t] and dvdy[i,j,t] and u[i,j,t] and v[i,j,t]) is not np.ma.masked:    
                Grad = np.array([[dudx[i, j,t], dudy[i, j,t]], [dvdx[i, j,t], dvdy[i, j,t]]])
                S = 0.5*(Grad + np.transpose(Grad))
                s1[i,j,t] = np.linalg.eig(S)[0].min()
    
            else:
                s1[i,j,t] = np.ma.masked
                




lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

lon, lat = np.meshgrid(lon,lat)
#"""
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )
from mpl_toolkits.axes_grid1 import make_axes_locatable
figwidth = 5+3/8
FigSize=(figwidth, 0.32*figwidth)
plt.figure(figsize=FigSize)
ax1 = plt.subplot(121)

s1p=-3600*s1[:,:,-1]
cs=m.contourf(lon,lat,s1p,levels=np.linspace(s1p.min(),s1p.max(),301),latlon=True)

m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+1,0.2)
meridians = np.arange(round(lon_max,0),lon_min-1,-0.2)
m.drawparallels(parallels,labels=[0,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,0],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=plt.colorbar(cs, cax=cax, orientation='vertical',format="%.1f");

ax1 = plt.subplot(122)

divp = abs(3600*div[:,:,-1])
cs=m.contourf(lon,lat,divp,levels=np.linspace(divp.min(),divp.max(),301),latlon=True)

m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+1,0.2)
meridians = np.arange(round(lon_max,0),lon_min-1,-0.2)
m.drawparallels(parallels,labels=[0,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,0],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar=plt.colorbar(cs, cax=cax, orientation='vertical',format="%.1f");

plt.savefig('s1_v_div_mv_v2.png', transparent=False, bbox_inches='tight',dpi=300)
#plt.savefig('s1_v_div_mv_v2.eps', transparent=False, bbox_inches='tight')

#'''
    
    
