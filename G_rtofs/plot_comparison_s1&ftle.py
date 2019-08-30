from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

grid_spacing = 200 # m

tdim = 25
xdim = 263
ydim = 240

figwidth = 5+3/8#6
FigSize=(figwidth, ydim/xdim*figwidth)

root = Dataset('windagedata.nc','r')
vars = root.variables
#Wind Velocity
u = vars['eastward_vel'][:]
v = vars['northward_vel'][:]
#Water Vapor Flux, Vertically Integrated
#u = vars['UQ'][:,:,:]
#v = vars['VQ'][:,:,:]
lat = vars['lat'][:]
lon = vars['lon'][:]
root.close()
u=u[-1,:,:]
v=v[-1,:,:]


x = np.linspace(0,grid_spacing*(xdim-1),xdim)
dx = x[1]-x[0]
#y = np.linspace(-grid_spacing*(ydim-1)/2,grid_spacing*(ydim-1)/2,ydim)
y = np.linspace(0,grid_spacing*(ydim-1),ydim)
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)

dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)


s1 = np.ma.empty([ydim,xdim])
s2 = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
            s2[i,j] = 3600*np.max(np.linalg.eig(S)[0])

        else:
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked

plt.close('all')
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)
parallels = np.arange(round(lat_min,0),lat_max+1,0.2)
meridians = np.arange(round(lon_max,0),lon_min-1,-0.2)
lev=301
height=(ydim-1)*dy#+10*1000
width=(xdim-1)*dx#+10*1000
lon, lat = np.meshgrid(lon,lat)
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )
#'''

fig = plt.figure(figsize=FigSize)
sub = plt.subplot(221)
cs = m.contourf(lon,lat,-s1,levels=np.linspace(np.min(-s1,axis=None),np.max(-s1,axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.title("s$_{1}$ field",**titlefont)#fontsize=18)
plt.annotate('A', xy=(0.91, 0.02), xycoords='axes fraction')

ncfile="windage_ftle.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
#lon, lat = np.meshgrid(lon,lat)
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)
lon, lat = np.meshgrid(lon,lat)
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )

sub = plt.subplot(222)
t=7
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
#plt.title("Integration time = -{0:02d}.{1:02g} hrs".format(hrs,mins/60),**titlefont)
plt.annotate('B', xy=(0.91, 0.02), xycoords='axes fraction')

sub = plt.subplot(223)
t=13
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
#plt.title("Integration time = -{0:02d}.{1:02g} hrs".format(hrs,mins/60),**titlefont)
plt.annotate('C', xy=(0.91, 0.02), xycoords='axes fraction')

sub = plt.subplot(224)
t=25
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
plt.annotate('D', xy=(0.91, 0.02), xycoords='axes fraction')

plt.savefig('s1_Backward-Time_FTLE_Comparison_v2.png', transparent=False, bbox_inches='tight',dpi=300)
plt.savefig('s1_Backward-Time_FTLE_Comparison_v2.eps', transparent=False, bbox_inches='tight')
