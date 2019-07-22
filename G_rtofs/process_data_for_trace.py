import mapping_functions as mf
from netCDF4 import Dataset
import numpy as np
from hdf5storage import savemat
from scipy.interpolate import griddata

#25N
#-90E
xdim=400
ydim=int(np.ceil(xdim*4/5))
tdim = 25
u_out = np.empty((tdim,ydim,xdim))
v_out = np.empty((tdim,ydim,xdim))
t = np.empty((tdim))

x_out = np.linspace(-900*1000,1100*1000,xdim)
y_out = np.linspace(-800*1000,800*1000,ydim)
xx,yy = np.meshgrid(x_out,y_out)
lon, lat = mf.m2lonlat(-90,25,xx,yy,std_lat1=20,std_lat2=30)
llon = np.linspace(lon[:,1].max(),lon[:,-1].min(),400)
llat = np.linspace(lat[1,:].max(),lat[-1,:].min(),320)
llon,llat = np.meshgrid(llon,llat)
xx,yy = mf.lonlat2m(-90,25,llon,llat,std_lat1=20,std_lat2=30)
Xi= list(zip(xx.ravel(),yy.ravel()))
del xx, yy


for i in range(25):
    print(i)
    root = Dataset('data/07_08_2019/rtofs_glo_2ds_f{0:03d}_1hrly_prog.nc'.format(i),'r')
    vars = root.variables #dictionary, all variables in dataset\
    t[i] = vars['MT'][0]
    lon = vars['Longitude'][:-1,:]-360
    lat = vars['Latitude'][:-1,:]
    u = vars['u_velocity'][:,:,:-1,:].squeeze()
    v = vars['v_velocity'][:,:,:-1,:].squeeze()
    root.close() 
    
    x_in,y_in = mf.lonlat2m(-90,25,lon,lat,std_lat1=20,std_lat2=30)
    del lon, lat
    
    u = u[x_in>=-1000*1000]
    v = v[x_in>=-1000*1000]
    y_in = y_in[x_in>=-1000*1000]
    x_in = x_in[x_in>=-1000*1000]
    
    u = u[x_in<=1200*1000]
    v = v[x_in<=1200*1000]
    y_in = y_in[x_in<=1200*1000]
    x_in = x_in[x_in<=1200*1000]
    
    u = u[y_in>=-900*1000]
    v = v[y_in>=-900*1000]
    x_in = x_in[y_in>=-900*1000]
    y_in = y_in[y_in>=-900*1000]
    
    u = u[y_in<=900*1000]
    v = v[y_in<=900*1000]
    x_in = x_in[y_in<=900*1000]
    y_in = y_in[y_in<=900*1000]
    
    v[u.mask] = np.nan
    u[u.mask] = np.nan
    
    u[v.mask] = np.nan
    v[v.mask] = np.nan
    
    points = list(zip(x_in,y_in))
    
    uu = griddata(points,u,Xi,method='linear')
    vv = griddata(points,v,Xi,method='linear')
    del points

    u_out[i,:,:] = np.reshape(uu,(ydim,xdim))
    v_out[i,:,:] = np.reshape(vv,(ydim,xdim))
    del uu,vv

land = np.zeros((320,400))

llon=llon[1,:]
llat=llat[:,1]

land[np.isnan(u_out[0,:,:])] = 1
u_out[np.isnan(u_out)]=999
v_out[np.isnan(v_out)]=999

dataset = Dataset('gulf.nc', mode='w', format='NETCDF4_CLASSIC') 
lat = dataset.createDimension('lat',llat.shape[0])
lon = dataset.createDimension('lon',llon.shape[0])
time = dataset.createDimension('time', None)

times = dataset.createVariable('time',np.float64,('time',),fill_value=999)
lats = dataset.createVariable('lat',np.float64,('lat'),fill_value=999)
lons = dataset.createVariable('lon',np.float64,('lon',),fill_value=999)
uo = dataset.createVariable('eastward_vel',np.float64,('time','lat','lon',),fill_value=999)
vo = dataset.createVariable('northward_vel',np.float64,('time','lat','lon',),fill_value=999)
lando = dataset.createVariable('land',np.int,('lat','lon',))


lons.standard_name = 'longitude'
lons.units = 'degree_east'
lons.positive = 'east'
lons._CoordinateAxisType = 'Lon'
lons.axis = 'X'
lons.coordsys = 'geographic'
lons[:] = llon

lats.standard_name = 'latitude'
lats.units = 'degree_north'
lats.positive = 'up'
lats._CoordinateAxisType = 'Lat'
lats.axis = 'Y'
lats.coordsys = 'geographic'
lats[:] = llat

times.standard_name = 'time'
times.long_name = 'time'
times.units = 'days since 1901-01-01 00:00:00 UTC'
times.calendar = 'gregorian'
times._CoordinateAxisType = 'Time'
times[:] = t[:]

uo.standard_name = 'surface_eastward_sea_water_velocity'
uo.long_name = 'surface_eastward_sea_water_velocity'
uo.units = 'meter second-1'
uo.coordsys = 'geographic'
uo.positive = 'toward east'
uo.coordinates = 'Longitude Latitude datetime'
uo[:] = u_out

vo.standard_name = 'surface_northward_sea_water_velocity'
vo.long_name = 'surface_northward_sea_water_velocity'
vo.units = 'meter second-1'
vo.coordsys = 'geographic'
vo.positive = 'toward north'
vo.coordinates = 'Longitude Latitude datetime'
vo[:] = v_out


lando.standard_name = 'land'
lando.units = 'Boolean'
lando.coordinates = 'Longitude Latitude'
lando[:] = land

dataset.close()










