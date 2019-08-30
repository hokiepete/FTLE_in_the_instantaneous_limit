import numpy as np
from scipy import interpolate
import mapping_functions as mf
import csv    
from hdf5storage import loadmat


f = loadmat('rtfos_vel.mat')
lon_0 = f['lon_0']
lat_0 = f['lat_0']
std_lat1 = f['std_lat1']
std_lat2 = f['std_lat2']
lat = f['lat']
lon = f['lon']
x = f['x']
y = f['y']
u = f['u']
v = f['v']
t = f['t']
t = (t-t.min())*24*3600

ydim, xdim, tdim = u.shape

u=u[:,:,0]
v=v[:,:,0]
dx = x[1]-x[0]
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)

dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)


s1 = np.ma.empty([ydim,xdim])
eigenMin = np.ma.empty([ydim,xdim,2])
s2 = np.ma.empty([ydim,xdim])
eigenMax = np.ma.empty([ydim,xdim,2])
#det = np.ma.empty([ydim,xdim])

J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if sum(np.isnan([dudx[i,j], dudy[i,j], dvdx[i,j], dvdy[i,j], u[i,j], v[i,j]]))==0:
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors_temp = np.linalg.eig(S)
            idx = eigenValues.argsort()
            eigenMin[i,j,:] = eigenVectors_temp[:,idx[0]]
            #eigenMin[i,j,0] = eigenVectors_temp[1,idx[0]]
            s1[i,j] = eigenValues[idx[0]]
            eigenMax[i,j,:] = eigenVectors_temp[:,idx[-1]]
            #eigenMax[i,j,0] = eigenVectors_temp[1,idx[-1]]
            s2[i,j] = eigenValues[idx[-1]]

        else:
            s1[i,j] = np.nan
            eigenMin[i,j,0] = np.nan
            eigenMin[i,j,1] = np.nan
            s2[i,j] = np.nan
            eigenMax[i,j,0] = np.nan
            eigenMax[i,j,1] = np.nan
            


s1 = np.ma.masked_where(np.isnan(s1),s1)
eigenMin = np.ma.masked_where(np.isnan(eigenMin),eigenMin)
s2 = np.ma.masked_where(np.isnan(s2),s2)
eigenMax = np.ma.masked_where(np.isnan(eigenMax),eigenMax)


adfdy,adfdx = np.gradient(s1,dy,dx,edge_order=2)
adfdydy,adfdydx = np.gradient(adfdy,dy,dx,edge_order=2)
adfdxdy,adfdxdx = np.gradient(adfdx,dy,dx,edge_order=2)
adirdiv = np.ma.empty([ydim,xdim])
aconcav = np.ma.empty([ydim,xdim])


rdfdy,rdfdx = np.gradient(s2,dy,dx,edge_order=2)
rdfdydy,rdfdydx = np.gradient(rdfdy,dy,dx,edge_order=2)
rdfdxdy,rdfdxdx = np.gradient(rdfdx,dy,dx,edge_order=2)
rdirdiv = np.ma.empty([ydim,xdim])
rconcav = np.ma.empty([ydim,xdim])
"""
for i in range(ydim):
    for j in range(xdim):
        if sum(np.isnan([adfdx[i,j], adfdy[i,j], adfdxdy[i,j], adfdydy[i,j], adfdxdx[i,j], adfdydx[i,j]]))==0:    
            adirdiv[i,j] = np.dot([adfdx[i,j],adfdy[i,j]],eigenMin[i,j,:])
            aconcav[i,j] = np.dot(np.dot([[adfdxdx[i,j],adfdxdy[i,j]],[adfdydx[i,j],adfdydy[i,j]]],eigenMin[i,j,:]),eigenMin[i,j,:])
        #print aconcav[i,j]
        else:
            adirdiv[i,j] = np.nan
            aconcav[i,j] = np.nan

        if (rdfdx[i,j] and rdfdy[i,j] and rdfdxdy[i,j] and rdfdydy[i,j] and rdfdxdx[i,j] and rdfdydx[i,j]) is not np.nan:    
            rdirdiv[i,j] = np.dot([rdfdx[i,j],rdfdy[i,j]],eigenMax[i,j,:])
            rconcav[i,j] = np.dot(np.dot([[rdfdxdx[i,j],rdfdxdy[i,j]],[rdfdydx[i,j],rdfdydy[i,j]]],eigenMax[i,j,:]),eigenMax[i,j,:])
        #print aconcav[i,j]
        else:
            rdirdiv[i,j] = np.nan
            rconcav[i,j] = np.nan
"""

for i in range(ydim):
    for j in range(xdim):
        if (adfdx[i,j] and adfdy[i,j] and adfdxdy[i,j] and adfdydy[i,j] and adfdxdx[i,j] and adfdydx[i,j]) is not np.ma.masked:    
            adirdiv[i,j] = np.dot([adfdx[i,j],adfdy[i,j]],eigenMin[i,j,:])
            aconcav[i,j] = np.dot(np.dot([[adfdxdx[i,j],adfdxdy[i,j]],[adfdydx[i,j],adfdydy[i,j]]],eigenMin[i,j,:]),eigenMin[i,j,:])
        #print aconcav[i,j]
        else:
            adirdiv[i,j] = np.ma.masked
            aconcav[i,j] = np.ma.masked

        if (rdfdx[i,j] and rdfdy[i,j] and rdfdxdy[i,j] and rdfdydy[i,j] and rdfdxdx[i,j] and rdfdydx[i,j]) is not np.ma.masked:    
            rdirdiv[i,j] = np.dot([rdfdx[i,j],rdfdy[i,j]],eigenMax[i,j,:])
            rconcav[i,j] = np.dot(np.dot([[rdfdxdx[i,j],rdfdxdy[i,j]],[rdfdydx[i,j],rdfdydy[i,j]]],eigenMax[i,j,:]),eigenMax[i,j,:])
        #print aconcav[i,j]
        else:
            rdirdiv[i,j] = np.ma.masked
            rconcav[i,j] = np.ma.masked

adirdiv = np.ma.masked_where(np.isnan(adirdiv),adirdiv)
aconcav = np.ma.masked_where(np.isnan(aconcav),aconcav)
adirdiv = np.ma.masked_where(aconcav<0,adirdiv)



rdirdiv = np.ma.masked_where(np.isnan(rdirdiv),rdirdiv)
rconcav = np.ma.masked_where(np.isnan(rconcav),rconcav)
rdirdiv = np.ma.masked_where(rconcav>0,rdirdiv)

import matplotlib.pyplot as plt
import os
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
#plt.register_cmap(name='co', data=mf.co())
plt.close('all')


height=y.max()-y.min()
width=x.max()-x.min()

m = Basemap(
        width=width,
        height=height,
        rsphere=(6378137.00,6356752.3142),
        resolution='l',
        area_thresh=5000.,
        projection='lcc',
        lat_1=std_lat1,
        lat_2=std_lat2,
        lat_0=lat_0,
        lon_0=lon_0
        )
parallels = np.arange(round(lat.min(),0),lat.max()+2,2)
meridians = np.arange(round(lon.max(),0),lon.min()-2,-2)

"""
aridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='blue')
rridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
#plt.close('all')
pp = aridge.collections[0].get_paths()
for p in [13,17,40,46,47,55]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    x,y = m(x,y,inverse=True)
    with open('attracting_ridge_{:02d}.txt'.format(p), 'w', newline='') as f:
        writer = csv.writer(f)
        for row in zip(x,y):
            writer.writerow(row)
        f.close()
    
pp = rridge.collections[0].get_paths()
for p in [8,18,19,53,56,60]:
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    x,y = m(x,y,inverse=True)        
    with open('repelling_ridge_{:02d}.txt'.format(p), 'w', newline='') as f:
        writer = csv.writer(f)
        for row in zip(x,y):
            writer.writerow(row)
        f.close()
    

m.drawcoastlines()
#plt.savefig('__First.png')

ax = plt.gca()
def format_coord(x, y):
    return 'x={0[0]:.4f}, y={0[1]:.4f}'.format(m(x, y, inverse = True))
ax.format_coord = format_coord
r=15#km
theta = np.arange(0,2*np.pi,0.1)
xx = r*np.cos(theta)
yy = r*np.sin(theta)
#x = x+967075/500
#y = y-694557/1000
latlon=[(-80.3136,33.5022),(-79.6674,34.3957),(-80.6516,34.4237),(-81.6697,36.3055),(-84.1812,34.3397),(-81.3643,31.3464),(-84.2830,31.4044),(-87.9021,30.7102),(-82.7258,34.0120)]
x=[]
y=[]
for i, latlon in enumerate(latlon):
    print(i)
    print(latlon[0])
    print(latlon[1])
    w,z = mf.km2lonlat(latlon[0],latlon[1],xx,yy,true_lat1,true_lat2)
    x.append(w)
    y.append(z)
x = [val for sublist in x for val in sublist]
y = [val for sublist in y for val in sublist]
m.scatter(x,y,latlon=True)

with open('tracers.txt', 'w', newline='') as f:
    writer = csv.writer(f)
    for row in zip(x,y):
        writer.writerow(row)
    f.close()

#"""
'''
fig = plt.figure(1)
plt.subplot(121)
#lon,lat = np.meshgrid(lon,lat)
cs = m.contourf(lon,lat,s1,levels=np.linspace(np.min(s1,axis=None),np.max(s1,axis=None),301),latlon=True)
#plt.colorbar()
ridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='cyan')
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

plt.subplot(122)
#lon,lat = np.meshgrid(lon,lat)
cs = m.contourf(lon,lat,s2,levels=np.linspace(np.min(s2,axis=None),np.max(s2,axis=None),301),latlon=True)
#plt.colorbar()
ridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='magenta')
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

fig = plt.figure(2)
#m.quiver(lon[::5,::5],lat[::5,::5],u[::5,::5],v[::5,::5],latlon=True)
ridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='blue')
ridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
m.drawcountries()
#ridge = m.contour(lon,lat,det,levels=0,latlon=True)
ax = plt.gca()
def format_coord(x, y):
    return 'x={0[0]:.4f}, y={0[1]:.4f}'.format(m(x, y, inverse = True))
ax.format_coord = format_coord
r=15#km
theta = np.arange(0,2*np.pi,0.1)
x = r*np.cos(theta)
y = r*np.sin(theta)
#x = x+967075/500
#y = y-694557/1000
x,y = mf.km2lonlat(-90,25,x,y,20,30)
m.scatter(x,y,latlon=True)
import csv
v = np.stack((x,y),axis=1)
with open("tracers.txt", "w") as f:
    writer = csv.writer(f)
    for row in v:
        writer.writerow(row)
    f.close()

m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#'''
'''
v = ridge.collections[0].get_paths()[40].vertices
x,y = m(v[:,0],v[:,1],inverse=True)
v = np.stack((x,y),axis=1)
with open("ridge.txt", "w", newline='') as f:
    writer = csv.writer(f)
    for row in v:
        writer.writerow(row)
    f.close()


#'''



adirdiv = np.ma.masked_where(-s1<=np.percentile(-s1[~np.isnan(s1)],65),adirdiv)
rdirdiv = np.ma.masked_where(s2<=np.percentile(s2[~np.isnan(s2)],65),rdirdiv)

import os
import glob

files = glob.glob('*.png')
for f in files:
    os.remove(f)

aridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='blue')
rridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
m.drawcoastlines()
plt.savefig('{:04d}_all.png'.format(000))
plt.close('all')

#size_a = []
theta = np.linspace(0,2*np.pi,101)
#7-8-2019
#px = [225830,350674,494974,626303,771414,917335]
#py = [503251,237349,228432,327334,147364,422994]
#7-25-2019
px = [48292,152058,486867,718720,727638,442280,489299]
py = [355708,240592,607017,581076,412455,294097,260048]

r = 25000

x = [r*np.cos(theta)+px_i for px_i in px]
y = [r*np.sin(theta)+py_i for py_i in py]
#x = r*np.cos(theta)+px
#y = r*np.sin(theta)+py
m.scatter(x,y)
from hdf5storage import savemat
x_o = np.asarray(x)#np.array([x[0],x[1],x[2],x[3],x[4],x[5],x[6]])
y_o = np.asarray(y)#np.array([y[0],y[1],y[2],y[3],y[4],y[5],y[6]])
savemat('tracers.mat',{'x':x_o,'y':y_o})
a=[36,48,54,65,83,137]
r=[55,64,85,107,127]
pp = aridge.collections[1].get_paths()
qq = rridge.collections[1].get_paths()
for p2 in a:#,255range(len(pp)):
    v = pp[p2].vertices
    x = v[:,0]
    y = v[:,1]
    savemat('a_{:03d}.mat'.format(p2),{'x':x,'y':y})
    #size_a.append(x.size)
    m.plot(x,y,'b-')#, latlon=True)
    
        #size_r=[]
for q2 in r:#,309]: #:range(len(qq)):
    v = qq[q2].vertices
    x = v[:,0]
    y = v[:,1]
    savemat('r_{:03d}.mat'.format(q2),{'x':x,'y':y})
    #size_r.append(x.size)
    m.plot(x,y,'r-')#, latlon=True)
m.drawcoastlines()
plt.title('{:04d}'.format(p2))
plt.savefig('attracting{:03d}.png'.format(p2))
#plt.close('all')
#"""


"""
a_thresh = 12
r_thresh = 12
print('attracting')
pp = aridge.collections[1].get_paths()
qq = rridge.collections[1].get_paths()
for p in range(len(pp)):
    v = pp[p].vertices
    xx = v[:,0]
    yy = v[:,1]
    if xx.size >= a_thresh:
        for p2 in range(len(pp)):
            v = pp[p2].vertices
            x = v[:,0]
            y = v[:,1]
            #size_a.append(x.size)
            if x.size >= a_thresh:
                m.plot(x,y,'b-')#, latlon=True)
                #m.drawcoastlines()
                #plt.title('{:04d}'.format(p))
                #plt.savefig('{:03d}_attracting.png'.format(p))
                #plt.close('all')
    
        #size_r=[]
        for q2 in range(len(qq)): #[58,75,96,118,185,220,309,315]:
            v = qq[q2].vertices
            x = v[:,0]
            y = v[:,1]
            #size_r.append(x.size)
            if x.size >= r_thresh:
                m.plot(x,y,'r-')#, latlon=True)
                #m.drawcoastlines()
                #plt.title('{:04d}'.format(p))
                #plt.savefig('{:03d}_repelling.png'.format(p))
                #plt.close('all')
        m.plot(xx,yy,'g-')
        plt.title('{:04d}'.format(p))
        plt.savefig('attracting{:03d}.png'.format(p))
        plt.close('all')

print('repelling')        
for q in range(len(qq)):
    v = qq[q].vertices
    xx = v[:,0]
    yy = v[:,1]
    if xx.size >= r_thresh:
        for p2 in range(len(pp)):
            v = pp[p2].vertices
            x = v[:,0]
            y = v[:,1]
            #size_a.append(x.size)
            if x.size >= a_thresh:
                m.plot(x,y,'b-')#, latlon=True)
                #m.drawcoastlines()
                #plt.title('{:04d}'.format(p))
                #plt.savefig('{:03d}_attracting.png'.format(p))
                #plt.close('all')
    
        #size_r=[]
        for q2 in range(len(qq)):
            v = qq[q2].vertices
            x = v[:,0]
            y = v[:,1]
            #size_r.append(x.size)
            if x.size >= r_thresh:
                m.plot(x,y,'r-')#, latlon=True)
                #m.drawcoastlines()
                #plt.title('{:04d}'.format(p))
                #plt.savefig('{:03d}_repelling.png'.format(p))
                #plt.close('all')
        m.plot(xx,yy,'g-')
        plt.title('{:04d}'.format(q))
        plt.savefig('repelling{:03d}.png'.format(q))
        plt.close('all')

#"""
"""
a_thresh = 12
r_thresh = 12
print('attracting')
pp = aridge.collections[1].get_paths()
qq = rridge.collections[1].get_paths()
for p in [36,48,54,65,83,137]:#range(len(pp)):
    v = pp[p].vertices
    xx = v[:,0]
    yy = v[:,1]
    m.plot(xx,yy,'r-')

print('repelling')        
for q in [55,64,85,107,127]:#range(len(qq)):
    v = qq[q].vertices
    xx = v[:,0]
    yy = v[:,1]
    m.plot(xx,yy,'b-')

plt.savefig('repelling{:03d}.png'.format(q))

"""









#plt.close('all')
#'''    
     
'''
r=15#km
theta = np.arange(0,2*np.pi,0.1)
x = r*np.cos(theta)
y = r*np.sin(theta)
x,y = mf.km2lonlat(-80.3136,33.5022,x,y,true_lat1,true_lat2)
m.scatter(x,y,latlon=True)
#'''




'''
t=1
hrs, mins = np.divmod((t-1)*10,60)
plt.title("Integration time = {0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
plt.savefig('SE_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
'''
