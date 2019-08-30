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

f = loadmat('tracers_out_20.mat')
x_t = f['fx']
y_t = f['fy']
x_t=np.delete(x_t,(3),axis=0)
y_t=np.delete(y_t,(3),axis=0)
a_fx=[]
a_fy=[]
r_fx=[]
r_fy=[]
for file in ['a_036','a_048','a_054','a_065','a_083','r_055','r_064','r_085','r_127']:#,'a_137','r_107'
    f = loadmat(file+'_out.mat')
    if file[0]=='a':
        a_fx.append(f['fx'])
        a_fy.append(f['fy'])
    else:
        r_fx.append(f['fx'])
        r_fy.append(f['fy'])
    del f

a_len = len(a_fx)
r_len = len(r_fx)
'''
for t in range(73):
    #for i in range(6):
    m.scatter(x_t[:,:,t],y_t[:,:,t])
    for i in range(a_len):
        m.plot(a_fx[i][:,t],a_fy[i][:,t],'b')
    for i in range(r_len):
        m.plot(r_fx[i][:,t],r_fy[i][:,t],'r')
    m.drawcoastlines()
    plt.savefig('{:04d}.png'.format(t))
    plt.close('all')
#'''

plt.figure(figsize=[5+3/8,5+3/8])    
plt.subplot(221)
t=0
for i in range(a_len):
    m.plot(a_fx[i][:,t],a_fy[i][:,t],'b')
for i in range(r_len):
    m.plot(r_fx[i][:,t],r_fy[i][:,t],'r')
m.scatter(x_t[:,:,t],y_t[:,:,t],color=marker_color,s=marker_size)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('A', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(222)
t=12
for i in range(a_len):
    m.plot(a_fx[i][:,t],a_fy[i][:,t],'b')
for i in range(r_len):
    m.plot(r_fx[i][:,t],r_fy[i][:,t],'r')
m.scatter(x_t[:,:,t],y_t[:,:,t],color=marker_color,s=marker_size)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('B', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(223)
t=24
for i in range(a_len):
    m.plot(a_fx[i][:,t],a_fy[i][:,t],'b')
for i in range(r_len):
    m.plot(r_fx[i][:,t],r_fy[i][:,t],'r')
m.scatter(x_t[:,:,t],y_t[:,:,t],color=marker_color,s=marker_size)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('C', xy=(0.91, 0.02), xycoords='axes fraction')
plt.subplot(224)
t=48
for i in range(a_len):
    m.plot(a_fx[i][:,t],a_fy[i][:,t],'b')
for i in range(r_len):
    m.plot(r_fx[i][:,t],r_fy[i][:,t],'r')
m.scatter(x_t[:,:,t],y_t[:,:,t],color=marker_color,s=marker_size)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.annotate('D', xy=(0.91, 0.02), xycoords='axes fraction')
plt.savefig('ocean_ilcs_tracers.png')
plt.savefig('ocean_ilcs_tracers.eps')
#plt.close('all')
    
    
    
    
    
    
    
    
    
    
    
    
    