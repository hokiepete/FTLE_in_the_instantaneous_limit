import hdf5storage as h
from numpy import linspace
import matplotlib.pyplot as p
import matplotlib

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
p.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

f=h.loadmat('dg_eulerian_data.mat')
x=f['x']
y=f['y']
s1=f['s1']
c1=f['c1']
t0 = 0
tf = -0.96
time = linspace(t0,tf,101)
f=h.loadmat('dg_sigma_big.mat')
sig=f['sigma']

i = 31
p.close('all')
width =5+3/8
p.figure(figsize=(width,width*1.5))
p.subplot(311)
ft=sig[i,:,:]
p.contourf(x,y,ft,levels=linspace(ft.min(),ft.max(),301))
p.xticks([])
#p.xlabel('x',**labelfont)
p.ylabel('y',**labelfont)
cbar=p.colorbar()        
cbar.set_ticks(linspace(0,1,6))

p.subplot(312)
approx=-s1-time[i]*c1
p.contourf(x,y,approx,levels=linspace(approx.min(),approx.max(),301))
p.ylabel('y',**labelfont)
p.xlabel('x',**labelfont)
cbar=p.colorbar()        
cbar.set_ticks(linspace(0,1,6))

p.subplot(313)
approx=-s1
p.contourf(x,y,approx,levels=linspace(approx.min(),approx.max(),301))
p.xticks([])
p.ylabel('y',**labelfont)
cbar=p.colorbar()        
cbar.set_ticks(linspace(0,1,6))

p.savefig('dg_ftle_vs_s1_v2.png', transparent=False, bbox_inches='tight',pad_inches=0.03)


"""
i = 31
p.close('all')
width =5+3/8
p.figure(figsize=(width,width/3))
p.subplot(121)
ft=sig[i,:,:]
p.contourf(x,y,ft,levels=linspace(ft.min(),ft.max(),301))
p.xlabel('x',**labelfont)
p.ylabel('y',**labelfont)
cbar=p.colorbar()        
cbar.set_ticks(linspace(0,1,6))
#cbar.set_yticklabels(linspace(0,1,11))
p.subplot(122)
approx=-s1-time[i]*c1
p.contourf(x,y,approx,levels=linspace(approx.min(),approx.max(),301))
p.yticks([])
p.xlabel('x',**labelfont)
cbar=p.colorbar()        
cbar.set_ticks(linspace(0,1,6))
p.savefig('dg_ftle_vs_s1.png', transparent=False, bbox_inches='tight',pad_inches=0.03)
"""