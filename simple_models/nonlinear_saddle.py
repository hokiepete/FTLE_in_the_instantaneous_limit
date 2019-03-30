# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 22:20:26 2019

@author: pnola
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

dim =7
bd =0.5
x = np.linspace(-bd,bd,dim)
x,y = np.meshgrid(x,x)
xs,ys = np.meshgrid(x,x)
yy = np.concatenate((y[0,:],y[-1,:],[0,0]))
xx = np.concatenate((x[0,:],x[-1,:],[-0.01,0.01]))
sp =np.stack((xx,yy),axis=1)

yy = y[0,:]
xx = x[0,:]
sp1 =np.stack((yy,xx),axis=1)
sp3 =np.stack((xx,yy),axis=1)


yy = y[-1,:]
xx = x[-1,:]
sp2 =np.stack((yy,xx),axis=1)


plt.close('all')
plt.quiver(x,y,x,-y-y**3,pivot='tail')#,scale=10,width=0.01)

dim =15
bdx =1
bdy =1
x = np.linspace(-bdx,bdx,dim+2)
y = np.linspace(-bdy,bdy,dim)
x,y = np.meshgrid(x,y)

#plt.streamplot(x,y,x,-y-y**3,start_points=sp,integration_direction='forward',maxlength=1,minlength=0.25)
#plt.streamplot(x,y,x,-y-y**3,start_points=sp2)
#plt.streamplot(x,y,x,-y-y**3,start_points=sp3)
plt.plot([0,0.2],[0,0],color='C0')
plt.streamplot(x,y,x,-y-y**3,start_points=np.stack((xs.ravel(),ys.ravel()),axis=1),arrowstyle='-',maxlength=1,integration_direction='both')
#plt.quiver(xs,ys,xs,-ys-ys**3,pivot='tail',scale=10,width=0.01)

plt.xlim([-0.6,0.6])
plt.ylim([-0.5,0.5])
#plt.axis('equal')
plt.savefig('nl_saddle.eps', transparent=False, bbox_inches='tight')
