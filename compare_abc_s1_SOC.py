# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 11:28:01 2018

@author: pnola
def kevrekidis(t,Y, eps=0.01):
  return np.array([-Y[0]-Y[1]+2, 1/eps*(Y[0]**3-Y[1])])

def ex11(t,Y):
  return -np.array([np.tanh(Y[0]**2/4)+Y[0],Y[0]+2*Y[1]])

def rotHoop(t,Y, eps=0.1, gamma=2.3):
  return np.array([Y[1], 1/eps*(np.sin(Y[0])*(gamma*np.cos(Y[0])-1)-Y[1])])


def verhulst(t,Y, eps=0.01):
  return [1, 1/eps*(Y[0]*Y[1]-Y[1]**2)]


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from velocities import double_gyre as vel_func
plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
dimx = 301
#dimy = int(np.ceil(dimx/2.5))
dimy = int(np.ceil(dimx/2))
t0 = 0
tf =0.1#days
dimt =2
t = np.linspace(t0,tf,dimt)
dt = t[1]-t[0]
x = np.linspace(0,2,dimx)
y = np.linspace(0,1,dimy)
dx = x[1]-x[0]
dy = y[1]-y[0]

#[x,y,t]=np.meshgrid(x,y,t)
[x,y]=np.meshgrid(x,y)

dudy = np.empty([dimy,dimx,dimt])
dudx = np.empty([dimy,dimx,dimt])
dvdy = np.empty([dimy,dimx,dimt])
dvdx = np.empty([dimy,dimx,dimt])
for tt in range(dimt):
    u,v = vel_func(t[tt],[x,y])
    dudy[:,:,tt],dudx[:,:,tt] = np.gradient(u,dy,dx)#,axis=(0,1))
    dvdy[:,:,tt],dvdx[:,:,tt] = np.gradient(v,dy,dx)#,axis=(0,1))


dudydt = np.gradient(dudy,dt,axis=2)
dudxdt = np.gradient(dudx,dt,axis=2)
dvdydt = np.gradient(dvdy,dt,axis=2)
dvdxdt = np.gradient(dvdx,dt,axis=2)

time_step=0
dudy=dudy[:,:,time_step]
dudx=dudx[:,:,time_step]
dvdy=dvdy[:,:,time_step]
dvdx=dvdx[:,:,time_step]
dudydt=dudydt[:,:,time_step]
dudxdt=dudxdt[:,:,time_step]
dvdydt=dvdydt[:,:,time_step]
dvdxdt=dvdxdt[:,:,time_step]

s1 = np.ma.empty([dimy,dimx])
s2 = np.ma.empty([dimy,dimx])
corr1 = np.ma.empty([dimy,dimx])
corr2 = np.ma.empty([dimy,dimx])
for i in range(dimy):
    for j in range(dimx):
        if (dudxdt[i,j] and dudydt[i,j] and dvdxdt[i,j] and dvdydt[i,j]) is not np.ma.masked:    
            Grad = np.array([[dudx[i,j], dudy[i,j]], [dvdx[i,j], dvdy[i,j]]])
            Grad_dt = np.array([[dudxdt[i,j], dudydt[i,j]], [dvdxdt[i,j], dvdydt[i,j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            B = (0.5*(Grad_dt + np.transpose(Grad_dt)+np.dot(np.transpose(Grad),Grad)))
            eigenValues, eigenVectors_temp = np.linalg.eig(S)
            idx = eigenValues.argsort()
            eigenMin = eigenVectors_temp[:,idx[0]]
            eigenMax = eigenVectors_temp[:,idx[1]]
            s1[i,j] = eigenValues[idx[0]]
            s2[i,j] = eigenValues[idx[1]]
            corr1[i,j] = -s1[i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin)
            corr2[i,j] = -s2[i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax)
            #sigb[i,j] = -s1[i,j] - (tf-t0)*(-s1[t,i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin))
            #sigf[i,j] = s2[i,j] - (tf-t0)*(-s2[t,i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax))
        else:
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked
            corr1[i,j] = np.ma.masked
            corr2[i,j] = np.ma.masked


fig=plt.figure(figsize=(4,4))
#fig, ax = plt.subplots()
plt.subplot(211)
#cs=m.pcolormesh(lon,lat,s1,latlon=True)
#cs=plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
cs=plt.pcolormesh(x,y,s1)
#plt.colorbar()
plt.colorbar( fraction=0.04, pad=0.04)
plt.title('$s_{1}$')


plt.subplot(212)
#cs=m.pcolormesh(lon,lat,corr1,latlon=True)
#cs=plt.contourf(x,y,corr1,levels=np.linspace(corr1.min(axis=None),corr1.max(axis=None),301))
cs=plt.pcolormesh(x,y,corr1)
plt.colorbar( fraction=0.04, pad=0.04)
plt.title('$-s_{1}^{2}+0.5 \\langle \\langle \\xi_{s_{1}}, B \\rangle, \\xi_{s_{1}}\\rangle$')

#plt.axis('tight')