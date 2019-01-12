# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 21:38:52 2018

@author: pnola
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
A = pd.read_csv('Correlation_and_stats_dg.csv')
t_cor = np.linspace(0,3,101)
y = A['s1']
height=2
plt.figure(figsize=(5/3*height,height))
plt.plot(x,y)
plt.ylabel('Correlation coefficient',**labelfont)
plt.xlabel('Backward-time integration length',**labelfont)

plt.xlabel('Non-dimensional Time')
plt.title('$s_{1}$ vs Backward-time FTLE correlation, Double Gyre')
plt.axis('tight')
dt_cor = t_cor[1]-t_cor[0]
dy_cor = np.gradient(y_cor,dt_cor,edge_order=2)
plt.figure(figsize=(8,6))
plt.plot(t_cor,dy_cor)

"""
xdim = 301
#dimy = int(np.ceil(dimx/2.5))
ydim = int(np.ceil(xdim/2))
tdim = 101
x = np.linspace(0,2,xdim)
dx = x[1]-x[0]
y = np.linspace(0,1,ydim)
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)
time = np.linspace(0,-2,tdim)
dt = time[1]-time[0]
from velocities import double_gyre as vel_func

u = np.empty([tdim,ydim,xdim])
v = np.empty([tdim,ydim,xdim])

for time_step,t in enumerate(time):
    uu,vv = vel_func(t,[x,y])
    u[time_step,:,:]= uu
    v[time_step,:,:]= vv
    
dudt = np.gradient(u,dt,axis=0)
dvdt = np.gradient(v,dt,axis=0)

dudy,dudx = np.gradient(u,dy,dx,axis=(1,2))
dvdy,dvdx = np.gradient(v,dy,dx,axis=(1,2))

Du = dudt+u*dudx+v*dudy
Dv = dvdy+u*dvdx+v*dvdy

dDudy,dDudx = np.gradient(Du,dy,dx,axis=(1,2))
dDvdy,dDvdx = np.gradient(Dv,dy,dx,axis=(1,2))


s1 = np.ma.empty([tdim,ydim,xdim])
s2 = np.ma.empty([tdim,ydim,xdim])
corr1 = np.ma.empty([tdim,ydim,xdim])
corr2 = np.ma.empty([tdim,ydim,xdim])
for t in range(tdim):
    print(t)
    for i in range(ydim):
        for j in range(xdim):
            if (dDudx[t,i,j] and dDudy[t,i,j] and dDvdx[t,i,j] and dDvdy[t,i,j]) is not np.ma.masked:    
                Grad_v = np.array([[dudx[t,i,j], dudy[t,i,j]], [dvdx[t,i,j], dvdy[t,i,j]]])
                Grad_D = np.array([[dDudx[t,i,j], dDudy[t,i,j]], [dDvdx[t,i,j], dDvdy[t,i,j]]])
                S = 0.5*(Grad_v + np.transpose(Grad_v))
                B = (0.5*(Grad_D + np.transpose(Grad_D)+np.matmul(np.transpose(Grad_v),Grad_v)))
                eigenValues, eigenVectors_temp = np.linalg.eig(S)
                idx = eigenValues.argsort()
                eigenMin = eigenVectors_temp[:,idx[0]]
                eigenMax = eigenVectors_temp[:,idx[1]]
                s1[t,i,j] = eigenValues[idx[0]]
                s2[t,i,j] = eigenValues[idx[1]]
                corr1[t,i,j] = -s1[t,i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin)
                corr2[t,i,j] = -s2[t,i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax)
                #sigb[t,i,j] = -s1[t,i,j] - (tf-t0)*(-s1[t,i,j]**2+0.5*np.dot(np.dot(eigenMin,B),eigenMin))
                #sigf[t,i,j] = s2[t,i,j] - (tf-t0)*(-s2[t,i,j]**2+0.5*np.dot(np.dot(eigenMax,B),eigenMax))
            else:
                s1[t,i,j] = np.ma.masked
                s2[t,i,j] = np.ma.masked
                corr1[t,i,j] = np.ma.masked
                corr2[t,i,j] = np.ma.masked
                
np.savez('dg_eulerian_data.npz',s1=s1,s2=s2,corr1=corr1,corr2=corr2)
#"""
with np.load('dg_eulerian_data.npz') as F:
    s1 = F['s1']
    corr1 = F['corr1']
time = np.linspace(0,2,s1.shape[0])
s1_mean = s1.mean(axis=(1,2))
corr1_mean = corr1.mean(axis=(1,2))

plt.figure(figsize=(8,6))
plt.pcolormesh(s1[-1,:,:])
plt.colorbar()
plt.axis('tight')

plt.figure(figsize=(8,6))
plt.plot(time,s1_mean)
plt.ylabel('$hr^{-1}$')
plt.xlabel('Time')
plt.title('$s_{1}$ spatial mean')
plt.axis('tight')
plt.savefig('s1_dg.png')


#dy = np.gradient(y,dt)
plt.figure(figsize=(8,6))
plt.plot(t_cor,dy_cor,label='d/dt correlation')
plt.plot(time,corr1_mean,label='correction term, x48')
plt.legend()
plt.axis('tight')
plt.savefig('CorCoef_vs_correction_notnormalized_dg.png')


dt = time[1]-time[0]
dy = np.gradient(s1_mean,dt)
plt.figure(figsize=(8,6))
plt.plot(time,dy)

dy = np.gradient(corr1_mean,dt)
plt.figure(figsize=(8,6))
plt.plot(time,dy)

ss=s1_mean-corr1_mean
ss = ss - ss.min()
ss = ss/ss.max()
s1_mean = s1_mean - s1_mean.min()
s1_mean = s1_mean/s1_mean.max()
corr1_mean = corr1_mean - corr1_mean.min()
corr1_mean = corr1_mean/corr1_mean.max()
plt.figure(figsize=(8,6))
#plt.plot(time,s1_mean,label='normalize s1, spatial mean')
plt.plot(time,corr1_mean,label='normalize correction term, spatial mean')
#plt.plot(x,y,label='s1 vs FTLE correlation')
dy_cor = dy_cor - dy_cor.min()
dy_cor = dy_cor/dy_cor.max()
plt.plot(t_cor,dy_cor,label='d/dt correlation')
plt.legend()
plt.xlabel('Time')
plt.axis('tight')
plt.savefig('correlation_vs_correction_normalized_dg.png')

"""
data = [s1_mean,corr1_mean,ss,y_cor,dy_cor]
name=['s1_mean','correction_mean','s1_mean-dt*correction_mean','correlation','d/dt correlation']
Alldata = pd.DataFrame(np.transpose(data),columns=name)
Alldata.corr().to_csv('s1+correction_correlation_data_dg.csv',mode='w')

plt.figure()
plt.scatter(y_cor,ss)
plt.xlabel('Correlation Coeff')
plt.ylabel('s1_mean-dt*correction_mean, normalized')
plt.savefig('CorCoef_vs_s1-Tcorrection_dg.png')

plt.figure()
plt.scatter(y_cor,s1_mean)
plt.xlabel('Correlation Coeff')
plt.ylabel('s1_mean, normalized')
plt.savefig('CorCoef_vs_s1_dg.png')

plt.figure()
plt.scatter(y_cor,corr1_mean)
plt.xlabel('Correlation Coeff')
plt.ylabel('correction spatial mean, normalized')
plt.savefig('CorCoef_vs_correction_dg.png')


plt.figure()
plt.scatter(y_cor[:2],corr1_mean[:2])
plt.xlabel('Correlation Coeff')
plt.ylabel('correction spatial mean, normalized')
plt.axis('equal')
plt.savefig('CorCoef_vs_correctionv2_dg.png')

#"""


'''
ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()
ftle_mean = ftle.mean(axis=(1,2))
time = np.linspace(0,24,ftle.shape[0])
plt.figure(figsize=(8,6))
plt.plot(time,ftle_mean[::-1])
plt.ylabel('$hr^{-1}$')
plt.xlabel('Time')
plt.title('FTLE spatial mean')
plt.axis('tight')

#'''