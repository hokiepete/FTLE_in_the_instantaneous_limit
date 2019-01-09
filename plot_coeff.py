# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 21:38:52 2018

@author: pnola
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

plt.close('all')
A = pd.read_csv('Correlation_and_stats_dg2.csv')
x = np.linspace(0,10,101)
y = A['s1']

plt.figure(figsize=(8,6))
plt.plot(x,y)
plt.ylabel('Pearson correlation coefficient')
plt.xlabel('Non-dimensional Time')
plt.title('$s_{1}$ vs Backward-time FTLE correlation, Double Gyre')
plt.axis('tight')
dt = x[1]-x[0]
dy = np.gradient(y)
plt.figure(figsize=(8,6))
plt.plot(x,dy)


A = pd.read_csv('Correlation_and_stats_v2.2.csv')
x = np.linspace(0,24,145)
y = A['s1']
height=2
plt.figure(figsize=(5/3*height,height))
plt.plot(x,y)
plt.ylabel('Correlation coefficient',**labelfont)
plt.xlabel('Backward-time integration length',**labelfont)
#plt.title('$s_{1}$ vs backward-time FTLE correlation',**titlefont)
plt.axis('tight')
plt.xticks([0,6,12,18,24],**tickfont)
plt.yticks(**tickfont)
plt.savefig('correlation_timeseries.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)
plt.savefig('correlation_timeseries.png', transparent=False, bbox_inches='tight',pad_inches=0.03)
"""
dt = x[1]-x[0]
dy = np.gradient(y,dt)
plt.figure(figsize=(8,6))
plt.plot(x,dy)

from netCDF4 import Dataset
root = Dataset('hosiendata_wind_velocity.nc','r')
vars = root.variables
#Wind Velocity
u = vars['eastward_vel'][:]
v = vars['northward_vel'][:]
root.close()
[tdim,ydim,xdim]=u.shape
grid_spacing = 12000 #m
dt=3600 #1 hrs 2 sec
x = np.linspace(0,grid_spacing*(xdim-1),xdim)
dx = x[1]-x[0]
y = np.linspace(0,grid_spacing*(ydim-1),ydim)
dy = y[1]-y[0]
x, y = np.meshgrid(x,y)


dudy,dudx = np.gradient(u,dy,dx,axis=(1,2))
dvdy,dvdx = np.gradient(v,dy,dx,axis=(1,2))

dudydt = np.gradient(dudy,dt,axis=0)
dudxdt = np.gradient(dudx,dt,axis=0)
dvdydt = np.gradient(dvdy,dt,axis=0)
dvdxdt = np.gradient(dvdx,dt,axis=0)
'''
u=u[time_step,:,:]
v=v[time_step,:,:]
dudy=dudy[time_step,:,:]
dudx=dudx[time_step,:,:]
dvdy=dvdy[time_step,:,:]
dvdx=dvdx[time_step,:,:]
dudydt=dudydt[time_step,:,:]
dudxdt=dudxdt[time_step,:,:]
dvdydt=dvdydt[time_step,:,:]
dvdxdt=dvdxdt[time_step,:,:]
[ydim,xdim]=u.shape
'''



s1 = np.ma.empty([tdim,ydim,xdim])
s2 = np.ma.empty([tdim,ydim,xdim])
corr1 = np.ma.empty([tdim,ydim,xdim])
corr2 = np.ma.empty([tdim,ydim,xdim])
for t in range(tdim):
    print(t)
    for i in range(ydim):
        for j in range(xdim):
            if (dudxdt[t,i,j] and dudydt[t,i,j] and dvdxdt[t,i,j] and dvdydt[t,i,j]) is not np.ma.masked:    
                Grad = np.array([[dudx[t,i,j], dudy[t,i,j]], [dvdx[t,i,j], dvdy[t,i,j]]])
                Grad_dt = np.array([[dudxdt[t,i,j], dudydt[t,i,j]], [dvdxdt[t,i,j], dvdydt[t,i,j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                B = (0.5*(Grad_dt + np.transpose(Grad_dt)+np.dot(np.transpose(Grad),Grad)))
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
                
np.savez('hossien_wrf_eulerian_data.npz',s1=s1,s2=s2,corr1=corr1,corr2=corr2)
#"""
"""
with np.load('hossien_wrf_eulerian_data.npz') as F:
    s1 = 3600*F['s1'][::-1]
    corr1 = 3600**2*F['corr1'][::-1]
time = np.linspace(0,24,s1.shape[0])
s1_mean = s1.mean(axis=(1,2))
corr1_mean = corr1.mean(axis=(1,2))
'''
plt.figure(figsize=(8,6))
plt.pcolormesh(s1[-1,:,:])
plt.colorbar()
plt.axis('tight')

plt.figure(figsize=(8,6))
plt.plot(time,s1_mean)
plt.ylabel('$hr^{-1}$')
plt.xlabel('Hours')
plt.title('$s_{1}$ spatial mean')
plt.axis('tight')
plt.savefig('s1.png')
'''

plt.figure(figsize=(8,6))
plt.plot(x,y,label='correlation')
plt.plot(time,-100*corr1_mean,label='correction term')
plt.legend()


dt = x[1]-x[0]
dy = np.gradient(y,dt)
plt.figure(figsize=(8,6))
plt.plot(x,dy,label='d/dt correlation')
plt.plot(time,19*corr1_mean,label='correction term, x19')
plt.legend()
plt.axis('tight')
plt.savefig('CorCoef_vs_correction_scaled.png')

'''
dt = time[1]-time[0]
dy = np.gradient(s1_mean,dt)
plt.figure(figsize=(8,6))
plt.plot(time,dy)

dy = np.gradient(corr1_mean,dt)
plt.figure(figsize=(8,6))
plt.plot(time,dy)
'''
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
dt = x[1]-x[0]
dy = np.gradient(y,dt)
dy = dy - dy.min()
dy = dy/dy.max()
plt.plot(x,dy,label='d/dt correlation')
plt.legend()
plt.xlabel('hours')
plt.axis('tight')
plt.savefig('correlation_vs_correction_normalized.png')

data = [s1_mean,corr1_mean,ss,y[::6],dy[::6]]
name=['s1_mean','correction_mean','s1_mean-dt*correction_mean','correlation','d/dt correlation']
Alldata = pd.DataFrame(np.transpose(data),columns=name)
Alldata.corr().to_csv('s1+correction_correlation_data.csv',mode='w')

plt.figure()
plt.scatter(y[::6],ss)
plt.xlabel('Correlation Coeff')
plt.ylabel('s1_mean-dt*correction_mean, normalized')
plt.savefig('CorCoef_vs_s1-Tcorrection.png')

plt.figure()
plt.scatter(y[::6],s1_mean)
plt.xlabel('Correlation Coeff')
plt.ylabel('s1_mean, normalized')
plt.savefig('CorCoef_vs_s1.png')

plt.figure()
plt.scatter(y[::6],corr1_mean)
plt.xlabel('Correlation Coeff')
plt.ylabel('correction spatial mean, normalized')
plt.savefig('CorCoef_vs_correction.png')


plt.figure()
plt.scatter(y[::6][:2],corr1_mean[:2])
plt.xlabel('Correlation Coeff')
plt.ylabel('correction spatial mean, normalized')
plt.axis('equal')
plt.savefig('CorCoef_vs_correctionv2.png')
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
plt.xlabel('Hours')
plt.title('FTLE spatial mean')
plt.axis('tight')

#'''