# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:48:37 2019

@author: pnola
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def rmse(a,b):
    return np.sqrt(np.sum((a-b)**2))
    #return np.sum((a-b)**2)
    
plt.close('all')
f=np.load('wrf_eulerian_data.npz')
s1 = np.ma.masked_where(f['s1']==0.,f['s1'])
corr1 = np.ma.masked_where(f['corr1']==0.,f['corr1'])
time1 = f['time']
f.close()
t0=time1[0]
plt.figure()
plt.pcolormesh(-s1[-1,:,:])
plt.colorbar()

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
ftle = 1/24*vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
time2 = 24*vars['time'][:]#,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()

plt.figure()
plt.pcolormesh(ftle[-2,:,:])
plt.colorbar()

with open('RMSE_wrf.csv','w') as fd:
        fd.close()
#"""
error_with_s1=[]
error_with_corr=[]
xx=[]
for i in range(1,len(time2)):
    tf = -time2[i]
    T=tf-t0
    xx.append(abs(T))
    ftle_approx = -s1[-1,:,:]-T*corr1[-1,:,:]
    if i==0:
        ftle_true=-s1[-1,:,:]
    else:
        ftle_true = ftle[-i-1,:,:]
    error_with_s1.append(rmse(ftle_true,-s1[-1,:,:]))
    error_with_corr.append(rmse(ftle_true,ftle_approx))
    """
    with open('RMSE_wrf.csv','a') as fd:
        fd.write(str(error_with_s1)+','+str(error_with_corr)+',\n')
        fd.close()
    #""" 
with open('RMSE_wrf.csv','w') as fd:
    fd.write(str(error_with_s1)+'\n')
    fd.write(str(error_with_corr))
    fd.close()


plt.figure(figsize=(8,4))
plt.plot(xx,np.log(error_with_s1),label='-s1')
plt.plot(xx,np.log(error_with_corr),label='-s1-T*correction')
#plt.plot(xx,error_with_s1,label='-s1')
#plt.plot(xx,error_with_corr,label='-s1-T*correction')
plt.legend()
plt.xlabel('Integration Time')
plt.ylabel('log RMSE')
plt.savefig('RMSE_timeseries_wrf.png', transparent=False, bbox_inches='tight',pad_inches=0.03)

'''
    plt.figure(figsize=(8,4))
    plt.subplot(121)
    plt.pcolormesh(ftle_approx)
    plt.colorbar()
    plt.subplot(122)
    plt.pcolormesh(ftle_true)
    plt.colorbar()
    plt.savefig('{0:03d}'.format(i))
    plt.close('all')
    '''
    #data_out.to_csv('Correction_correlation_dg.csv',mode='a',index=False, header=False)#,float_format='%1.3f')
#"""