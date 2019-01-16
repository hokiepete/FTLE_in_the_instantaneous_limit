# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:48:37 2019

@author: pnola
"""

import numpy as np
import matplotlib.pyplot as plt

def rmse(a,b):
    return np.sqrt(np.sum((a-b)**2))
    #return np.sum((a-b)**2)
    
plt.close('all')
f=np.load('dg_eulerian_data.npz')
s1 = f['s1']
corr1 = f['corr1']
time1 = f['time']
f.close()
t0=time1[0]


plt.figure(figsize=(8,4))
plt.subplot(121)
plt.pcolormesh(-s1[0,:,:])
plt.colorbar()

plt.subplot(122)
plt.pcolormesh(corr1[0,:,:])
plt.colorbar()
plt.title('Numerical')
#f=np.load('dg_ftle_data_short.npz')
f=np.load('dg_ftle_data_long.npz')
ftle = f['ftle']
time2 = f['time']
f.close()

error_with_s1=[]
error_with_corr=[]
xx=[]
for i in range(1,len(time2)):
    tf = -time2[i]
    T=tf-t0
    xx.append(abs(T))
    ftle_approx = -s1[0,:,:]-T*corr1[0,:,:]
    if i==0:
        ftle_true=-s1[0,:,:]
    else:
        ftle_true = ftle[i,:,:]
    error_with_s1.append(rmse(ftle_true,-s1[0,:,:]))
    error_with_corr.append(rmse(ftle_true,ftle_approx))
    """
    plt.figure(figsize=(8,4))
    plt.subplot(121)
    plt.pcolormesh(ftle_approx)
    plt.colorbar()
    plt.subplot(122)
    plt.pcolormesh(ftle_true)
    plt.colorbar()
    plt.savefig('{0:03d}'.format(i))
    plt.close('all')
    """
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
plt.savefig('RMSE_long_timeseries_dg.png', transparent=False, bbox_inches='tight',pad_inches=0.03)

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