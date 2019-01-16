# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:48:37 2019

@author: pnola
"""

import numpy as np
import pandas as pd
from netCDF4 import Dataset
f=np.load('wrf_eulerian_data.npz')
s1 = np.ma.masked_where(f['s1']==0.,f['s1'])
corr1 = np.ma.masked_where(f['corr1']==0.,f['corr1'])
time1 = f['time']
f.close()
t0=time1[0]

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
ftle = 1/24*vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
time2 = 24*vars['time'][:]#,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()

with open('Correction_correlation_wrf.csv','w') as fd:
        fd.write('1.0,\n')
        fd.close()

for i in range(1,len(time2)):
    tf = -time2[i]
    T=tf-t0
    ftle_approx = -s1[-1,:,:]-T*corr1[-1,:,:]
    ftle_true = ftle[-i-1,:,:]
    data = pd.DataFrame(np.transpose([ftle_approx.ravel(),ftle_true.ravel()]),columns=['approx','true'])
    correlation = data.corr()
    data_out = correlation['approx'][1]
    with open('Correction_correlation_wrf.csv','a') as fd:
        fd.write(str(data_out)+',\n')
        fd.close()
    #data_out.to_csv('Correction_correlation_dg.csv',mode='a',index=False, header=False)#,float_format='%1.3f')
