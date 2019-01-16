# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:48:37 2019

@author: pnola
"""

import numpy as np
import pandas as pd
f=np.load('dg_eulerian_data.npz')
s1 = f['s1']
corr1 = f['corr1']
time1 = f['time']
f.close()
t0=time1[0]

f=np.load('dg_ftle_data_long.npz')
#f=np.load('dg_ftle_data_short.npz')
ftle = f['ftle']
time2 = f['time']
f.close()

with open('Correction_correlation_dg.csv','w') as fd:
        fd.write('1.0,\n')
        fd.close()

for i in range(1,len(time1)):
    tf = time2[i]
    T=tf-t0
    ftle_approx = -s1[0,:,:]-T*corr1[0,:,:]
    ftle_true = ftle[i,:,:]
    data = pd.DataFrame(np.transpose([ftle_approx.ravel(),ftle_true.ravel()]),columns=['approx','true'])
    correlation = data.corr()
    data_out = correlation['approx'][1]
    with open('Correction_correlation_dg.csv','a') as fd:
        fd.write(str(data_out)+',\n')
        fd.close()
    #data_out.to_csv('Correction_correlation_dg.csv',mode='a',index=False, header=False)#,float_format='%1.3f')
