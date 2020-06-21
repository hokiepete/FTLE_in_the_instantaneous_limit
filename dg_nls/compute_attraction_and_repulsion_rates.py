# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 16:28:35 2019

@author: pnola

This script shows how to compute the attraction and repulsion rate fields (s_1 & s_n) for a 2 dimensional flow.

This script Assumes that the flow is presented in a MESHGRID format, i.e. arrays are structued as [i,j,...]
where i is the index for rows or y-dimension and j is the index for columns or x-dimension.

Note! while this script is written for a 2-dimensional flow, it can be easily extended to n-dimensions.
"""
import numpy as np
import pickle
from velocity_fields import bickley_jet, nonlinear_saddle

#set up your velocity field and domain data
size = 101
x = np.linspace(-0.5,0.5,size)
dx = x[1]-x[0]
y = np.linspace(-0.5,0.5,size)
dy = y[1]-y[0]
x,y=np.meshgrid(x,y)
u,v = nonlinear_saddle(0,[x,y])


#This line isn't necessarary, BUT it allows the script to be easily reused for
#different sized data sets
[ydim,xdim] = u.shape


#Calculate the gradients of the velocity field
dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)

#Initialize arrays for the attraction rate and repullsion rate
#Using masked arrays can be very useful when dealing with geophysical data and
#data with gaps in it.
s1 = np.ma.empty([ydim,xdim])
s2 = np.ma.empty([ydim,xdim])

# For each point in the domain (i,j)
for i in range(ydim):
    for j in range(xdim):
        #Make sure the data is not masked, masked gridpoints do not work with
        #Python's linalg module
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j]) is not np.ma.masked:
            #If the data is not masked, compute s_1 and s_2
            Grad = np.array([[dudx[i,j], dudy[i,j]], [dvdx[i,j], dvdy[i,j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors = np.linalg.eig(S)
            idx = eigenValues.argsort()
            s1[i,j] = eigenValues[idx[0]]
            s2[i,j] = eigenValues[idx[-1]] #The use of -1 here allows this to be more easily extended to n-dimensions.
            
        else:
            #If the data is masked, then mask the grid point in the output.
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked

sigma = s1
with open(f'ftle_0.pickle','wb') as file:
    pickle.dump(sigma,file)

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.pcolormesh(x,y,-s1)

plt.figure(2)
plt.pcolormesh(x,y,s2)

