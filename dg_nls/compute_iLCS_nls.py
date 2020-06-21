# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 16:28:35 2019

@author: pnola

This script shows how to compute the attracting and repelling iLCS for a 2 dimensional flow.

This script Assumes that the flow is presented in a MESHGRID format, i.e. arrays are structued as [i,j,...]
where i is the index for rows or y-dimension and j is the index for columns or x-dimension.

Note! while this script is written for a 2-dimensional flow, it can be easily extended to n-dimensions.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from velocity_fields import bickley_jet, nonlinear_saddle, double_gyre
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
plt.close('all')
#set up your velocity field and domain data
size = 100
x = 0.5*np.linspace(-1,1,size)
dx = x[1]-x[0]
y = 0.5*np.linspace(-1,1,size)
dy = y[1]-y[0]
x,y=np.meshgrid(x,y)
u,v = nonlinear_saddle(0,[x,y])


#This line isn't necessarary, BUT it allows the script to be easily reused for
#different sized data sets
[ydim,xdim] = u.shape

#FIRST we need to calculate the attraction and repulsion rates, along with their
#associated eigenvectors

#Calculate the gradients of the velocity field
dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)

#Initialize arrays for the attraction rate, repullsion rate, and eigenvectors
#Using masked arrays can be very useful when dealing with geophysical data and
#data with gaps in it.
s1 = np.ma.empty([ydim,xdim])
Xi1 = np.ma.empty([ydim,xdim,2])
s2 = np.ma.empty([ydim,xdim])
Xi2 = np.ma.empty([ydim,xdim,2])

# For each point in the domain (i,j)
for i in range(ydim):
    for j in range(xdim):
        #Make sure the data is not masked, masked gridpoints do not work with
        #Python's linalg module
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j]) is not np.ma.masked:
            #If the data is not masked, compute s_1, s_2 and eigenvectors
            Grad = np.array([[dudx[i,j], dudy[i,j]], [dvdx[i,j], dvdy[i,j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors = np.linalg.eig(S)
            idx = eigenValues.argsort()
            s1[i,j] = eigenValues[idx[0]]
            s2[i,j] = eigenValues[idx[-1]] #The use of -1 here allows this to be more easily extended to n-dimensions.
            Xi1[i,j,:] = eigenVectors[:,idx[0]]
            Xi2[i,j,:] = eigenVectors[:,idx[-1]]
            
        else:
            #If the data is masked, then mask the grid point in the output.
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked
            Xi1[i,j,0] = np.ma.masked
            Xi1[i,j,1] = np.ma.masked
            Xi2[i,j,0] = np.ma.masked
            Xi2[i,j,1] = np.ma.masked


#Now we need to calculate the directional derivative of the attraction and repulsion rates
#W.R.T. the direction of stretching, along with the concavity of the rates.

#Calculate gradients of the s1 field
ds1dy,ds1dx = np.gradient(s1,dy,dx,edge_order=2)
ds1dydy,ds1dydx = np.gradient(ds1dy,dy,dx,edge_order=2)
ds1dxdy,ds1dxdx = np.gradient(ds1dx,dy,dx,edge_order=2)

#Calculate gradients of the s2 field
ds2dy,ds2dx = np.gradient(s2,dy,dx,edge_order=2)
ds2dydy,ds2dydx = np.gradient(ds2dy,dy,dx,edge_order=2)
ds2dxdy,ds2dxdx = np.gradient(ds2dx,dy,dx,edge_order=2)

#initialize arrays
s1_directional_derivative = np.ma.empty([ydim,xdim])
s1_concavity = np.ma.empty([ydim,xdim])
s2_directional_derivative = np.ma.empty([ydim,xdim])
s2_concavity = np.ma.empty([ydim,xdim])

#Loop through each position in the flow
for i in range(ydim):
    for j in range(xdim):
        #Check that the gridpoint is not masked
        if (ds1dx[i,j] and ds1dy[i,j] and ds1dxdy[i,j] and ds1dydy[i,j] and ds1dxdx[i,j] and ds1dydx[i,j]) is not np.ma.masked:    
            #If it's not masked compute the directional derivative and the concavity
            s1_directional_derivative[i,j] = np.dot([ds1dx[i,j],ds1dy[i,j]],Xi1[i,j,:])
            s1_concavity[i,j] = np.dot(np.dot([[ds1dxdx[i,j],ds1dxdy[i,j]],[ds1dydx[i,j],ds1dydy[i,j]]],Xi1[i,j,:]),Xi1[i,j,:])
        else:
            s1_directional_derivative[i,j] = np.ma.masked
            s1_concavity[i,j] = np.ma.masked
        
        #Check that the gridpoint is not masked
        if (ds2dx[i,j] and ds2dy[i,j] and ds2dxdy[i,j] and ds2dydy[i,j] and ds2dxdx[i,j] and ds2dydx[i,j]) is not np.ma.masked:    
            #If it's not masked compute the directional derivative and the concavity
            s2_directional_derivative[i,j] = np.dot([ds2dx[i,j],ds2dy[i,j]],Xi2[i,j,:])
            s2_concavity[i,j] = np.dot(np.dot([[ds2dxdx[i,j],ds2dxdy[i,j]],[ds2dydx[i,j],ds2dydy[i,j]]],Xi2[i,j,:]),Xi2[i,j,:])
        else:
            s2_directional_derivative[i,j] = np.ma.masked
            s2_concavity[i,j] = np.ma.masked

#Now we're going to look for iLCS using the directional derivative field,
#So apply the first condition for iLCS to the field,
#i.e. s1<0 and s2>0
s1_directional_derivative = np.ma.masked_where(s1>0,s1_directional_derivative)
s2_directional_derivative = np.ma.masked_where(s2<0,s2_directional_derivative)

#You may find additional thresholding useful if there are a lot of features in the field
#s1_directional_derivative = np.ma.masked_where(-s1<=np.percentile(-s1,85),s1_directional_derivative)
#s2_directional_derivative = np.ma.masked_where(s2<=np.percentile(s2,85),s2_directional_derivative)

#Now we apply the third condition, i.e. concavity 
s1_directional_derivative = np.ma.masked_where(s1_concavity<0,s1_directional_derivative)
s2_directional_derivative = np.ma.masked_where(s2_concavity>0,s2_directional_derivative)

#Now we can extract iLCS
plt.figure(figsize=[5+3/8,5+3/8]) 
attracting_ilcs = plt.contour(x,y,s1_directional_derivative,levels=[0],colors='blue',linewidths=4)
s2_directional_derivative = np.ma.masked_where(abs(x)>0.03,s2_directional_derivative)
repelling_ilcs = plt.contour(x,y,s2_directional_derivative,levels=[0],colors='red',linewidths=4)
i=10
plt.quiver(x[::i,::i],y[::i,::i],u[::i,::i],v[::i,::i])
plt.ylabel('y',**labelfont)
plt.xlabel('x',**labelfont)

plt.savefig('nls_iles_v1.png', transparent=False, bbox_inches='tight',dpi=300)
plt.savefig('nls_iles_v1.eps', transparent=False, bbox_inches='tight')

