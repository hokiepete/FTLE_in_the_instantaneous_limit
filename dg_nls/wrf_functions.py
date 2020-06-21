# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 12:10:58 2018

A function to unstager WRF/NAM data
@author: pnola
"""

# This function will unstager velocity field data from the WRF/NAM models
def unstagger(X,axis=0):
    dim = len(X.shape)
    if dim == 1:
        return (X[0:-1] + X[1:])/2
    elif dim == 2:
        if axis == 0:
            return (X[0:-1,:] + X[1:,:])/2
        elif axis == 1:
            return (X[:,0:-1] + X[:,1:])/2
    elif dim == 3:
        if axis == 0:
            return (X[0:-1,:,:] + X[1:,:,:])/2
        elif axis == 1:
            return (X[:,0:-1,:] + X[:,1:,:])/2
        elif axis == 2:
            return (X[:,:,0:-1] + X[:,:,1:])/2
    elif dim == 4:
        if axis == 0:
            return (X[0:-1,:,:,:] + X[1:,:,:,:])/2
        elif axis == 1:
            return (X[:,0:-1,:,:] + X[:,1:,:,:])/2
        elif axis == 2:
            return (X[:,:,0:-1,:] + X[:,:,1:,:])/2
        elif axis == 3:
            return (X[:,:,:,0:-1] + X[:,:,:,1:])/2
        