# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 20:29:21 2018

@author: pnola
"""

from velocities import abc
import numpy as np
t=0
x = np.linspace(0,2*np.pi,9)
x,y,z=np.meshgrid(x,x,x)
u,v,w=abc(t,[x,y,z])

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.quiver(x,y,z,u,v,w,normalize=True)