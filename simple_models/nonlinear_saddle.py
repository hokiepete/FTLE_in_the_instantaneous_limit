# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 22:20:26 2019

@author: pnola
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

dim =7
bd =0.5
x = np.linspace(-bd,bd,dim)
x,y = np.meshgrid(x,x)
plt.close('all')
plt.quiver(x,y,x,-y-y**3)
plt.axis('equal')
plt.savefig('nl_saddle.eps', transparent=False, bbox_inches='tight')
