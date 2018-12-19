# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 21:38:52 2018

@author: pnola
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

plt.figure(figsize=(8,6))
plt.plot(x,y)
plt.ylabel('Pearson correlation coefficient')
plt.xlabel('Hours')
plt.title('$s_{1}$ vs Backward-time FTLE correlation, WRF model: SE US')
plt.axis('tight')

dt = x[1]-x[0]
dy = np.gradient(y)
plt.figure(figsize=(8,6))
plt.plot(x,dy)