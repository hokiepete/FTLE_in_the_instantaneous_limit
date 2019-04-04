import matplotlib.pyplot as plt
import matplotlib
import scipy.io as sio
import numpy as np
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['lines.linewidth']=1
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
plt.close('all')

width = 4.5#5+3/8
height = width/1.61803399
figSize = (width,height)

f = sio.loadmat('wrf_plot_corr_data.mat')
t = f['time'][0][0:90]
rmse1 = f['corr1'][0][0:90]

fig = plt.figure(figsize=figSize)
ax = fig.add_subplot(111)
ax.plot(-t,rmse1,'b-')
ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
plt.ylabel('Pearson correlation coefficient',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
plt.savefig('wrf_corr.png', transparent=False, bbox_inches='tight',pad_inches=0.03)
plt.savefig('wrf_corr.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)
