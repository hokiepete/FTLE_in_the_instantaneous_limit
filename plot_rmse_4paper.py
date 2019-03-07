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

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])#,axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

width = 4.5#5+3/8
height = width/1.61803399
figSize = (width,height)
'''
f = sio.loadmat('dg_plot_data.mat')
end = 105
t = -f['time'][0][0:end]
rmse2 = f['rmse2'][0][0:end]
rmse1 = f['rmse1'][0][0:end]
rmse3 = f['rmse3'][0][0:end]

fig = plt.figure(figsize=figSize)
ax = fig.add_subplot(111)
ax.plot(t,rmse1,'b-')
ax.plot(t,rmse2,'m-')
ax.plot(t,rmse3,'k-')
ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
plt.ylim([0,0.175])
ax1=add_subplot_axes(ax,rect = [0.09,0.54,0.37,0.37])
#ax1.plot(t[0:-1],rmse1[0:1],'r-')

ax1.plot(t[0:5],rmse2[0:5],'m-')
ax1.plot(t[0:5],rmse3[0:5],'k-')

ax1.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
#plt.xlim([0,0.02])
plt.savefig('dg_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)
plt.savefig('dg_rmse.png', transparent=False, bbox_inches='tight',pad_inches=0.03)

'''
f = sio.loadmat('wrf_plot_data.mat')
t = f['time'][0]
#rmse2 = f['rmse2'][0][0:-60]
#rmse1 = f['rmse1'][0][0:-60]
rmse1 = f['corr1'][0]
#rmse3 = f['rmse3'][0][0:-60]

fig = plt.figure(figsize=figSize)
ax = fig.add_subplot(111)
ax.plot(t,rmse1,'b-')
#ax.plot(-t,rmse2,'m-')
#ax.plot(-t,rmse3,'k-')
ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
#plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.ylabel('Pearson correlation coefficient',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
"""
ax1=add_subplot_axes(ax,rect = [0.09,0.54,0.37,0.37])
#ax1.plot(t[0:-1],rmse1[0:1],'r-')
ax1.plot(-t[0:10],rmse2[0:10],'m-')
ax1.plot(-t[0:10],rmse3[0:10],'k-')
ax1.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
"""
#plt.xlim([0,400])
plt.savefig('wrf_rmse.png', transparent=False, bbox_inches='tight',pad_inches=0.03)
plt.savefig('wrf_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)
#'''
'''
f = sio.loadmat('dg3d_plot_data.mat')
t = f['time'][0]
rmse2 = f['rmse_corrected'][0]
rmse1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,rmse1,'b-')
plt.plot(t,rmse2,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
plt.xlim([0,0.7])
plt.savefig('dg3d_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)


f = sio.loadmat('abc_plot_data.mat')
t = f['time'][0]
rmse2 = f['rmse_corrected'][0]
rmse1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,rmse1,'b-')
plt.plot(t,rmse2,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
plt.xlim([0,0.7])
plt.savefig('abc_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)

del f
#'''
#plt.close('all')