import matplotlib.pyplot as plt
import matplotlib
import scipy.io as sio

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
plt.close('all')

width = 5+3/8
height = width/1.61803399
figSize = (width,height)

f = sio.loadmat('dg_plot_data.mat')
t = -f['time'][0]
cor = f['rmse_corrected'][0]
s1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,s1,'b-')
plt.plot(t,cor,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.axis('tight')
plt.xlim([0,0.7])
plt.savefig('dg_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)


f = sio.loadmat('wrf_plot_data.mat')
t = f['time'][0]
cor = f['rmse_corrected'][0]
s1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,s1,'b-')
plt.plot(t,cor,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.xlim([0,0.7])
plt.savefig('wrf_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)


f = sio.loadmat('dg3d_plot_data.mat')
t = f['time'][0]
cor = f['rmse_corrected'][0]
s1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,s1,'b-')
plt.plot(t,cor,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.xlim([0,0.7])
plt.savefig('dg3d_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)


f = sio.loadmat('abc_plot_data.mat')
t = f['time'][0]
cor = f['rmse_corrected'][0]
s1 = f['rmse_uncorrected'][0]

plt.figure(figsize=figSize)
plt.plot(t,s1,'b-')
plt.plot(t,cor,'m-')
plt.ylabel('FTLE field root mean-squared error',**labelfont)
plt.xlabel('$|T|$',**labelfont)
plt.xlim([0,0.7])
plt.savefig('abc_rmse.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)

del f
