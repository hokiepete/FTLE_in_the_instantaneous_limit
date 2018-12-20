from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy import integrate
from velocities import double_gyre as vel_func

dimx = 301
#dimy = int(np.ceil(dimx/2.5))
dimy = int(np.ceil(dimx/2))
dimt = 101
#dimy = 0.5*dimx
t0 = 0
tf = -2 #days

x = np.linspace(0,2,dimx)
y = np.linspace(0,1,dimy)
dx = x[1]-x[0]
dy = y[1]-y[0]

want_time = np.linspace(t0,tf,dimt)
fu = np.empty([dimt,dimy,dimx])
fv = np.empty([dimt,dimy,dimx])
for i, yy in enumerate(y):
    for j, xx in enumerate(x):
        print('integrating x={0}, y={1}'.format(xx,yy))
        sol = integrate.solve_ivp(vel_func,[t0,tf],(xx,yy),rtol=1e-8,atol=1e-8)
        x_flow = sol.y[0,:]
        y_flow = sol.y[1,:]
        int_time = sol.t
        f = interpolate.interp1d(int_time,x_flow)
        fu[:,i,j] = f(want_time)
        f = interpolate.interp1d(int_time,y_flow)
        fv[:,i,j] = f(want_time)
        
dfudy,dfudx = np.gradient(fu,dy,dx,axis=(1,2),edge_order=2)
dfvdy,dfvdx = np.gradient(fv,dy,dx,axis=(1,2),edge_order=2)
del fu,fv
#JF = np.empty([dimy,dimx])
ftle = np.empty([dimt,dimy,dimx])
for t in range(dimt):
    for i in range(dimy):
        for j in range(dimx):
            print('FTLE calculation x={0}, y={1}, t={2}'.format(x[j],y[i],want_time[t]))
            JF = np.array([[dfudx[t,i,j],dfudy[t,i,j]],[dfvdx[t,i,j],dfvdy[t,i,j]]])
            C = np.dot(JF.T, JF)
            lam=np.max(np.linalg.eig(C)[0])
            if lam>=1:
                ftle[t,i,j]=1.0/(2.0*abs(tf-t0))*np.log(lam)
            else:
                ftle[t,i,j]=0
                #sigma[i,j]=1/(2.0*abs(tf-t0))*np.log(lam)
del dfudy,dfudx,dfvdy,dfvdx

x,y = np.meshgrid(x,y)

u,v = vel_func(0,[x,y])

dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)


s1 = np.ma.empty([dimy,dimx])
J = np.array([[0, 1], [-1, 0]])
for i in range(dimy):
    for j in range(dimx):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = np.min(np.linalg.eig(S)[0])

        else:
            #s1[i,j] = np.nan #np.ma.masked
            s1[i,j] =  999999

s1 = np.ma.masked_where(s1==999999,s1)            
s1 = s1 - s1.max(axis=None)
s1 = s1/s1.min(axis=None)
s1 = s1.filled(np.nan)
data = [s1.ravel()]
name = ['s1']
import matplotlib.pyplot as plt
for tt,time in enumerate(want_time):
    if tt == 0:
        continue
    
    #f = ftle[timelen-1-tt,:,:] - ftle[timelen-1-tt,:,:].min(axis=None)
    f = ftle[tt,:,:] - ftle[tt,:,:].min(axis=None)
    f = f/f.max(axis=None)
    f =  np.ma.filled(f,np.nan)
    data.append(f.ravel())
    name.append('{0:2.3f}'.format(time))
    plt.figure()
    plt.pcolormesh(x,y,f)
    plt.savefig('ftle_{:}.png'.format(tt))
    plt.close('all')

#Alldata = pd.DataFrame(np.transpose([s1.ravel(),ftle_2hr.ravel(),ftle_4hr.ravel(),ftle_6hr.ravel()]),columns=['s1','2hr','4hr','6hr'])
Alldata = pd.DataFrame(np.transpose(data),columns=name)
Alldata.corr().to_csv('Correlation_and_stats_dg.csv',mode='w')
B=Alldata.describe()
B.to_csv('Correlation_and_stats_dg.csv',mode='a')

plt.close('all')
plt.figure(1)
plt.pcolormesh(s1)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(f)
plt.colorbar()





