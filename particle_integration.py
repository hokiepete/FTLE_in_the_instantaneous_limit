# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 11:28:01 2018

@author: pnola
def kevrekidis(t,Y, eps=0.01):
  return np.array([-Y[0]-Y[1]+2, 1/eps*(Y[0]**3-Y[1])])

def ex11(t,Y):
  return -np.array([np.tanh(Y[0]**2/4)+Y[0],Y[0]+2*Y[1]])

def rotHoop(t,Y, eps=0.1, gamma=2.3):
  return np.array([Y[1], 1/eps*(np.sin(Y[0])*(gamma*np.cos(Y[0])-1)-Y[1])])


def verhulst(t,Y, eps=0.01):
  return [1, 1/eps*(Y[0]*Y[1]-Y[1]**2)]


"""

import numpy as np
import scipy.integrate as sint
import matplotlib.pyplot as plt
#from velocities import van_der_pol_oscillator as vel_func
from velocities import verhulst as vel_func
from velocities import co
plt.close('all')
dimx = 101
#dimy = int(np.ceil(dimx/2.5))
#dimy = int(np.ceil(dimx/2))
dimy = dimx
t0 = 0
tf = -5 #days

#x = np.linspace(0,20000,dimx)
#y = np.linspace(-4000,4000,dimy)
#x = np.linspace(-10,10,dimx)
#y = np.linspace(-10,10,dimy)
x = np.linspace(-1,1,dimx)
y = np.linspace(-1,1,dimy)
dx = x[1]-x[0]
dy = y[1]-y[0]
yout=np.empty([len(y)*len(x),2])
yout2=np.empty([len(y)*len(x),2])

x,y = np.meshgrid(x,y)
for k,y0 in enumerate(zip(x.ravel(),y.ravel())):
    print(k)
    sol = sint.solve_ivp(vel_func,[t0,tf],y0,rtol=1e-8,atol=1e-8)
    yout[k,:] = sol.y[:,-1]
    sol = sint.solve_ivp(vel_func,[t0,-1/24],y0,rtol=1e-8,atol=1e-8)
    yout2[k,:] = sol.y[:,-1]

fu,fv = zip(*yout)
fu = np.reshape(fu,[dimy,dimx])
fv = np.reshape(fv,[dimy,dimx])
dfudy,dfudx = np.gradient(fu,dy,dx,edge_order=2)
dfvdy,dfvdx = np.gradient(fv,dy,dx,edge_order=2)
del fu,fv
#JF = np.empty([dimy,dimx])
sigma = np.empty([dimy,dimx])
for i in range(dimy):
    for j in range(dimx):
        JF = np.array([[dfudx[i,j],dfudy[i,j]],[dfvdx[i,j],dfvdy[i,j]]])
        C = np.dot(JF.T, JF)
        lam=np.max(np.linalg.eig(C)[0])
        if lam>=1:
            sigma[i,j]=1.0/(2.0*abs(tf-t0))*np.log(lam)
        else:
            sigma[i,j]=0
            #sigma[i,j]=1/(2.0*abs(tf-t0))*np.log(lam)
del dfudy,dfudx,dfvdy,dfvdx

fu,fv = zip(*yout2)
fu = np.reshape(fu,[dimy,dimx])
fv = np.reshape(fv,[dimy,dimx])
dfudy,dfudx = np.gradient(fu,dy,dx,edge_order=2)
dfvdy,dfvdx = np.gradient(fv,dy,dx,edge_order=2)
del fu,fv
#JF = np.empty([dimy,dimx])
sigma2 = np.empty([dimy,dimx])
for i in range(dimy):
    for j in range(dimx):
        JF = np.array([[dfudx[i,j],dfudy[i,j]],[dfvdx[i,j],dfvdy[i,j]]])
        C = np.dot(JF.T, JF)
        lam=np.max(np.linalg.eig(C)[0])
        if lam>=1:
            sigma2[i,j]=1.0/(2.0*abs(-1/24-t0))*np.log(lam)
        else:
            sigma2[i,j]=0
            #sigma[i,j]=1/(2.0*abs(tf-t0))*np.log(lam)
del dfudy,dfudx,dfvdy,dfvdx


plt.figure(1)
plt.subplot(321)
plt.contourf(x,y,sigma2,levels=np.linspace(sigma2.min(axis=None),sigma2.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('$\\sigma$, integration time = -1hr')

plt.subplot(322)
plt.contourf(x,y,sigma,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('$\\sigma$, integration time  t = -5days')

#plt.gca().set_aspect('equal', adjustable='box', anchor='C')

u,v = vel_func(0,[x,y])

dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)


s1 = np.ma.empty([dimy,dimx])
rhodot = np.ma.empty([dimy,dimx])
nudot = np.ma.empty([dimy,dimx])
J = np.array([[0, 1], [-1, 0]])
for i in range(dimy):
    for j in range(dimx):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            if np.dot(Utemp, Utemp) != 0:
                Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
                rhodot[i, j] = 3600*np.dot(Utemp, np.dot(np.dot(np.transpose(J), np.dot(S, J)), Utemp))/np.dot(Utemp, Utemp)
                nudot[i, j] = 3600*np.dot(Utemp,np.dot(np.trace(S)*np.identity(2) - 2*S, Utemp))/np.dot(Utemp, Utemp)
            else:
                rhodot[i,j] = np.ma.masked
                nudot[i,j] = np.ma.masked
                s1[i,j] = np.ma.masked

        else:
            rhodot[i,j] = np.ma.masked
            nudot[i,j] = np.ma.masked
            s1[i,j] = np.ma.masked
plt.register_cmap(name='co', data=co())
rholevel = np.min(np.abs([rhodot.min(axis=None),rhodot.max(axis=None)]))
nulevel = np.min(np.abs([nudot.min(axis=None),nudot.max(axis=None)]))
plt.subplot(323)
plt.contourf(x,y,rhodot,levels=np.linspace(rhodot.min(axis=None),rhodot.max(axis=None),301),\
             vmin=-2/3*rholevel,vmax=2/3*rholevel,cmap='co')
plt.colorbar()#format=ticker.FuncFormatter(mf.fmt))
plt.title('$\dot{\\rho}$')

plt.subplot(324)
plt.contourf(x,y,nudot,levels=np.linspace(nudot.min(axis=None),nudot.max(axis=None),301),\
             vmin=-2/3*nulevel,vmax=2/3*nulevel,cmap='PiYG')
plt.colorbar()
plt.title('$\dot{\\nu}$')

rdim = rhodot.shape
for i in range(rdim[0]):
    for j in range(rdim[1]):
        if (rhodot[i,j]>0 and nudot[i,j]>0) or (rhodot[i,j]<0 and nudot[i,j]<0):
            'nothing'
            #rhodot = rhodot[(rhodot>0 and nudot>0) or (rhodot<0 and nudot<0)]
        else:
            rhodot[i,j] = np.ma.masked

plt.subplot(325)
plt.contourf(x,y,rhodot,levels=np.linspace(rhodot.min(axis=None),rhodot.max(axis=None),301),\
             vmin=-2/3*rholevel,vmax=2/3*rholevel,cmap='co')
plt.colorbar()#format=ticker.FuncFormatter(mf.fmt))
plt.title('$\dot{\\rho}$ filtered by $\dot{\\nu}$')


plt.subplot(326)
plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
plt.colorbar()
plt.title('s$_{-1}$')

fig=plt.figure(4)
plt.plot(sigma2.ravel())