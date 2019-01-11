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
from velocities import double_gyre as vel_func
#from velocities import auto_double_gyre as vel_func
from velocities import co
plt.close('all')
dimx = 301
#dimy = int(np.ceil(dimx/2.5))
dimy = int(np.ceil(dimx/2))
#dimy = 0.5*dimx
t0 = 0
tf =-0.1#days
dimt =2
t = np.linspace(t0,tf,dimt)
dt = t[1]-t[0]
x = np.linspace(0,2,dimx)
y = np.linspace(0,1,dimy)
dx = x[1]-x[0]
dy = y[1]-y[0]

#[x,y,t]=np.meshgrid(x,y,t)
[x,y]=np.meshgrid(x,y)

dudy = np.empty([dimy,dimx,dimt])
dudx = np.empty([dimy,dimx,dimt])
dvdy = np.empty([dimy,dimx,dimt])
dvdx = np.empty([dimy,dimx,dimt])
for tt in range(dimt):
    u,v = vel_func(t[tt],[x,y])
    dudy[:,:,tt],dudx[:,:,tt] = np.gradient(u,dy,dx)#,axis=(0,1))
    dvdy[:,:,tt],dvdx[:,:,tt] = np.gradient(v,dy,dx)#,axis=(0,1))


dudydt = np.gradient(dudy,dt,axis=2)
dudxdt = np.gradient(dudx,dt,axis=2)
dvdydt = np.gradient(dvdy,dt,axis=2)
dvdxdt = np.gradient(dvdx,dt,axis=2)

s1 = np.ma.empty([dimy,dimx,dimt])
b1 = np.ma.empty([dimy,dimx,dimt])
s2 = np.ma.empty([dimy,dimx,dimt])
b2 = np.ma.empty([dimy,dimx,dimt])
for t in range(dimt):
    for i in range(dimy):
        for j in range(dimx):
            if (dudxdt[i,j,t] and dudydt[i,j,t] and dvdxdt[i,j,t] and dvdydt[i,j,t]) is not np.ma.masked:    
                Grad = np.array([[dudx[i,j,t], dudy[i,j,t]], [dvdx[i,j,t], dvdy[i,j,t]]])
                Grad_dt = np.array([[dudxdt[i,j,t], dudydt[i,j,t]], [dvdxdt[i,j,t], dvdydt[i,j,t]]])
                S = 0.5*(Grad + np.transpose(Grad))
                D = (0.5*(Grad_dt + np.transpose(Grad_dt)+np.matmul(np.transpose(Grad),Grad)))
                s1[i,j,t] = np.min(np.linalg.eig(S)[0])
                b1[i,j,t] = np.min(np.linalg.eig(D)[0])
                s2[i,j,t] = np.max(np.linalg.eig(S)[0])
                b2[i,j,t] = np.max(np.linalg.eig(D)[0])
            else:
                s1[i,j,t] = np.ma.masked
                b1[i,j,t] = np.ma.masked
                s2[i,j,t] = np.ma.masked
                b2[i,j,t] = np.ma.masked
    
plt.register_cmap(name='co', data=co())


s1=s1[:,:,0]
b1=b1[:,:,0]
s2=s2[:,:,0]
b2=b2[:,:,0]
y=y[:,:]
x=x[:,:]


plt.figure()

plt.subplot(221)
#plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
plt.pcolormesh(x,y,-3600*s1)
plt.colorbar()
plt.title('s$_{1}$')

plt.subplot(222)
#plt.contourf(x,y,sigma,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.pcolormesh(x,y,-3600*b1)#,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('b$_{1}$')

plt.subplot(223)
#plt.contourf(x,y,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301))
plt.pcolormesh(x,y,3600*s2)
plt.colorbar()
plt.title('s$_{2}$')

plt.subplot(224)
#plt.contourf(x,y,sigma,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.pcolormesh(x,y,3600*b2)#,levels=np.linspace(sigma.min(axis=None),sigma.max(axis=None),301))
plt.colorbar()
#pu,pv = vel_func(0,[x[::3,::3],y[::3,::3]])
#plt.quiver(x[::3,::3],y[::3,::3],pu,pv)
plt.title('b$_{2}$')

#plt.gca().set_aspect('equal', adjustable='box', anchor='C')




'''

A = 0.1
w = np.pi*0.2
e = 0.25
a = e*np.sin(w*t)
b = 1-2*e*np.sin(w*t)
f = a*x**2+b*y
dfdx = 2*a*x+b    
u =-np.pi*A*np.sin(np.pi*f)*np.cos(y*np.pi)    
v = np.pi*A*np.cos(np.pi*f)*np.sin(y*np.pi)*dfdx


sol = sint.solve_ivp(vel_func,[t0,-1/24],y0,rtol=1e-8,atol=1e-8)
    yout2[k,:] = sol.y[:,-1]
    
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


'''
