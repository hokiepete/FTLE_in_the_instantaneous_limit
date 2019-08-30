import numpy as np
import scipy.integrate as sint
import matplotlib.pyplot as plt
plt.close()

def 

#Initalize variables
xdim = 61
ydim = 31
x = np.linspace(0,2,xdim)
y = np.linspace(0,1,ydim)
flow_map = np.empty([len(y)*len(x),2])
dx = x[1]-x[0]
dy = y[1]-y[0]
x,y = np.meshgrid(x,y)

t0 = 0
tf = 10
backward_time = [tf,t0] # AT t_0
forward_time = [t0,tf] # AT t_0

#define velocity field
from velocity_fields import double_gyre

#integrate velocity field
for k,y0 in enumerate(zip(x.ravel(),y.ravel())):
    sol = sint.solve_ivp(double_gyre,backward_time,y0,rtol=1e-8,atol=1e-8)
    flow_map[k,:] = sol.y[:,-1]

fx,fy = zip(*flow_map)
fx = np.reshape(fx,[ydim,xdim])
fy = np.reshape(fy,[ydim,xdim])

#Calculate flow map gradients
dfxdy,dfxdx = np.gradient(fx,dy,dx,edge_order=2)
dfydy,dfydx = np.gradient(fy,dy,dx,edge_order=2)
del fx,fy

sigma = np.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        #Calculate Cauchy-Green tensor, C
        JF = np.array([[dfxdx[i,j],dfxdy[i,j]],[dfydx[i,j],dfydy[i,j]]])
        C = np.dot(JF.T, JF)
        
        #Calculate FTLE, sigma
        lam=np.max(np.linalg.eig(C)[0])
        if lam>=1:
            sigma[i,j]=1.0/(2.0*abs(tf-t0))*np.log(lam)
        else:
            sigma[i,j]=0
            
del dfxdy,dfxdx,dfydy,dfydx

#Plot results
plt.pcolormesh(x,y,sigma)
plt.colorbar()
plt.gca().set_aspect('equal', adjustable='box', anchor='C')
pu,pv = double_gyre(0,[x[::3,::3].ravel(),y[::3,::3].ravel()])
plt.quiver(x[::3,::3].ravel(),y[::3,::3].ravel(),pu,pv)
