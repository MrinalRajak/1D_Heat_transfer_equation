

#1D Heat transfer equation
#simplest expression of the computational algorithm using Forward Euler
#method and explicit python loops.
#For this method F<=0.5 for stability to solve 1D diffusion equation

import numpy as np
import matplotlib.pyplot as plt
L=0.8
T=0.1
alpha=0.001
dt=0.001
F=0.018 #constant number
Nt=int(round(T/float(dt)))
print(Nt)

#mesh points in time (Nt+1) in number i.e 0 to Nt
t=np.linspace(0.0001,Nt*dt,Nt+1)
dx=np.sqrt(alpha*dt/F)
Nx=int(round(L/dx))
print(Nx)
#mesh points in space (Nx+1) in number i.e 0 to Nx
x=np.linspace(0,L,Nx+1)
s='t%f'
#Total no.of space steps neccessary is (Nx+1) as it runs from 0 to Nx
u=np.zeros(Nx+1)
u_n=np.zeros(Nx+1)
#initial distribution at t=0
sig=0.1

def I(x):
    return (1./np.sqrt(2*np.pi*sig))*np.exp(-(x-L/2)**2/(2*sig**2))

#set initial condition at  u(x,0)=I(x)
#initial distribution defined over the entire range 0 to Nx

for i in range(0,Nx+1):
    u_n[i]=I(x[i])
for t in range(0,Nt):
    #compute u at inner mesh points runs between 1 to (Nx-1),boundary
    #condition take care of u at boundaries,i.e u[0] and u[Nx]
    for i in range(1,Nx):
        u[i]=u_n[i]+F*(u_n[i-1]-2*u_n[1]+u_n[1+i])
    #insert boundary condition predefined according to the problem have
    #both edges kept at zeros as initial conditions in space and time
    u[0]=0
    u[Nx]=0
    #set variables before next time step
    u_n,u=u,u_n
    if(t%1==0):
        plt.plot(x,u,label=s%t)
        plt.legend()
plt.show()


























