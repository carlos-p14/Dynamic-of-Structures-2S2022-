# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 19:39:50 2022

@author: carlo
"""
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def Newmark(Type, u0, v0, k, c, m, p, dt, td):
    if Type == 'Average':
    # Newmark Average Acceleration Method
        beta = 1/4
        gama = 1/2
    elif Type == 'Linear':
    # Newmark Linear Acceleration Method
        beta = 1/6
        gama = 1/2
    
    an = (p[0] - c*v0 -k*u0)/m # aceleração inicial [m/s²]
    
    un = u0
    vn = v0
    
    pn = np.zeros(int(td/dt))
    
    pn[0:len(p)] = p

    a0 = 1/(beta*dt**2)
    a1 = gama/(beta*dt)
    a2 = 1/(beta*dt)
    a3 = 1/(2*beta)-1
    a4 = gama/beta -1
    a5 = (gama/(2*beta)-1)*dt
    a6 = (1-gama)*dt
    a7 = gama*dt
    
    k_ = k+a0*m+a1*c
    time = np.arange(0,td+dt,dt)
    
    ut = np.zeros_like(time)
    at = np.zeros_like(time)
    vt = np.zeros_like(time)
    
    ut[0] = u0
    at[0] = an
    vt[0] = vn
    
    for i in range(len(time)-1):
        P = pn[i] + m*(a0*ut[i] + a2*vt[i] + a3*at[i]) + c*(a1*ut[i]+a4*vt[i]+a5*at[i])
        
        ut[i+1] = P/k_
        at[i+1] = a0*(ut[i+1]-ut[i])-a2*vt[i]-a3*at[i]
        vt[i+1] = vt[i]+a6*at[i]+a7*at[i+1]

    return ut, vt, at    
    
Type = 'Linear'
td = 3
dt = 0.05

m = 15e3
k = 2000e3
w = np.sqrt(k/m)
xi = 0
c = xi*2*m*w

v0 = 0
u0 = 0
p0 = 0
pn = [0,
        28000,
        51000,
        71000,
        86000,
        97000,
        100000,
        90000,
        71000,
        50000,
        29000,
        11000]

Ut, Vt, At = Newmark(Type, u0, v0, k, c, m, pn, dt, td)
time = np.arange(0,td+dt,dt)

plt.figure(1)
plt.plot(time,Ut, label = f"Dt = {dt} s")
plt.title('Newmark Linear Acceleration Method')
plt.xlabel('Tempo [s]')
plt.ylabel('Resposta Dinâmica - u(t) [m]')
plt.legend()
plt.grid(True)
plt.show()

d = {'At [m/s²]': At, 'Vt [m/s]': Vt, 'Ut [m]': Ut}
df = pd.DataFrame(data=d)

df.to_csv('Ex95.csv',index=False)