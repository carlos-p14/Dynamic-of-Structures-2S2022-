# -*- coding: utf-8 -*-
"""
@author: carlos

Código para elaboração da lista 5 de IC956 - Dinâmica das Estrturas
"""
#%% Exercício 9.1 - Dynamic Of Structures (Patrick Paultre)
from matplotlib import pyplot as plt
import numpy as np

# dados e condições iniciais do problema
m = 1 # massa [kg]
k = (2*np.pi)**2 # rigidez [N/m]
c = 0 # amortecimento [%]
final_t = 2 # tempo final [s]

u0 = 1 # deslocamento inicial [m]
v0 = 0 # velocidade inciial [m/s]

w = np.sqrt(k/m) # frequencia natural [rad/s]
T = 2*np.pi/w # periodo [s]

dt = [0.05*T, 0.1*T, 0.2*T] # passo de tempo [s]

# Resposta Analítica
A = u0
B = v0/w

for Dt in dt:
    time = np.arange(0,final_t+Dt,Dt)
    ut = [A*np.cos(w*t) + B*np.sin(w*t) for t in time]
    
    plt.figure(1)
    plt.plot(time,ut, label = f"Dt = {Dt} s")
    plt.title('Resposta Analítica')
    plt.xlabel('Tempo [s]')
    plt.ylabel('Resposta Dinâmica - u(t) [m]')
    plt.legend()
    plt.grid(True)
    plt.show()

# Newmark Average Acceleration Method
#beta = 1/4
#gama = 1/2

# Newmark Linear Acceleration Method
beta = 1/6
gama = 1/2

an = -k*u0/m # aceleração inicial [m/s²]

un = u0
vn = v0

for Dt in dt:
    a0 = 1/(beta*Dt**2)
    a1 = gama/(beta*Dt)
    a2 = 1/(beta*Dt)
    a3 = 1/(2*beta)-1
    a4 = gama/beta -1
    a5 = (gama/(2*beta)-1)*Dt
    a6 = (1-gama)*Dt
    a7 = gama*Dt
    
    k_ = k+a0*m+a1*c
    time = np.arange(0,final_t+Dt,Dt)
    
    ut2 = np.zeros_like(time)
    at = np.zeros_like(time)
    vt = np.zeros_like(time)
    
    ut2[0] = u0
    at[0] = an
    vt[0] = vn
    
    for i in range(len(time)-1):
        pn = m*(a0*ut2[i] + a2*vt[i] + a3*at[i]) + c*(a1*ut2[i]+a4*vt[i]+a5*at[i])
        
        ut2[i+1] = pn/k_
        at[i+1] = a0*(ut2[i+1]-ut2[i])-a2*vt[i]-a3*at[i]
        vt[i+1] = vt[i]+a6*at[i]+a7*at[i+1]
        
    plt.figure(2)
    plt.plot(time,ut2, label = f"Dt = {Dt} s")
    plt.title('Newmark Linear Acceleration Method - Algoritmo Próprio')
    plt.xlabel('Tempo [s]')
    plt.ylabel('Resposta Dinâmica - u(t) [m]')
    plt.legend()
    plt.grid(True)
    plt.show()
    
# Newmark Avarage Acceleration Method
# Resposta do LAS
ut3 = [
       [1.000E+000,
        9.518E-001,
        8.120E-001,
        5.939E-001,
        3.187E-001,
        1.273E-002,
        -2.945E-001,
        -5.733E-001,
        -7.969E-001,
        -9.437E-001,
        -9.997E-001,
        -9.593E-001,
        -8.266E-001,
        -6.142E-001,
        -3.427E-001,
        -3.818E-002,
        2.700E-001,
        5.522E-001,
        7.812E-001,
        9.350E-001,
        9.987E-001,
        9.662E-001,
        8.407E-001,
        6.341E-001,
        3.665E-001,
        6.361E-002,
        -2.454E-001,
        -5.308E-001,
        -7.651E-001,
        -9.257E-001,
        -9.971E-001,
        -9.725E-001,
        -8.542E-001,
        -6.536E-001,
        -3.901E-001,
        -8.900E-002,
        2.207E-001,
        5.091E-001,
        7.484E-001,
        9.157E-001,
        9.948E-001],
               
       [1.000E+000,
        8.203E-001,
        3.459E-001,
        -2.528E-001,
        -7.607E-001,
        -9.952E-001,
        -8.722E-001,
        -4.357E-001,
        1.573E-001,
        6.938E-001,
        9.810E-001,
        9.157E-001,
        5.214E-001,
        -6.027E-002,
        -6.203E-001,
        -9.574E-001,
        -9.505E-001,
        -6.021E-001,
        -3.732E-002,
        5.409E-001,
        9.247E-001],
       
       [1.000E+000,
        4.339E-001,
        -6.234E-001,
        -9.750E-001,
        -2.227E-001,
        7.817E-001,
        9.011E-001,
        2.320E-004,
        -9.009E-001,
        -7.820E-001,
        2.222E-001]
       ]

for u, Dt in zip(ut3,dt):
    time = np.arange(0,final_t+Dt,Dt)
    
    plt.figure(3)
    plt.plot(time,u, label = f"Dt = {Dt} s")
    plt.title('Newmark Average Acceleration Method - LAS')
    plt.xlabel('Tempo [s]')
    plt.ylabel('Resposta Dinâmica - u(t) [m]')
    plt.legend()
    plt.grid(True)
    plt.show()
    
ut4 = [
       [1.000E+000,
        9.515E-001,
        8.105E-001,
        5.909E-001,
        3.139E-001,
        6.389E-003,
        -3.017E-001,
        -5.805E-001,
        -8.030E-001,
        -9.474E-001,
        -9.999E-001,
        -9.553E-001,
        -8.179E-001,
        -6.011E-001,
        -3.260E-001,
        -1.916E-002,
        2.895E-001,
        5.701E-001,
        7.953E-001,
        9.433E-001,
        9.997E-001,
        9.590E-001,
        8.252E-001,
        6.113E-001,
        3.380E-001,
        3.194E-002,
        -2.773E-001,
        -5.595E-001,
        -7.875E-001,
        -9.390E-001,
        -9.993E-001,
        -9.625E-001,
        -8.324E-001,
        -6.214E-001,
        -3.500E-001,
        -4.471E-002,
        2.650E-001,
        5.489E-001,
        7.795E-001,
        9.345E-001,
        9.987E-001],
       
       [1.000E+000,
        8.148E-001,
        3.278E-001,
        -2.807E-001,
        -7.851E-001,
        -9.988E-001,
        -8.425E-001,
        -3.741E-001,
        2.328E-001,
        7.535E-001,
        9.951E-001,
        8.681E-001,
        4.195E-001,
        -1.844E-001,
        -7.201E-001,
        -9.890E-001,
        -8.916E-001,
        -4.639E-001,
        1.356E-001,
        6.849E-001,
        9.805E-001],
       
       [1.000E+000,
        3.749E-001,
        -7.188E-001,
        -9.140E-001,
        3.346E-002,
        9.391E-001,
        6.707E-001,
        -4.361E-001,
        -9.978E-001,
        -3.121E-001,
        7.637E-001]
       ]

for u, Dt in zip(ut4,dt):
    time = np.arange(0,final_t+Dt,Dt)
    
    plt.figure(4)
    plt.plot(time,u, label = f"Dt = {Dt} s")
    plt.title('Newmark Linear Acceleration Method - LAS')
    plt.xlabel('Tempo [s]')
    plt.ylabel('Resposta Dinâmica - u(t) [m]')
    plt.legend()
    plt.grid(True)
    plt.show()