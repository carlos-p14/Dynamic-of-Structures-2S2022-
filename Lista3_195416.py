# -*- coding: utf-8 -*-
"""
@author: carlos
01/10/2022

Código para elaboração da terceira lista de IC956 - Dinâmica das Estruturas
"""
#%% Módulos usados
from matplotlib import pyplot as plt
import numpy as np
import math
#%% Exercício 1
# Dados:
m = 1000 # massa [kg]
k = 157913.4 # rigidez [N/m]
nc = 5 # nº de ciclos
tf = 30 # tempo total [s]
E = 5/100 # amortecimento [%]
p0 = 1000 # força inicial [N]
w_ = math.pi # frequência de excitação [rad/s]

step = 0.05 # passo [s]

# Cálculo das váriaveis do problema:
w = np.sqrt(k/m) # frequência natural [rad/s]
beta = w_/w # razão de frequência
T = 2*math.pi/w_ # periodo [s]

tc = np.arange(0,nc*T+step,step) # tempo de atuação da força
tl = np.arange(nc*T,tf+step,step) # tempo de vibração livre

# constantes da solução
G1 = (p0/k)*(-2*E*beta)/((1-beta**2)**2+(2*E*beta)**2)
G2 = (p0/k)*(1-beta**2)/((1-beta**2)**2+(2*E*beta)**2)

# resposta dinâmica
# deslocamento 
up = [G1*np.cos(w_*t) + G2*np.sin(w_*t) for t in tc]
# velocidade
vp = [-w_*G1*np.sin(w_*t) + w_*G2*np.cos(w_*t) for t in tc]

# resposta livre
wd = w*np.sqrt(1-E**2) # frequencia amortecida

# condições de contorno
u0 = up[-1]
v0 = vp[-1]

A = u0
B = (E*w*u0 + v0)/wd

ut = [np.exp(-E*w*(t-nc*T))*(A*np.cos(wd*(t-nc*T))+B*np.sin(wd*(t-nc*T))) for t in tl]

# plotando resposta
plt.figure(1)
plt.plot(tc,up, label='Resposta Forçada')
plt.plot(tl,ut, label='Resposta Livre')
plt.grid(True)
plt.legend()
plt.xlabel('Tempo [s]')
plt.ylabel('Deslocamento u(t) [m]')
plt.title('Exercício 1 - Resposta Dinâmica')
plt.show()

#%% Exercício 2 - 4.1 Dynamics of Structures
ust = 0.05 # deslocamento estático [m] (p0/k)
E = 0 # taxa de amortecimento [%]
beta = [0.2, 0.9, 1.1, 1.8, 3.0] # razões de frequência

G = [ust/np.sqrt((1-B**2)**2+(2*E*B)**2) for B in beta]

plt.figure(2)
plt.plot(beta,G,marker='o')
plt.grid(True)
plt.xlabel('Beta')
plt.ylabel('|G|')
plt.title('Amplitude da Resposta Dinâmica')
plt.show()

#%% Exercício 3 - 4.8 Dynamics of Structures
m = 0.5e3/9.81 # mass [kg]
A = 5e-3 # amplitude de deslocamento [m]
wsup = 1800*2*math.pi/60 # frequência [rad/s]
R = 80/100 

# cálculo do beta e da frequência do aparelho
beta = np.sqrt((2-R)/(1-R))
wap = wsup/beta

# cálculo da rigidez
Kt = wap**2*m
ki = Kt/4 # considerando molas em paralelo

#%% Exercício 4 - 5.1 Dynamics of Structures
# Dados do exercício
w_ = [0,
      3.05,
      6.10,
      9.15,
      12.20,
      15.25,
      18.30,
      21.35,
      24.40,
      27.45,
      30.50,
      33.55,
      36.60,
      39.64,
      42.69,
      45.74,
      48.79,
      51.84,
      54.89,
      57.94,
      60.99,
      64.04
      ]

u_p = [0,
       0.103,
       0.425,
       1.007,
       1.938,
       3.385,
       5.692,
       9.642,
       17.410,
       37.004,
       78.413,
       48.549,
       31.443,
       24.250,
       20.448,
       18.129,
       16.581,
       15.482,
       14.665,
       14.037,
       13.541,
       13.141       
       ]

# encontrando maior frequência
u_pmax = max(u_p)
w = w_[u_p.index(u_pmax)] 

# e a frequência de meia largura de banda
usqrt2 = u_pmax/np.sqrt(2)

# interpolando valores manualmente (numpy tem a função inter)
wa = -1*((u_pmax - usqrt2)/(u_pmax - u_p[u_p.index(u_pmax)-1])*(w-w_[u_p.index(u_pmax)-1])-w) 
wb = -1*((u_pmax - usqrt2)/(u_pmax - u_p[u_p.index(u_pmax)+1])*(w-w_[u_p.index(u_pmax)+1])-w) 

# encontrando a taxa de amortecimento
E = (wb-wa)/(2*w)

# plotando gráfico da resposta
plt.figure(3)
plt.plot(w_,u_p, marker='o', label = r'$\frac{\overline{u_0}}{p_0}$')
plt.grid(True)
plt.plot([0, max(w_)], [usqrt2, usqrt2], ls='--', label = r'$\frac{u_{max}}{\sqrt{2}}$' f' = {usqrt2:.2f} g/kN')
plt.plot([wa,wa], [0,usqrt2], color='black', ls='--', label = r'$\overline{\omega}_a$' f' = {wa:.2f} rad/s')
plt.plot([wb,wb], [0,usqrt2], color='red', ls='--', label = r'$\overline{\omega}_b$' f' = {wb:.2f} rad/s')
plt.xlabel(r'$\overline{\omega}$' f' [rad/s]')
plt.ylabel(r'$\frac{\overline{u_0}}{p_0}\times10^3$' f' [g/kN]')
plt.title('Ensaio de Vibração')
plt.legend()
plt.show()

plt.figure(4)
plt.plot(w_,u_p, label = r'$\frac{\overline{u_0}}{p_0}$')
plt.grid(True)
plt.scatter(w, u_pmax, color='r', label=r'$\omega$' f' = {w:.2f} rad/s')
plt.xlabel(r'$\overline{\omega}$' f' [rad/s]')
plt.ylabel(r'$\frac{\overline{u_0}}{p_0}\times10^3$' f' [g/kN]')
plt.title('Ensaio de Vibração')
plt.legend()
plt.show()

#%% Exercício 5 - 5.4 Dynamics of Structures
tests = [20, 50, 80]

# Dados do exercício
ED = [9371, 33987, 64978]
Fmax = [6443, 10283, 13514]
Fmin = [-6465, -10286, -13570]
uFmax = [1.855, 4.741, 7.671]
uFmin = [-1.752, -4.623, -7.514]
umax = [1.885, 4.793, 7.728]
umin = [-1.780, -4.668, -7.570]

h = 10 # mm
f = 4 # Hz
w_ = 2*math.pi*f

k = []
Eeq = []
ceq = []

# calculando a rigidez k, taxa de amortecimento e coeficiente de amortecimento
# para cada nível do ensaio
for i in range(len(tests)):
    ki = (Fmax[i] - Fmin[i])/(uFmax[i] - uFmin[i])
    Eeqi = ED[i]/(2*math.pi*uFmax[i]**2*ki)
    ceqi = ED[i]/(math.pi*w_*uFmax[i]**2)
    
    k.append(ki)
    Eeq.append(Eeqi)
    ceq.append(ceqi)
