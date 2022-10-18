#%% Script para elaboração da lista 4 de IC956
# Dinâmica das Estruturas I - 2S2022 
# Carlos Henrique Chama Puga - 195416

# Importando os  Módulos usados
from matplotlib import pyplot as plt
import numpy as np
import math

#%% Exercício 1: Dynamics of Structures (6.3)
p0 = 5 # Definindo p0 [N]
T0 = 2 # e o período [s]
w_ = 2*math.pi/T0 # frequência de excitação [rad/s]

# settings do tempo [s] 
step = 0.01
tmax = 10+step

# Número de termos considerados na série de Fourier
n_therms = 100

# Criando listas dentro dos intervalos
n = np.arange(1,n_therms)
t = np.arange(0,tmax,step)

PT = []

# calculando a0
# a0, an e bn foram integrados usando o Mathematica
a0 = p0/2

for time in t:
    # calculando pt
    pt = a0
    for N in n:
        # e bn (nesse caso, an = 0)
        bn = -p0/(N*math.pi)*math.sin((N*w_*time))
        pt += bn

    PT.append(pt)

# Plotando a série de Fourier para esse caso
plt.figure(1)
plt.plot(t,PT)
plt.grid()
plt.title('Série de Fourier')
plt.xlabel('tempo [s]')
plt.ylabel('p(t) [N]')
plt.show()

#%% Exercício 2: Dynamics of Structures 7.1
