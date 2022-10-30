from matplotlib import pyplot as plt
import numpy as np
import math
#import tras modulos, bibliotecas e users scripts
p0 = 5
T0 = 2

tmax = 10
step = 0.1

n_therms = 5

n = np.arange(1,n_therms)
t = np.arange(0,tmax,step)

PT = []
a0 = p0/2

for time in t:
    pt = a0
    for N in n:
        bn = -p0/(N*math.pi)*math.sin((N*2*math.pi/T0*time))
        pt += bn

    PT.append(pt)

plt.figure(1)
plt.plot(t,PT)
plt.grid()
plt.show()