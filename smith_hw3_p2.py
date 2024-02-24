"""import numpy as np
import matplotlib.pyplot as plt
def mad(n):
    M=float(0)
    for l in range(-n,n+1):
        for k in range(-n,n+1):
            for j in range(-n,n+1):
                M+=((-1)**(j+k+l))/(np.sqrt(j**2+k**2+l**2))
    return(M)
Ms=[]
n=range(0,100)
for i in n:
    Ms.append(mad(n))
plt.plot(n,Ms)
plt.show()"""

"""import math

nmax = 13
M = 0.
for n in range(1,nmax,2):
    for m in range(n,nmax,2):
        fac = 2
        if n == m:
            fac = 1
        M += fac * math.cosh(0.5*math.pi*math.hypot(n, m))**-2
print(-12*math.pi*M)"""
import numpy as np
nmax=15
def mad(n):
    M=0.
    for l in range(-n,n+1):
        for k in range(-n,n+1):
            for j in range(-n,n+1):
                M+=((-1)**(j+k+l))/(np.sqrt(j**2+k**2+l**2))
    return M
print (mad(5))
        
        