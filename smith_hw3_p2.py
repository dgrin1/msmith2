import numpy as np
import matplotlib.pyplot as plt
def mad(n):
    M=0.
    for l in range(-n,n+1):
        for k in range(-n,n+1):
            for j in range(-n,n+1):
                if l!=0 and k!=0 and j!=0:
                    M+=((-1)**(j+k+l))/(np.sqrt(j**2+k**2+l**2))
                if l==0 and k==0 and j==0:
                    pass
    return M
Ms=[]
x=range(0,100)
for n in x:
    Ms.append(mad(n))
plt.plot(x,Ms)
plt.show()
        
        