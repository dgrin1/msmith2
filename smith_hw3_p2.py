#import libraries
import numpy as np
import matplotlib.pyplot as plt
#defining madelung function
def mad(n):
#set M as integer 0
    M=0.
#create nested for loops in summation range for different variables
    for l in range(-n,n+1):
        for k in range(-n,n+1):
            for j in range(-n,n+1):
#exclude cases where denominator would =0 for Madelung calculation
                if l!=0 and k!=0 and j!=0:
                    M+=((-1)**(j+k+l))/(np.sqrt(j**2+k**2+l**2))
#don't calculate when denomination =0
                if l==0 and k==0 and j==0:
                    pass
    return M
#empty list for values of the Madelung function
Ms=[]
#range of x values to calculate
x=range(0,100)
#append Madelung for each value in the range
for n in x:
    Ms.append(mad(n))
#plot, title, show
plt.plot(x,Ms)
plt.title("The Madelung Constant")
plt.show()