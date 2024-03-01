#import libraries
import numpy as np
import matplotlib.pyplot as plt
#%% The first function
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
#set fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
#plot, title, show
plt.plot(x,Ms)
plt.title("The Madelung Constant")
plt.show()
#%% The second function
#import needed library
import math
#set M,n values
M=0.
n=1
#while loop for the summation
while True:
    Mold=M
#the summation, over odd values
    for m in range(1,n+1,2):
        s=2
#for when n=m in summation index
        if n==m:
            s=1
#creating value we are summing over
        M+=s*math.cosh(0.5*math.pi*np.sqrt(m**2+n**2))**-2
#break when machine precision met
    if abs(M-Mold)<10e-4:
        break
#keep sum over odd values
    n+=2
M2=12*np.pi*M
print(M2)
