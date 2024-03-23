#%% Problem 2: Exercise 6.10 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#constants
accuracy=10.0e-6
x1=1.0
error=1.0

#loop to find value where x converges
while error>accuracy:
    x1,x2=1-np.e**(-2*x1),x1
    error=np.abs((x1-x2)/(1-(1/2)*np.e**(2*x2)))

#printing the value for c=2
print(x1)

#%% part (b)

#accuracy from the problem
accuracy=10.0e-6

#setting up lists for plotting
y=[]
cvals=np.arange(0,3,0.01)

#loop to find values of c
for c in cvals:
    x1=1.0
    error=1.0
#run loop until values converge
    while error>accuracy:
        x1,x2=1-np.e**(-c*x1),x1
        error=np.abs((x1-x2)/(1-(1/2)*np.e**(c*x2)))
#append to empty list
    y.append(x1)

#creating the graph
plt.plot(cvals,y,label="$x=1-e^{-cx}$")
plt.ylim(-0.1,1.1)
plt.title("Plot of x values versus c values")
plt.xlabel("c")
plt.ylabel("x")
plt.legend(loc=2)
plt.show()