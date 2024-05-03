#Homework 11: Exercise 10.9 - Mary Smith

#importing libraries
import numpy as np
import math
import random as ran
import matplotlib.pyplot as plt

#importing latex fonts
#plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#dimensions of square array,steps
N=20
steps=250000

#create array to store spins
s=np.ones([N,N],int)

#initial values
def initial(s):
    ini=2*np.random.randint(2, size=(N,N))-1
    return ini

#apply inital values
s=initial(s)

#energy equation
def energy(n):
    e=0
    for i in range(len(n)):
        for j in range(len(n)):
            element=s[i,j]
#interactions between 4 different particles
            interactions=s[(i+1)%N,j]+s[i,(j+1)%N]+s[(i-1)%N,j]+s[i,(j-1)%N]
            e+=-element*interactions
    return e/4

#call the energy and print
en=energy(s)
print(en)

#%% part (b)

#given value of constants
J=1
T=1
beta=1

#restate 2D array to store spins
s=np.ones([N,N],int)

#initial values
def initial(s):
    ini=2*np.random.randint(2, size=(N,N))-1
    return ini

#apply inital values
s=initial(s)

#empty list for energy values, call energy
Evals=[]
E=energy(s)

#main loop
for k in range(steps):

    # Choose the particle and the move
    i = ran.randrange(N)
    j = ran.randrange(N)
    if ran.random()<0.5:
        dn = 1
        dE = J*s[i,j]*(1/2)
    else:
        dn = -1
        dE = -J*s[i,j]*(1/2)

    # Decide whether to accept the move
    if s[i,j]>1 or dn==1:
        if ran.random()<math.exp(-beta*dE):
            s[i,j] += dn
            E += dE

    Evals.append(E)

#%% part (c)

#set steps
steps=1000000

#restate 2D array to store spins
s=np.ones([N,N],int)

#initial values
def initial(s):
    ini=2*np.random.randint(2, size=(N,N))-1
    return ini

#apply inital values
s=initial(s)

#empty list for plotting, inital M value
Mvals=[]
M=0

#main loop
for k in range(steps):

    # Choose the particle and the move
    i = ran.randrange(N)
    j = ran.randrange(N)
    if ran.random()<0.5:
        dn = 1
        dM = 1*(1/2)
    else:
        dn = -1
        dM = -1*(1/2)

    # Decide whether to accept the move
    if s[i,j]>1 or dn==1:
        if ran.random()<math.exp(-beta*dM):
            s[i,j] += dn
            M += dM

    Mvals.append(M)

#create the graph
plt.plot(Mvals,label=r"M")
plt.title(r"Magnetization of the Ising Model")
plt.xlabel(r"Steps")
plt.ylabel(r"Total magnetization (M)")
plt.legend(loc=0)
plt.show()

#%% part (e)

#given value of constants
J=1
T=1
beta=1

#restate 2D array to store spins
s=np.ones([N,N],int)

#initial values
def initial(s):
    ini=2*np.random.randint(2, size=(N,N))-1
    return ini

#apply inital values
s=initial(s)

