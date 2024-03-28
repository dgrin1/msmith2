#%% Problem 1: Exercise 8.4 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#set constants
g=9.8 #m/s^2
l=0.1 #m
T=2*np.pi*np.sqrt(l/g) #s

#define the pendulum function
def f(r,t):
#for the vector r
    theta=r[0]
    omega=r[1]
#setting simultaneous differential equations
    ftheta=omega
    fomega=-(g/l)*np.sin(theta)
#returning values as an array
    return np.array([ftheta,fomega],float)

#set endpoints, N and h
a=0.0
b=T*10
N=1000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
thetapoints=[]
omegapoints=[]

#initial conditions
r=np.array([179*(2*np.pi)/(360),0],float)

#apply Runge-Kutta 4th order method
for t in tpoints:
    thetapoints.append(r[0])
    omegapoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

#generate the graph
plt.plot(tpoints,thetapoints, label=r"$\theta$")
plt.plot(tpoints,omegapoints, label=r"$\omega$")
plt.xlabel("t (s)")
plt.ylabel(r"$\theta$")
plt.title("The nonlinear pendulum")
plt.legend(loc=1)
plt.show()

#%% Exercise 8.5

#%% part (a)

#set constants
g=9.8 #m/s^2
l=0.1 #m
T=2*np.pi*np.sqrt(l/g) #s
c=2 #s^-1
oomega=5 #s^-1

#define the pendulum function
def f(r,t):
#for the vector r
    theta=r[0]
    omega=r[1]
#setting simultaneous differential equations
    ftheta=omega
    fomega=-(g/l)*np.sin(theta) +c*np.cos(theta)*np.sin(oomega*t)
#returning values as an array
    return np.array([ftheta,fomega],float)

#set endpoints, N and h
a=0.0
b=100
N=1000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
thetapoints=[]
omegapoints=[]

#initial conditions
r=np.array([0,0],float)

#apply Runge-Kutta 4th order method
for t in tpoints:
    thetapoints.append(r[0])
    omegapoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

#generate the graph
plt.plot(tpoints,thetapoints, label=r"$\theta$")
plt.xlabel("Time (s)")
plt.ylabel(r"$\theta$")
plt.title("The nonlinear pendulum")
plt.legend(loc=1)
plt.show()

#%% part (b)

def fb(r,t,oomega):
#for the vector r
    theta=r[0]
    omega=r[1]
#setting simultaneous differential equations
    ftheta=omega
    fomega=-(g/l)*np.sin(theta) +c*np.cos(theta)*np.sin(oomega*t)
#returning values as an array
    return np.array([ftheta,fomega],float)

#set endpoints, N and h
a=0.0
b=100
N=1000
h=(b-a)/N

#create oomegavals, empty maxes list
maxes=[]
oomegavals=np.linspace(0,100,N)

#changing big omega to see how the function responds
for oomega in oomegavals:
#empty theta, omega lists
    thetapoints=[]
    omegapoints=[]
#initial conditions
    r=np.array([0,0],float)
#apply Runge-Kutta 4th order method
    for t in tpoints:
        thetapoints.append(r[0])
        omegapoints.append(r[1])
        k1 = h*fb(r,t,oomega)
        k2 = h*fb(r+0.5*k1,t+0.5*h,oomega)
        k3 = h*fb(r+0.5*k2,t+0.5*h,oomega)
        k4 = h*fb(r+k3,t+h,oomega)
        r += (k1+2*k2+2*k3+k4)/6
    maxes.append(max(thetapoints))

#make graph
plt.plot(oomegavals,maxes, label=r"$\theta$")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude (rad)")
plt.title("Driven Pendulum Resonance")
plt.legend(loc=0)
plt.show()

