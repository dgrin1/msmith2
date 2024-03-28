#%% Problem 2: Exercise 8.3 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#set constants
sigma=10
rr=28
bb=8/3

#define function
def f(r,t):
    x=r[0]
    y=r[1]
    z=r[2]
    fx= sigma*(y-x)
    fy=rr*x-y-x*z
    fz=x*y-bb*z
    return np.array([fx,fy,fz],float)

#define endpoints
a=0.0
b=50.0
N=1000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
xpoints=[]
ypoints=[]
zpoints=[]

#setting initial conditions
r=np.array([0.0,1.0,0.0],float)

#apply the Runge-Kutta 4th order
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    zpoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

#generate the graph for y vs t
plt.plot(tpoints,ypoints,label="y")
plt.xlabel("t")
plt.ylabel("y")
plt.title("Lorentz Equation for y")
plt.legend(loc=0)
plt.show()

#%% part (b)

#set constants
sigma=10
rr=28
bb=8/3

#define function
def f(r,t):
    x=r[0]
    y=r[1]
    z=r[2]
    fx= sigma*(y-x)
    fy= rr*x-y-x*z
    fz= x*y-bb*z
    return np.array([fx,fy,fz],float)

#define endpoints
a=0.0
b=50.0
N=5000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
xpoints=[]
ypoints=[]
zpoints=[]

#setting initial conditions
r=np.array([0.0,1.0,0.0],float)

#apply the Runge-Kutta 4th order
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    zpoints.append(r[2])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

#generate the graph for z vs x
plt.plot(xpoints,zpoints,label="z vs x")
plt.xlabel("x")
plt.ylabel("z")
plt.title("Lorentzian z vs x")
plt.legend(loc=0)
plt.show()