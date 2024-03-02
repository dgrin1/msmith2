#%%Code taken from Mark Newman, <mejn@umich.edu>
from numpy import ones,copy,cos,tan,pi,linspace

def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w
#%% Exercise 5.9

#%% part (a)
#import libararies
import numpy as np
import scipy.constants as sc
#set values given in problem
V=0.001 #m^3
rho=6.022e28 #m^-3
theta=428. #K
#define cv(T) function
def cv(T):
#internal integral function
    def cint(x):
        c=((x**(4))*np.exp(x))/((np.exp(x)-1)**2)
        return c
#using Newman's Gaussian quadrature
    x,w=gaussxwab(50,0,theta/T)
#50 points
    N=50
    s=0.0
#find value of integral w/ gauss
    for k in range(N):
        s+=w[k]*cint(x[k])
#finds value of full cv(T) function
    cv=9*V*rho*sc.Boltzmann*((T/theta)**3)*s
    return cv

#%% part (b)
#import needed library
import matplotlib.pyplot as plt
#given domain of T
T=linspace(5,500,1000)
#empty list of values for cv(T)
cvs=[]
#appends list for the number of T values in domain
for k in range(len(T)):
    cvs.append(cv(T[k]))
#latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
#title, labels
plt.title("Debye's Heat Capacity of a Solid")
plt.xlabel("Temperature (K)")
plt.ylabel("Heat capacity (J/K)")
#plot and show
plt.plot(T,cvs)
plt.show()