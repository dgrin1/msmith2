#%% Exercise 5.13 - Mary Smith
#importing libraries
import math
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#defining the Hermite polynomial
def H(n,x):
#given values of H
    if n==0: 
        return 1
    elif n==1:
        return 2*x
#to keep n values greater than or equal to 0
    elif n<0:
        pass
#Hermite polynomial
    else:
        h=2*x*H(n-1,x)-2*(n-1)*H(n-2,x)
    return h

#defining the wavefunction
def ps(n,x):
#keeping n values greater than or equal to 0
    if n<0:
        pass
#the wavefunction
    else:
        p=(1/(math.sqrt((2**(n))*(math.factorial(n))*(math.sqrt(math.pi)))))*(np.exp(((-x**2)/2))*(H(n,x)))
    return p

#domain of x
x=np.linspace(-4,4,100)

#plotting the wavefunction for given values of n
plt.plot(x,ps(0,x),label="n=0")
plt.plot(x,ps(1,x),label="n=1")
plt.plot(x,ps(2,x),label="n=2")
plt.plot(x,ps(3,x),label="n=3")

#title, labels, legend, show
plt.title("Wavefunction of the Quantum Harmonic Oscillator at various values of n")
plt.xlabel("Values of x")
plt.ylabel("Wavefunction $\psi_n(x)$")
plt.legend(loc=1)
plt.show()

#%% part (b)

#setting domain of x to given values
x=np.linspace(-10,10,500)

#creating the plot
plt.plot(x,ps(30,x),label="n=30")

#title, labels, legend, show
plt.title("Wavefunction of the Quantum Harmonic Oscillator")
plt.xlabel("Values of x")
plt.ylabel("Wavefunction $\psi_n(x)$")
plt.legend(loc=1)
plt.show()

#%% part (c)

#import Newman's code
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
#%%

#defining root-mean-square function
def delx(n):
#integrand function
    def exint(n,x):
        exs=(x**2)*(np.abs(ps(n,x))**2)
        return exs
#Newman's Gaussian quadrature function
    x,w=gaussxwab(50,-122,122)
    s=0.0
    N=50
#for loop to use Newman's code
    for k in range(N):
        s+=w[k]*exint(n,x[k])
#square root of for loop result
    dels=np.sqrt(s)
    return dels

#find the uncertainty for n=5
print(delx(5))
