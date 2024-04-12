#%% Problem 2: Exercise 8.14 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
#plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (a)

#define constants
m=9.1094e-31 #kg
hbar = 1.0546e-34 #J s
e = 1.6022e-19 #C
Vn=50*e #V
a=1e-11 #m
N=1000
h=(20*a)/N

#define the potential function
def V(x):
    v=(Vn*(x**2))/(a**2)
    return v

#solving the differential eqs
def f(r,x,E):
    psi=r[0]
    phi=r[1]
    fpsi=phi
    fphi = ((2*m)/(hbar**2))*(V(x)-E)*psi
    return np.array([fpsi,fphi],float)

#find psi for given E
def solve(E):
#initial conditions
   psi=0.0
   phi=1.0
   r=np.array([psi,phi],float)
#apply Runge-Kutta 4th order in for loop
   for x in np.arange(-10*a,10*a,h):
       k1 = h*f(r,x,E)
       k2 = h*f(r+0.5*k1,x+0.5*h,E)
       k3 = h*f(r+0.5*k2,x+0.5*h,E)
       k4 = h*f(r+k3,x+h,E)
       r += (k1+2*k2+2*k3+k4)/6
   return r[0]
   

#find energies using secant method


#values provided from book
E1=0.0
E2=e
psi2=solve(E1)

#define energy function to generate guesses for the secant method
def Eguessh(n):
    k=(2*Vn)/(a**2)
    omega=np.sqrt(k/m)
    en=(n+0.5)*hbar*omega
    return en

#set target accuracy
target=e/1000

#find ground state energy
while abs(E1-E2)>target:
    psi1,psi2=psi2,solve(E2)
    E1,E2=E2,E2-(psi2*(E2-E1)/(abs(psi2-psi1)))
Eground=E2/e
print("The ground state energy is E1=",Eground,"eV")

#find first excited state energy
E3=Eguessh(1)
psi3=solve(E2)

#set target accuracy
target=E2/1000

while abs(E2-E3)>target:
    psi2,psi3=psi3,solve(E3)
    E2,E3=E3,E3-(psi3*(E3-E2)/(abs(psi3-psi2)))
Efes=E3/e
print("The first excited state energy is E2=",Efes,"eV")

#find second excited state energy
E4=Eguessh(2)
psi4=solve(E3)

#set target accuracy
target=E3/1000

while abs(E3-E4)>target:
    psi3,psi4=psi4,solve(E4)
    E3,E4=E4,E4-(psi4*(E4-E3)/(abs(psi4-psi3)))
Eses=E3/e
print("The second excited state energy is E3=",Eses,"eV")

#verify spacings are equal to my accuracy

#calculate spacing between states
spacing1=Efes-Eground
spacing2=Eses-Efes

#verify spacings are equal
if abs(spacing1-spacing2)<=1/100:
    print("The energy states are equally spaced.")
if abs(spacing1-spacing2)>1/100:
    print("The energy states are not equally spaced. The spacing is",abs(spacing1-spacing2))

#%% part(b)

#define constants
m=9.1094e-31 #kg
hbar = 1.0546e-34 #J s
e = 1.6022e-19 #C
Vn=50*e #V
a=1e-11 #m
N=1000
h=(20*a)/N

#define the potential function
def Va(x):
    v=(Vn*(x**4))/(a**4)
    return v

#solving the differential eqs
def fa(r,x,E):
    psi=r[0]
    phi=r[1]
    fpsi=phi
    fphi = ((2*m)/(hbar**2))*(Va(x)-E)*psi
    return np.array([fpsi,fphi],float)

#find psi for given E
def solve(E):
#initial conditions
   psi=0.0
   phi=1.0
   r=np.array([psi,phi],float)
#apply Runge-Kutta 4th order in for loop
   for x in np.arange(-10*a,10*a,h):
       k1 = h*fa(r,x,E)
       k2 = h*fa(r+0.5*k1,x+0.5*h,E)
       k3 = h*fa(r+0.5*k2,x+0.5*h,E)
       k4 = h*fa(r+k3,x+h,E)
       r += (k1+2*k2+2*k3+k4)/6
   return r[0]
   

#find energies using secant method


#values provided from book
E1=Eguessh(1)
E2=Eguessh(2)
psi2=solve(E1)

#define energy function to generate guesses for the secant method
def Eguessh(n):
    k=(2*Vn)/(a**2)
    omega=np.sqrt(k/m)
    en=(n+(1/(2*n)))*hbar*omega
    return en

#set target accuracy
target=E2/1000

#find ground state energy
while abs(E1-E2)>target:
    psi1,psi2=psi2,solve(E2)
    E1,E2=E2,E2-(psi2*(E2-E1)/(abs(psi2-psi1)))
Eground=E2/e
print("The ground state energy is E1=",Eground,"eV")

#find first excited state energy
E3=Eguessh(3)
psi3=solve(E2)

#reset target accuracy
target=E2/1000

while abs(E2-E3)>target:
    psi2,psi3=psi3,solve(E3)
    E2,E3=E3,E3-(psi3*(E3-E2)/(abs(psi3-psi2)))
Efes=E3/e
print("The first excited state energy is E2=",Efes,"eV")

#find second excited state energy
E4=Eguessh(4)
psi4=solve(E3)

#reset target accuracy
target=E3/1000

while abs(E3-E4)>target:
    psi3,psi4=psi4,solve(E4)
    E3,E4=E4,E4-(psi4*(E4-E3)/(abs(psi4-psi3)))
Eses=E3/e
print("The second excited state energy is E3=",Eses,"eV")

#verify spacings are equal to my accuracy

#calculate spacing between states
spacing1=Efes-Eground
spacing2=Eses-Efes

#verify spacings are equal
if abs(spacing1-spacing2)<=1/100:
    print("The energy states are equally spaced.")
if abs(spacing1-spacing2)>1/100:
    print("The energy states are not equally spaced. The spacing is",abs(spacing1-spacing2))

#%% part (c) joey's way

#generate lists for plotting
psi2=[]
psi3=[]
psi4=[]
xpoints=np.arange(-5*a,5*a,h)

#initial conditions for psi and phi
r2= np.array([0.0,1.0],float)
r3= np.array([0.0,1.0],float)
r4= np.array([0.0,1.0],float)

#populate empty lists using Runge-Kutta 4th order lists

for x in xpoints:
    psi2.append(r2[0])
    k1 = h*f(r2,x,E2)
    k2 = h*f(r2+0.5*k1,x+0.5*h,E2)
    k3 = h*f(r2+0.5*k2,x+0.5*h,E2)
    k4 = h*f(r2+k3,x+h,E2)
    r2 += (k1+2*k2+2*k3+k4)/6

for x in xpoints:
    psi3.append(r3[0])
    k1 = h*f(r3,x,E3)
    k2 = h*f(r3+0.5*k1,x+0.5*h,E3)
    k3 = h*f(r3+0.5*k2,x+0.5*h,E3)
    k4 = h*f(r3+k3,x+h,E3)
    r3 += (k1+2*k2+2*k3+k4)/6

for x in xpoints:
    psi4.append(r4[0])
    k1 = h*f(r4,x,E4)
    k2 = h*f(r4+0.5*k1,x+0.5*h,E4)
    k3 = h*f(r4+0.5*k2,x+0.5*h,E4)
    k4 = h*f(r4+k3,x+h,E4)
    r4 += (k1+2*k2+2*k3+k4)/6


#normalize values using trapezoidal rule


#calculate integral for psi2
psi2=np.array(psi2,float)
#probability density
p2=0.5*h*(psi2[0]**2)
for n in range(1,psi2.size-1):
    p2+=psi2[n]**2
p2+=0.5*h*(psi2[psi2.size-1]**2)

#divide by the square root of the integral, will be graphed later
psi2n=psi2/np.sqrt(p2)

#compute new integral to check it is normalized
p2norm=0.5*h*(psi2[0]**2)
for n in range(1,psi2.size-1):
    p2norm+=psi2n[n]**2
p2norm+=0.5*h*(psi2[psi2.size-1]**2)


#calculate integral for psi3
psi3=np.array(psi3,float)
#probability density
p3=0.5*h*(psi3[0]**2)
for n in range(1,psi3.size-1):
    p3+=psi3[n]**2
p3+=0.5*h*(psi3[psi3.size-1]**2)

#divide by the square root of the integral, will be graphed later
psi3n=psi3/np.sqrt(p3)

#compute new integral to check it is normalized
p3norm=0.5*h*(psi3[0]**2)
for n in range(1,psi3.size-1):
    p3norm+=psi3n[n]**2
p3norm+=0.5*h*(psi3[psi3.size-1]**2)


#calculate integral for psi4
psi4=np.array(psi4,float)
#probability density
p4=0.5*h*(psi4[0]**2)
for n in range(1,psi4.size-1):
    p4+=psi4[n]**2
p4+=0.5*h*(psi4[psi4.size-1]**2)

#divide by the square root of the integral, will be graphed later
psi4n=psi4/np.sqrt(p4)

#compute new integral to check it is normalized
p4norm=0.5*h*(psi4[0]**2)
for n in range(1,psi4.size-1):
    p4norm+=psi4n[n]**2
p4norm+=0.5*h*(psi4[psi4.size-1]**2)

plt.plot(xpoints,psi2n,label="ground")
plt.plot(xpoints,psi3n,label="1st")
plt.plot(xpoints,psi4n,label="2nd")
plt.legend(loc=0)
plt.show()