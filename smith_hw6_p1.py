#%% Problem 1: Exercise 6.9 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
#%% part (b)

#defining constants
L=5.0e-10 #m
a=1.6022e-18 #J
M=9.1094e-31 #kg
hbar=1.05457e-34 #J s

#defing function
def H(m,n):
    m=int(m)
    n=int(n)
#determining if numbers given are even or odd
    rm=int(m%2)
    rn=int(n%2)
#gives values for integrals for when m doesn't =n
    if m!=n:
#if they are both even or both odd
        if rm==rn:
            h=0
#if one is even and one is odd
        elif rm!=rn:
            h=-((2*a)/(L**2))*(((2*L)/(np.pi))**2)*((m*n)/((m**2)-(n**2))**2)
#integral values when m=n
    elif m==n:
        h=(((hbar**2)*(n**2)*(np.pi**2))/(2*M*(L**2)))+(a/2)
    return h


#%% part (c)

#creating empty matrix
Hvalsc=np.zeros([10,10],float)

#calculating H for each position in matrix by indices
for i in range(10):
    for j in range(10):
        m=i+1
        n=j+1
        Hvalsc[i,j]=H(m,n)

#calculating Eigenvalues
Ec=np.linalg.eigvalsh(Hvalsc)

#converting to eV
eigc=Ec*6.242e18 #eV
print(eigc)

#%% part (d)

#creating empty matrix
Hvalsd=np.zeros([100,100],float)

#calculating H for each position in matrix by summing through indices
for i in range(100):
    for j in range(100):
        m=i+1
        n=j+1
        Hvalsd[i,j]=H(m,n)

#calculating Eigenvalues
Ed=np.linalg.eigvalsh(Hvalsd)

#converting to eV
eigd=Ed*6.242e18 #eV
print(eigd)

#%% part (e)

#calling eigenvalues and eigenvectors
eigvals,eigvecs=np.linalg.eigh(Hvalsd)

#setting values of x to plot
x=np.linspace(0,L,100)

#psi function
def ps(state,x):
    p=0
#summing over range of all values of eigenvectors
    for j in range(eigvecs[state].size):
        p+=eigvecs[state,j]*np.sin((np.pi*x*j)/L)
    return p

#creating empty lists for plotting
p0=[]
p1=[]
p2=[]

#appending list for range of eigenvectors for probability densities
for i in range(100):
    p0.append((ps(0,x[i]))**2)
    p1.append((ps(1,x[i]))**2)
    p2.append((ps(2,x[i]))**2)

#creating the plot
plt.plot(x,p0,label="Ground state")
plt.plot(x,p1,label="1st excited state")
plt.plot(x,p2,label="2nd excited state")
plt.xlabel("x")
plt.ylabel("$|\psi(x)|^{2}$")
plt.legend(loc=1)
plt.show()