import numpy as np
import matplotlib.pyplot as plt
#%% Exercise 5.4

#%% part (a)
#defining j(m,x)
def j(m,x):
#defining internal integral function of theta
    def jint(theta):
        j=(1/np.pi)*np.cos(m*theta-x*np.sin(theta))
        return j
#number of bins, starting points, h
    N=1000
    a=0
    b=np.pi
    h=(b-a)/N
#breaking the summations in Simpson's rule into 2 for loops
    for1=0
#summing over 1st summation
    for k in range(1,int((N/2)+1)):
        for1+=jint(a+(2*k-1)*h)
    for2=0
#summing over second summation
    for k in range(1,int(N/2)):
        for2+=jint(a+2*k*h)
#calculates integral for entire Bessel function
    i=(1/3)*h*(jint(a)+jint(b)+4*for1+2*for2)
    return i
#creating domain of x
x=np.linspace(0,20,100)
#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
#plotting for various values of m
plt.plot(x,j(0,x),label='$J_{0}(x)$')
plt.plot(x,j(1,x),label='$J_{1}(x)$')
plt.plot(x,j(2,x),label='$J_{2}(x)$')
#making labels, legend, title
plt.xlabel("Distance (m)")
plt.ylabel("Bessel Function $J_{m}(x)$")
plt.legend(loc=1)
plt.title("The Bessel Function for Various m Values")
plt.show()
#%% part (b)
#setting lambda and k values
lambd=500*10.e-9 #m
k=(2*np.pi)/lambd #1/m
#creating domains for x and y
x=np.linspace(-1.e-6,1.e-6,200) #m
y=np.linspace(-1.e-6,1.e-6,200) #m
#vectorizing x and y
xvec,yvec=np.meshgrid(x,y)
#calculating r
r=np.sqrt(xvec**2+yvec**2) #m
#calculating intensity
inte=((j(1,k*r))/(k*r))**2
#latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
#label axes, title
plt.xlabel("Distance (m)")
plt.ylabel("Distance (m)")
plt.title("Intensity of the Circular Diffraction Pattern")
#plot and show
plt.imshow(inte)
plt.colorbar()
plt.show()
