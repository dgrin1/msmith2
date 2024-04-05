#%% Problem 1: Exercise 8.10 - Mary Smith

#importing libraries
import numpy as np
import matplotlib.pyplot as plt

#importing latex fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#%% part (b)

#define constants
G=6.6743e-11 #N m^2 kg^-2
M=1.989e30 #kg


#define function
def f(r,t):
#for vector r
    x=r[0]
    vx=r[1]
    y=r[2]
    vy=r[3]
#define radius
    rad=np.sqrt((x**2)+(y**2))
#setting simultaneous differential equations
    fx=vx
    fvx=-G*M*(x/rad**3)
    fy=vy
    fvy=-G*M*(y/rad**3)
#return values as an array
    return np.array([fx,fvx,fy,fvy],float)

#test periods
p= 50 #yr
T=p*365.25*24*3600 #s

#set endpoints, N and fixed h
a=0.0
b=T*2.0
N=100000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
xpoints=[]
vxpoints=[]
ypoints=[]
vypoints=[]

#initial conditions
r=np.array([4e12,0,0,500],float)

#apply Runge-Kutta 4th order method
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

#generate the graph
plt.plot(xpoints,ypoints,label=r"Comet position")
plt.xlabel(r"x Coordinate (m)")
plt.ylabel(r"y Coordinate (m)")
plt.title(r"Comet position using a fixed step size")
plt.legend(loc=0)
plt.show()

#%% part (c)

#define constants
G=6.6743e-11 #N m^2 kg^-2
M=1.989e30 #kg


#define function
def f(r,t):
#for vector r
    x=r[0]
    vx=r[1]
    y=r[2]
    vy=r[3]
#define radius
    rad=np.sqrt((x**2)+(y**2))
#setting simultaneous differential equations
    fx=vx
    fvx=-G*M*(x/rad**3)
    fy=vy
    fvy=-G*M*(y/rad**3)
#return values as an array
    return np.array([fx,fvx,fy,fvy],float)

#test periods
p= 50 #yr
T=p*365.25*24*3600 #s

#set endpoints and N, inital step size
a=0.0
b=T*2.0
N=100000
h0=(b-a)/N
h=h0
i=0
t=a
#accuracy
delta=3.16880878e-5 #m/s

#initial conditions
r=np.array([4e12,0,0,500],float)

#create empty lists
tpoints=[]
xpoints=[]
vxpoints=[]
ypoints=[]
vypoints=[]
dxarr=[]


#grab first values
xpoints.append(r[0])
vxpoints.append(r[1])
ypoints.append(r[2])
vypoints.append(r[3])
tpoints.append(t)
dxarr.append(0)

#applying adaptive Runge-Kutta 4th order
while t<b: #condition on time
	i+=1 #increment to count steps
	#RK4 used throughout

#arrays for different guesses
	r1=np.copy(r)
	r2=np.copy(r)
    
#one large step
	k1 = 2.*h*f(r2,t)
	k2 = 2.*h*f(r2+0.5*k1,t+h)
	k3 = 2.*h*f(r2+0.5*k2,t+h)
	k4 = 2.*h*f(r2+k3,t+2*h)
	r2 += (k1+2*k2+2*k3+k4)/6
	# Two small steps
	
	k1 = h*f(r1,t)
	k2 = h*f(r1+0.5*k1,t+0.5*h)
	k3 = h*f(r1+0.5*k2,t+0.5*h)
	k4 = h*f(r1+k3,t+h)
	r1 += (k1+2*k2+2*k3+k4)/6
	
	k1 = h*f(r1,t)
	k2 = h*f(r1+0.5*k1,t+0.5*h)
	k3 = h*f(r1+0.5*k2,t+0.5*h)
	k4 = h*f(r1+k3,t+h)
	r1 += (k1+2*k2+2*k3+k4)/6

	#Calculate rho value and assess error
	dx=r1[0]-r2[0]
	rho=30.*h*delta/abs(dx)
#	print(t,rho,delta,h,dx,r1[0],r2[0])
	#adjust step size
	if rho>=1.0001:
		t+=2*h
		h*=min(rho**0.25,1.001)
		r=np.copy(r1)
		xpoints.append(r[0])
		vxpoints.append(r[1])
		ypoints.append(r[2])
		vypoints.append(r[3])
		tpoints.append(t)
		dxarr.append(dx)
	else:
		h*=rho**0.25
		k1 = 2.*h*f(r,t)
		k2 = 2.*h*f(r+0.5*k1,t+h)
		k3 = 2.*h*f(r+0.5*k2,t+h)
		k4 = 2.*h*f(r+k3,t+2*h)
		r += (k1+2*k2+2*k3+k4)/6
		
		xpoints.append(r[0])
		vxpoints.append(r[1])
		ypoints.append(r[2])
		vypoints.append(r[3])
		tpoints.append(t)
		dxarr.append(dx)
        #k's use r1, r at bottom

#generate the graph
plt.plot(xpoints,ypoints,label=r"Comet position")
plt.xlabel(r"x Coordinate (m)")
plt.ylabel(r"y Coordinate (m)")
plt.title(r"Comet position using an adaptive step size")
plt.legend(loc=0)
plt.show()

#%% part (d)

#define constants
G=6.6743e-11 #N m^2 kg^-2
M=1.989e30 #kg


#define function
def f(r,t):
#for vector r
    x=r[0]
    vx=r[1]
    y=r[2]
    vy=r[3]
#define radius
    rad=np.sqrt((x**2)+(y**2))
#setting simultaneous differential equations
    fx=vx
    fvx=-G*M*(x/rad**3)
    fy=vy
    fvy=-G*M*(y/rad**3)
#return values as an array
    return np.array([fx,fvx,fy,fvy],float)

#test periods
p= 50 #yr
T=p*365.25*24*3600 #s

#set endpoints and N, inital step size
a=0.0
b=T
N=100000
h0=(b-a)/N
h=h0
i=0
t=a
#accuracy
delta=3.16880878e-5 #m/s

#initial conditions
r=np.array([4e12,0,0,500],float)

#create empty lists
tpoints=[]
xpoints=[]
vxpoints=[]
ypoints=[]
vypoints=[]
dxarr=[]


#grab first values
xpoints.append(r[0])
vxpoints.append(r[1])
ypoints.append(r[2])
vypoints.append(r[3])
tpoints.append(t)
dxarr.append(0)

#applying adaptive Runge-Kutta 4th order
while t<b: #condition on time
	i+=1 #increment to count steps
	#RK4 used throughout

#arrays for different guesses
	r1=np.copy(r)
	r2=np.copy(r)
    
#one large step
	k1 = 2.*h*f(r2,t)
	k2 = 2.*h*f(r2+0.5*k1,t+h)
	k3 = 2.*h*f(r2+0.5*k2,t+h)
	k4 = 2.*h*f(r2+k3,t+2*h)
	r2 += (k1+2*k2+2*k3+k4)/6
	# Two small steps
	
	k1 = h*f(r1,t)
	k2 = h*f(r1+0.5*k1,t+0.5*h)
	k3 = h*f(r1+0.5*k2,t+0.5*h)
	k4 = h*f(r1+k3,t+h)
	r1 += (k1+2*k2+2*k3+k4)/6
	
	k1 = h*f(r1,t)
	k2 = h*f(r1+0.5*k1,t+0.5*h)
	k3 = h*f(r1+0.5*k2,t+0.5*h)
	k4 = h*f(r1+k3,t+h)
	r1 += (k1+2*k2+2*k3+k4)/6

	#Calculate rho value and assess error
	dx=r1[0]-r2[0]
	rho=30.*h*delta/abs(dx)
#	print(t,rho,delta,h,dx,r1[0],r2[0])
	#adjust step size
	if rho>=1.0001:
		t+=2*h
		h*=min(rho**0.25,1.001)
		r=np.copy(r1)
		xpoints.append(r[0])
		vxpoints.append(r[1])
		ypoints.append(r[2])
		vypoints.append(r[3])
		tpoints.append(t)
		dxarr.append(dx)
	else:
		h*=rho**0.25
		k1 = 2.*h*f(r,t)
		k2 = 2.*h*f(r+0.5*k1,t+h)
		k3 = 2.*h*f(r+0.5*k2,t+h)
		k4 = 2.*h*f(r+k3,t+2*h)
		r += (k1+2*k2+2*k3+k4)/6
		
		k1 = 2.*h*f(r,t)
		k2 = 2.*h*f(r+0.5*k1,t+h)
		k3 = 2.*h*f(r+0.5*k2,t+h)
		k4 = 2.*h*f(r+k3,t+2*h)
		r += (k1+2*k2+2*k3+k4)/6
		
		xpoints.append(r[0])
		vxpoints.append(r[1])
		ypoints.append(r[2])
		vypoints.append(r[3])
		tpoints.append(t)
		dxarr.append(dx)
        #k's use r1, r at bottom

#generate the graph
plt.plot(xpoints,ypoints,'.',markersize=0.1,label=r"Comet position")
plt.xlabel(r"x Coordinate (m)")
plt.ylabel(r"y Coordinate (m)")
plt.title(r"Comet position using an adaptive step size")
plt.legend(loc=0)
plt.show()
