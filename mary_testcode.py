# Physics 304 Final Project - Mary Smith and Joey Carol
# This program is designed to simulate the trajectory of an intercontinental ballistic missile subject to several different forces including gravity, thrust, drag, and the coriolis and centrifugal forces. It plots the trajectory on a 3D cartesian coordinate grid.

#import constants
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#%% fixing mass

#earth constants
G=6.67430e-11 #Nm^2/kg^2
Me=5.97219e24 #kg
Re=6.378e6 #m
we=7.2921159e-5 #rad/s
p0=1204.0 #kg/m^3
H=8570.0 #m
gamma=5/3
R=8.3145 #J/mol K

#rocket constant
m0=171661.0 #kg
mf=36030.0 #kg
A=2.1904 #m^2
z=1120.0e3 #m

#define vectors not based on position
w=np.array([0,0,we],float)

#define drag coefficient function
def cd(s):
    if s<0:
        s=-s
    s=np.piecewise(s,[(s>=0 and s<=1),(s>1 and s<=1.2),(s>1.2 and s<=3),(s>3)], \
                   [0.2,(0.08/0.2)*s-0.02,-(0.08/0.02)*s+(1/3),0.2])
    return s

#atmosphere temperature function
def T0(h):
    t=np.piecewise(h,[h>25000+Re,(h>11000+Re and h<25000+Re),h<11000+Re], [-131.25+0.00299*(h+Re), \
                   -56.46,15.04-0.00649*(h+Re)])
    return t

i=0

# Create a general function to compute the net force on the  rocket
def f(r,t):
# Define the different components of the r vector fed to the function storing position and velocity
    xm=r[0]
    vx=r[1]
    ym=r[2]
    vy=r[3]
    zm=r[4]
    vz=r[5]
# Create position and velocity vectors
    rm=np.array([xm,ym,zm]) # 3D cartesian postion vector
    rmmag=np.sqrt(rm[0]**2+rm[1]**2+rm[2]**2) # Rocket distance from flat Earth origin
    v=np.array([vx,vy,vz]) # 3D cartesian velocity vector
    vmag=np.sqrt(v[0]**2+v[1]**2+v[2]**2) # Magnitude of rocket velocity
# Define terms relevant to the different forces and compute the individual forces
    atemp=T0(zm)+273.15 # Atmospheric temperature in K
    Te=atemp/((1+((gamma-1)/2)*Me**2)*(1-(1/(1-(gamma-1)/2))*Me**2)) # Exhaust temperature in K
    g=np.array([-Me*G*xm/(rmmag**3),-Me*G*ym/(rmmag**3),-Me*G*zm/(rmmag**3)],float) # Computing the force due to gravity in all components
    cmag=Me*np.sqrt(gamma*R*Te) # Magnitude of the exhaust velocity
    # Resolving the exhaust velocity into its components (they are opposite the direction of the normal velocity components)
    cx=-vx*(cmag/vmag)
    cy=-vy*(cmag/vmag)
    cz=-vz*(cmag/vmag)
    c=np.array([cx,cy,cz]) # Creating an exhaust velocity vector
    delv=c*np.log(m0/mf) # Computing the total change in velocity
    delvmag=np.sqrt(delv[0]**2+delv[1]**2+delv[2]**2) # Magnitude of change in velocity
    # Compute the force due to thrust on the rocket
    T=np.array(-c[0]*(-(delv/c[0]**2)*m0*np.exp(-delv*t/c[0])) \
               -c[1]*(-(delv/c[1]**2)*m0*np.exp(-delv*t/c[1])) \
            -c[2]*(-(delv/c[2]**2)*m0*np.exp(-delv*t/c[2])),float)
    T=np.nan_to_num(T,nan=0.0) # Change all nan values to 0.0
    print(i,T)
# Ensure the mass does not decrease below the mass of the rocket
    mi=m0*np.exp(-delvmag*t/cmag) # Compute current mass
    if mi>mf:
        m=mi
    else:
        m=mf # Mass cannot go below the rocket's own mass
    p=p0*np.exp((-zm-Re)/H) # Compute atmospheric density
    
# Find the total net force on the rocket by using first order differential equations
    fxm=vx # Compute first-order derivative of position
    fvx=g[0]-2*np.cross(w,v)[0]+np.cross(np.cross(w,rm),w)[0] + \
        (T[0]/m)-(1/(2*m))*p*cd(v[0])*A*v[0]*np.linalg.norm(v) # Compute first order derivative of velocity
    fym=vy
    fvy=g[1]-2*np.cross(w,v)[1]+np.cross(np.cross(w,rm),w)[1] + \
        (T[1]/m)-(1/(2*m))*p*cd(v[1])*A*v[1]*np.linalg.norm(v)
    fzm=vz
    fvz=g[2]-2*np.cross(w,v)[2]+np.cross(np.cross(w,rm),w)[2] + \
        (T[2]/m)-(1/(2*m))*p*cd(v[2])*A*v[2]*np.linalg.norm(v)
    return np.array([fxm,fvx,fym,fvy,fzm,fvz],float)
    # Return the derivatives in an array so that they can be used for the Runge-Kutta method


#set endpoints, N and h
a=0.0
b=1000
N=2000
h=(b-a)/N

#create t, empty lists
tpoints=np.arange(a,b,h)
xmpoints=[]
vxpoints=[]
ympoints=[]
vypoints=[]
zmpoints=[]
vzpoints=[]

# Define the initial conditions of the rocket incorporating Earth's rotation
xmi=0.0 #m
vxi=Re*we*np.cos((2*np.pi*66.5)/360) #m/s
ymi=0.0 #m
vyi=Re*we*np.sin((2*np.pi*66.5)/360) #m/s
zmi=Re #m
vzi=Re*we #m/s
r=np.array([xmi,vxi,ymi,vyi,zmi,vzi],float)

# Create a for loop to run the Runge-Kutta method over all the values of t
for t in tpoints:
	# Start the loop by appending the most recent values for position and velocity stored in the r vector
    xmpoints.append(r[0])
    vxpoints.append(r[1])
    ympoints.append(r[2])
    vypoints.append(r[3])
    zmpoints.append(r[4])
    vzpoints.append(r[5])
    i+=1
    # Use the Runge-Kutta4 method to step the position and velocity according to k1, k2, k3, and k4 using the big function defined above
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6 # Update the r vector to store the new information for position and velocity

#generate the graph
plt.figure(1)
plt.plot(xmpoints,ympoints)
plt.xlabel("x position")
plt.ylabel("y position")

plt.figure(2)
plt.plot(tpoints,ympoints)
plt.xlabel("time")
plt.ylabel("y position")

plt.figure(3)
plt.plot(tpoints,xmpoints)
plt.xlabel("time")
plt.ylabel("x position")

#testing 3d plot
fig=plt.figure()
ax=fig.add_subplot(projection='3d')
ax.scatter(xmpoints,ympoints,zmpoints,'.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
"""ax.set_xlim(-1.2,-0.4)
ax.set_ylim(4,6)
ax.set_zlim(-10,-6)"""
plt.show()