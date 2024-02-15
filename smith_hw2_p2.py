###Exercise 3.2

#importing necessary libraries
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show
import numpy as np
from numpy import sin,cos,pi
import math
#%%
#part (a)
#setting range of theta as from 0 to 2pi
theta=np.linspace(0,2*np.pi,100)
#x and y parametric functions
x=2*np.cos(theta)+np.cos(2*theta)
y=2*np.sin(theta)-np.sin(2*theta)
#labeling axes, making title
plt.xlabel("x")
plt.ylabel("y")
plt.title("parametric functions y versus x")
plot(x,y)
show()
#%%
#part (b)
#setting range of theta from 0 to 10pi
theta=np.linspace(0,10*np.pi,300)
#making r function for Galilean spiral
r=theta**2
#standard x and y parametric equations
x=r*np.cos(theta)
y=r*np.sin(theta)
#labeling axes, making title
plt.xlabel("x")
plt.ylabel("y")
plt.title("Galilean spiral")
plot(x,y)
plot(x,y)
show()
#%%
#part (c)
#setting range of theta from 0 to 24pi
theta=np.linspace(0,24*np.pi,2000)
#making r function for Fey's function
r=math.e**(np.cos(theta))-2*np.cos(4*theta)+(np.sin(theta/12))**5
#standard x and y parametric equations
x=r*np.cos(theta)
y=r*np.sin(theta)
#labeling axes, making title
plt.xlabel("x")
plt.ylabel("y")
plt.title("Fey's function")
plot(x,y)
show()