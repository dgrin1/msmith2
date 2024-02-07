###Exercise 2.2
#%%
###part (b)
import math
#defining variables
G=6.67e-11 #gravitational constant
M=5.97e24 #Earth's mass in kg
M=int(M)
R=6371e3 #Earth's radius in meters
#creating input
T=input("Enter desired value of T in seconds: ")
T=int(T)
#creating equation for the height of the satellite
h=((G*M*T**2)/(4*(math.pi)**2))**(1/3)-R
print("The altitude of the satellite is in meters",h)

