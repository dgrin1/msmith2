###Exercise 3.1
#%%
#part (a)
#necessary libraries
from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show
#importing data
data = loadtxt("sunspots.txt",float)
x=data[:,0]
y=data[:,1]
#labeling axes, title
plt.xlabel("Months since January 1749")
plt.ylabel("Number of sunspots")
plt.title("Sunspot number versus Months since January 1749")
#legend
plot(x,y,label='sunspot number')
plt.legend(loc=2)
show()
#%%
#part (b)
from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show
#importing data
data = loadtxt("sunspots.txt",float)
#data to the 1000th row for x,y columns
x=data[:1000,0]
y=data[:1000,1]
#labeling axes, title
plt.xlabel("Months since January 1749")
plt.ylabel("Number of sunspots")
plt.title("First 1000 sunspots versus Months since January 1749")
#legend
plot(x,y,label='sunspot number')
plt.legend(loc=2)
show()
#%%
#part (c)
from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show
#importing data
data = loadtxt("sunspots.txt",float)
#data to the 1000th row for x,y columns
x=data[:1000,0]
y=data[:1000,1]
#creating the sum in the average
r=5
m=0
s=0
for m in range(-r,r+m):
    s+=y+m
#running average as influenced by sum above
avg=(1/(2*r+1))*s
#labeling axes, title
plt.xlabel("Months since January 1749")
plt.ylabel("Number of sunspots")
plt.title("First 1000 sunspots versus Months since January 1749 with running average of data")
plot(x,y,label='sunspot number')
#legend
plot(x,avg,'m:',label='running average of the data')
plt.legend(loc=2)
show()
