import numpy as np
import matplotlib.pyplot as plt
"""plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})
plt.ion()"""
#%% SINE FUNCTION

#%% defining sine Bhaskara function
def sbhas(x):
        return (4*x*(180-x))/(40500-x*(180-x))

#%% double angle function
def double(x):
    return 2*x*np.sqrt(1-x**2)

#%% marysin function
def marysin(x):
#setting variables equal to x and 2x to initiate for loop
    h=x
    hold=2*x
#counter i
    i=0
#continue this loop until machine precision says the next bhaskara approximation equals the former
    while (np.abs(sbhas(h)-sbhas(hold)))>0.0001:
#half h until loop breaks, add one to counter
        hold=h
        h=h/2
        i+=1
    s=sbhas(h)
#in range of counter above double values back to original number
    for n in range(1,i+1):
        s=double(s)
    return s,i

#%% error function
def serror(x):
#call s,i from previous part
    s,i=marysin(x)
#go one value past machine precision determined above
    s2=sbhas(x/(2**(i+1)))
#in range of counter+1, double values back to original number
    for k in range(1,i+2):
        s2=double(s2)
#gives numerator of error
    return (np.abs(s2-s))
#divides numerator by denominator of error function
def srealerror(x):
    s=marysin(x)
    r=serror(x)/s
    return r

#%% plotting the function and calculating error
#create x linspace for the function
x=np.linspace(0,1800,300)
#empty lists for the values of sine and the error
Bsine=[]
esine=[]
#Bhaskara's approximation has bounds x<180
for n in x:
#trig functions are periodic every 360 degrees, subtract to use in this function
    while n>360:
        n=n-360
#keeps sine negative when x>180
    if n>180:
        s,i=marysin(n)
        s=-s
#cases where x<180 apply function normally
    else:
        s,i=marysin(n)
#append the lists
    Bsine.append(s)
    esine.append(srealerror(n))
#plot, title, show
plt.plot(x,Bsine)
plt.title("The Sine Approximation")
plt.show()

#%% COSINE FUNCTION

#%% define marycos function
def marycos(x):
#call values from sin fn
    s,i=marysin(x)
#cosine identity
    c=np.sqrt(1-s**2)
    return c, i
   
#%% plotting the function
#create x linspace for the function
x=np.linspace(0,1800,300)
#empty list for values of cosine
Bcos=[]
for n in x:
#trig functions periodic every 360 degrees, subtract to use this function
    while n>360:
        n=n-360
#fixing function to decrease and increase over correct bounds
    if n>90 and n<270:
        c,i=marycos(n)
        c=-c
#for x<90, call cosine fn
    else:
        c,i=marycos(n)
    Bcos.append(c)
#plot, title, show
plt.plot(x,Bcos)
plt.title("The Cosine Approximation")
plt.show()

#%% TANGENT FUNCTION

#%% define marytan function
def marytan(x):
#call sine, cosine
    s,i=marysin(x)
    c=marycos(x)
#define tangent
    t=s/c
    return t,i

#%% plotting the function
#create linspace for the function
x=np.linspace(0,1800,300)
#empty list for values of tangent
Btan=[]
for n in x:
#trig functions periodic every 360 degrees, subtract to use this function
    while n>360:
        n=n-360
#fixing function to have correct behavior over correct bounds
    if n>180:
        t,i=marytan(n)
        t=-t
#for x<180, call the tangent fn
    else:
        t,i=marytan(n)
#append the empty list
    Btan.append(t)
#plot, title, show
plt.plot(x,Btan)
plt.title("The Tangent Approximation")
plt.show()

#%% PLOT OF ERROR
#plot error calculated above
#can just use sine error, because sine is only point where approximations are used
#trig identites are used elsewhere
plt.plot(x,esine)
#title, show
plt.title("Fractional Numerical Error")
plt.show()