import numpy as np
import matplotlib.pyplot as plt
#%% SINE FUNCTION

#%% defining sine Bhaskara function
def sbhas(x):
        return (4*x*(180-x))/(40500-x*(180-x))

#%% double angle
def double(x):
    return 2*x*np.sqrt(1-x**2)

#%% marysin function
def marysin(x):
    h=x
    hold=2*x
    i=0
    while (np.abs(sbhas(h)-sbhas(hold)))>0.01:
        hold=h
        h=h/2
        i+=1
    s=sbhas(h)
    for n in range(1,i+1):
        s=double(s)
    return s,i

#%% finding error
def serror(x):
    s,i=marysin(x)
    s2=sbhas(x/(2**(i+1)))
    for k in range(1,i+2):
        s2=double(s2)
    return (np.abs(s2-s))
def srealerror(x):
    s=marysin(x)
    r=serror(x)/s
    return r

#%% plotting the function
x=np.linspace(0,1800,300)
Bsine=[]
esine=[]
for n in x:
    while n>360:
        n=n-360
    if n>180:
        s,i=marysin(n)
        s=-s
    else:
        s,i=marysin(n)
    Bsine.append(s)
    esine.append(srealerror(n))
plt.plot(x,Bsine)
plt.plot(x,esine)
plt.title("the sine function")
plt.show()

#%% COSINE FUNCTION

#%% defining cosine Bhaskara function
def cbhas(x):
    return (32400-4*x**2)/(32400+x**2)

#%% marycos function
def marycos(x):
    h=x
    hold=2*x
    i=0
    while (np.abs(cbhas(h)-cbhas(hold)))>0.01:
        hold=h
        h=h/2
        i+=1
    s=cbhas(h)
    for n in range(1,i+1):
        s=double(s)
    return s,i

#%% finding error
def cerror(x):
    s,i=marycos(x)
    s2=cbhas(x/(2**(i+1)))
    for k in range(1,i+2):
        s2=double(s2)
    return (np.abs(s2-s))
def crealerror(x):
    s=marycos(x)
    r=cerror(x)/s
    return r

#%% plotting the function
x=np.linspace(0,1800,300)
Bcos=[]
ecos=[]
for n in x:
    while n>360:
        n=n-360
    if n>180:
        s,i=marycos(n)
        s=-s
    else:
        s,i=marycos(n)
    Bcos.append(s)
    ecos.append(crealerror(n))
plt.plot(x,Bcos)
plt.plot(x,ecos)
plt.title("the cosine function")
plt.show()

#%% TANGENT FUNCTION

#%% defining tangent Bhaskara function
def tbhas(x):
    if cbhas(x)==0:
        t=0
    else:
        t=sbhas(x)/cbhas(x)
    return t

#%% marytan function
def marytan(x):
    h=x
    hold=2*x
    i=0
    while (np.abs(tbhas(h)-tbhas(hold)))>0.01:
        hold=h
        h=h/2
        i+=1
    s=tbhas(h)
    for n in range(1,i+1):
        s=double(s)
    return s,i

#%% finding error
def terror(x):
    s,i=marytan(x)
    s2=tbhas(x/(2**(i+1)))
    for k in range(1,i+2):
        s2=double(s2)
    return (np.abs(s2-s))
def trealerror(x):
    s=marytan(x)
    r=terror(x)/s
    return r

#%% plotting the function
x=np.linspace(0,1800,300)
Btan=[]
etan=[]
for n in x:
    while n>360:
        n=n-360
    if n>180:
        s,i=marytan(n)
        s=-s
    else:
        s,i=marytan(n)
    Btan.append(s)
    etan.append(trealerror(n))
plt.plot(x,Btan)
plt.plot(x,etan)
plt.show()

