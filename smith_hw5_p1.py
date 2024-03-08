#%% Problem 1 - Mary Smith
#defining example function to integrate over
def f(x):
    return x**2
#setting number of slices - needs to be an integer
N=int(10)
#setting example bounds a and b
a=0
b=10
#definition of h
h=(b-a)/N
#function to integrate using Simpson's rule
def I(f,N,a,b,h):
#creating the first for loop for the first summation
    for1=0
    for k in range(1,int((N/2)+1)):
        for1+=f(a+(2*k-1)*h)
#creating second for loop for the second summation
    for2=0
    for k in range(1,int(N/2)):
        for2+=f(a+2*k*h)
#Simpson's rule in full
    i=(1/3)*h*(f(a)+f(b)+4*for1+2*for2)
    return i
#print out integration of example function
print(I(f,N,a,b,h))