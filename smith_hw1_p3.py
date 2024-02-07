###Exercise 2.13 part (a)
#defining the Catalan numbers function to find Cn
def cat(n):
    if n==0:
        return 1
    else:
        return ((4*n-2)/(n+1))*cat(n-1)
#calculating and printing C100
print(cat(100))