# @Author: Federico Pantini
# Contact: fedpa35@gmail.com

#The script is an algorithm that works in 2 variables quadratic functions to find the global minimum
#using the conjugated Gradient method

#References: https://www.stat.cmu.edu/~ryantibs/convexopt-F13/scribes/lec10.pdf

#Notation:
#^t indicates the transposed operator in arrays
#* indicates multiplication sign except in x* indicating a point in R^2

import numpy
import math
import random

#initial random point x in R^2 
x=[random.random()*1000,random.random()*1000]

#function (f0 function,g0 gradient,h0 hessian matrix) f0 must match 1/2*x^tQx+c^tx form
#we can input any function as long it matches the required form with Q nxn with determinant defined positive
f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
#h0=[[10,-2],[-2,6]]

#we will need to express the function in a quadratic form 1/2*x^tQx+c^tx
#so we declare 2 arrays composed by c (2x1) and Q coefficients (2x2)
c=[-2,5]
Q=[[10,-2],[-2,6]]

#d0 taken using the antigradient
d0=[-g0[0],-g0[1]]

#at this point we know that an x* critical point exists in R^2 and it's the only point of minimum 
#since Q determinant is defined positive in the quadratic form (function is strictly convex)

#we now have all the data needed to start the algorithm and can procede with the iteration

k=0
while (k<=10**4):
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    d0Q=numpy.dot(d0,Q) 
    print("k=",k,"Gradient",numpy.linalg.norm(g0),"point",x)
    if (numpy.linalg.norm(g0)==0): #algorithm found the only critical point and must be global in the strictly convex function
        print("global minimum at",x,"after",k,"iterations")
        break
    alfak=-numpy.dot(g0,d0)/numpy.dot(d0Q,d0) #calculates alfak step to take from point xk
    x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]]
    b1=numpy.dot(g0,d0Q)/numpy.dot(d0Q,d0)
    d0=[g0[0]+numpy.dot(b1,d0)[0],g0[1]+numpy.dot(b1,d0)[1]]
    k+=1


