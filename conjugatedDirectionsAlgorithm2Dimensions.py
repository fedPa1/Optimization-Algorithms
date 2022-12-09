#Author: Federico Pantini fedpa35@gmail.com

#The script is an algorithm that works in 2 variables quadratic functions to find the global minimum
#using the conjugated directions method

#Why use conjugated directions method?
#(1)Accelarate the convergence rate of steepest decent.
#(2)Avoiding the high computational cost Newton’s method

#link to the article https://www.stat.cmu.edu/~ryantibs/convexopt-F13/scribes/lec10.pdf

import numpy
import math
import random
import numdifftools as nd

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
Q=[[10,-4],[0,6]]

#d0 taken as arbitrary initial direction,in this case using the antigradient

d0=[-g0[0],-g0[1]]

#at this point we know that an x* critical point exists in R^2 and it's the only point of minimum 
#since Q determinant is defined positive in the quadratic form (function is strictly convex)

#1) calculate dk+1 direction given dk and Q (already multiplied)

def calculateDirectionMethodD1(d0Q):
    d1=[0,0]
    # case 0=0 we can choose components of d1 arbitrary , we choose d1=(1,1) in this case
    if(d0Q[0]==0 and d0Q[1]==0):
        d1=[1,1]

    # case d0Q[0]=0 we can choose x arbitrary (1 in this case) and  y=0
    elif(d0Q[0]==0):
        d1=[1,0]

    # case d0Q[1]=0 we can choose y arbitrary (1 in this case) and  x=0
    elif(d0Q[1]==0):
        d1=[0,1]

    #case d0Q components are both not 0 , can be demonstrated by simple algebric steps that
    #d1[0]=(-d0Q[1]d1[1])/d0Q[0], so we choose d1[1] arbitrary 1 in this case and d1[0]=-d0Q[1]/d0Q[0]
    else:
        d1[1]=1
        d1[0]=(-d0Q[1]*d1[1])/d0Q[0]
    return d1

#2) we now have all the data needed to start the algorithm and can procede with the iteration

k=0
while (k<=10**4):
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    d0Q=numpy.dot(d0,Q) #matrix multiplication between d0^t and Q
    print("k=",k,"Gradient",numpy.linalg.norm(g0),"point",x,"direction",calculateDirectionMethodD1(d0Q))
    if (numpy.linalg.norm(g0)==0): #algorithm found the only critical point and must be global in the strictly convex function
        print("global minimum at",x,"after",k,"iterations")
        break
    alfak=-numpy.dot(g0,d0)/numpy.dot(d0Q,d0) #calculate alfak step to take from point xk
    x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]]
    d0=calculateDirectionMethodD1(d0Q) #calculate d1 based on d0
    k+=1


