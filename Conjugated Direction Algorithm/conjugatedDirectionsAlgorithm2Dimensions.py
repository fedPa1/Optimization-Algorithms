# @Author: Federico Pantini
# Contact: fedpa35@gmail.com

#The script is an algorithm that works in 2 variables quadratic functions to find the global minimum using the conjugated direction method

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
#so we declare 2 arrays composed by c (2x1) and Q coefficients (2x2) , Q corresponds to Hessian matrix in the quadratic form
c=[-2,5]
Q=[[10,-2],[-2,6]]

#d0 taken as arbitrary initial direction,in this case using the antigradient
d0=[-g0[0],-g0[1]]

#at this point we know that an x* critical point exists in R^2 and it's the only point of minimum 
#since Q determinant is defined positive in the quadratic form (function is strictly convex)

#1) calculate dk+1 conjugated direction given dk and Q (premultiplied)

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
    print("k=",k,"Gradient",numpy.linalg.norm(g0),"point",x,"direction",calculateDirectionMethodD1(numpy.dot(d0,Q)))
    if (numpy.linalg.norm(g0)<=1.e-7): #algorithm found the only critical point and must be global in the strictly convex function
                                   #you might want to set an arbitrary tolerance eps>0 s.t. norm(g0)<=eps since there can be cases where
                                   # norm(g0) will never converge to 0
        print("global minimum at",x,"after",k,"iterations")
        break 
    d0Q=numpy.dot(d0,Q)  
    
    alfak=-numpy.dot(g0,d0)/numpy.dot(d0Q,d0) #calculates alfak step to take from point xk, we use the same alfak in the original conjugated direction method
    
    x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]]  #calculates xk+1
                                           
    d0=calculateDirectionMethodD1(d0Q) #calculates next conjugate direction just comparing the previous one, this works because we are in R^2 and we just need 2 conjugated directions,
                                       #however, for n>2 we must compare each new direction with each already existing direction and this might not work unless
                                       # we define a new "calculate new dk+1 direction" that does the job for the dimension we need, or create a B basis of Q-conj directions in R^n
                                       # outside the loop and use for each k iteration the k direction in B
    k+=1


