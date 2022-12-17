# @Author: Federico Pantini 
# email:fedpa35@gmail.com

# The script computes the global minimum point of a 2 variable function (that admits a global minimum) using a modification of the Newton Method


import numpy
import random
import armijoRule

#initial random point x in R^2 

x=[random.random()*1000,random.random()*1000]

#function (f0 function,g0 gradient,Q hessian matrix), we can input any 2 variables function

f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]  #your input function,Q,g0,d0, in case you change the function remember to change the fAlpha in armijorule.py too
d0=[-g0[0],-g0[1]]
Q=[[10,-2],[-2,6]]
alfak=0
sk=[0,0]

k=0

print("Initial point",x)

while k<=10**7:
    f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    if (numpy.linalg.norm(g0)<=1.e-6):
        print("global minimum at",x,"after",k,"iterations")
        break
    DetQ=Q[0][0]*Q[1][1]-Q[0][1]*Q[1][0]

    c1=1.e-7 #arbitrary constant >0 we set the degree of tolerance choosing Newton method over the antigradient
         # to calculate next direction, low c1 => more iterations with Newton method
         #you can play around yourself with c1,c2 constants to see different results

    if (abs(DetQ)>c1): #DetQ must be !=0 

        sk=[-numpy.dot(numpy.linalg.inv(Q),g0)[0],-numpy.dot(numpy.linalg.inv(Q),g0)[1]]

        c2=1.e-7 #arbitrary constant >0 , this sets how much the direction should be descending ,
             # low c2 => more tolerance towards less descending direction, however a low tolerance doesn't mean
             # that it will need less iterations, it is used primarily to avoid potential too weak directions calculated
             #in the next alphak step

        if(numpy.dot(g0,sk)<=-c2):
            print("Newthon method used")
            d0=[sk[0],sk[1]]
            x=[x[0]+d0[0],x[1]+d0[1]] #its x+d0 because we already put a minus in sk (xk+1=xk-g0/Q)
        
        #the direction sk couldn't match expectations so we use antigradient path
        else:
            print("Gradient method used")
            d0=[-g0[0],-g0[1]]
            alfak=armijoRule.armijo(f0,x,d0,g0)
            x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]]

    #the direction sk couldn't match expectations so we use antigradient path
    else:
        print("Gradient method used")
        d0=[-g0[0],-g0[1]]
        alfak=armijoRule.armijo(f0,x,d0,g0)
        x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]]
    k+=1
    print("point x",k,x)
