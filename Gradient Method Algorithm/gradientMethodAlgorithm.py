#@Author: Federico Pantini
#contact: fedpa35@gmail.com

import numpy
import math
import random
import armijoRule

#initial random point x in R^2 

x=[random.random()*1000,random.random()*1000]

#function (f0 function,g0 gradient)
#we can input any 2 variables function

f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]

k=0
while (k<=10**4):
    f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1]
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    print("k=",k,"Gradient",numpy.linalg.norm(g0),"point",x)
    if (numpy.linalg.norm(g0)<=1.e-6): #algorithm found the only critical point and must be global in the strictly convex function
                                        #Generally we must put a certain tolerance to the gradient tending towards 0
        print("global minimum at",x,"after",k,"iterations")
        break
    d0=[-g0[0],-g0[1]] # we always use the antigradient as direction dk
    alfak=armijoRule.armijo(f0, x, d0, g0) #we use armijo to calculate the alphak step
    x=[x[0]+alfak*d0[0],x[1]+alfak*d0[1]] #new xk+1 point
    k+=1

