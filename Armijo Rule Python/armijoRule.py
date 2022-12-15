# @Author: Federico Pantini 

# email:fedpa35@gmail.com

# The script computes an alpha step that minimizes the given 2 variables function using the
# armijo's rule

import numpy 
import math
import random

x=[random.random()*100,random.random()*100] #initial random point

#we can input any 2 variable function,gradient,direction and point as long they match this format:
f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
#d0=[-g0[0],-g0[1]] #d0 must be descent to guarantee that an alpha step will be found
d0=[-g0[0],-g0[1]]


def armijo(f,x,d,g0): 
        print(x,d,numpy.dot(g0,d),numpy.linalg.norm(g0)**2)
        print(1.e-3*numpy.linalg.norm(g0))
        if(numpy.dot(g0,d)>=0):
            print("initial direction not descent, please choose a valid direction")
            return
        alfa=1 #alfa0 step >0
        delta=0.5 # 0<delta<1
        gamma=0.5 # 0<gamma<=0.5
        j=0
        while(j<=10**5):
            fAlfa=5*(x[0]+alfa*d[0])**2-2*(x[0]+alfa*d[0])+5*(x[1]+alfa*d[1])+3*(x[1]+alfa*d[1])**2-2*(x[0]+alfa*d[0])*(x[1]+alfa*d[1])
            if (fAlfa<=f+gamma*alfa*numpy.dot(g0,d)):
                print("after",j,"iterations found","alfa=",alfa)
                return alfa
                break
            alfa=alfa*delta
            j+=1
alfa=armijo(f0,x,d0,g0)
