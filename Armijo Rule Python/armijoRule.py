# @Author: Federico Pantini
# Contact: fedpa35@gmail.com

# The script computes an alpha step that minimizes the given 2 variables function using the
# armijo's rule

import numpy 
import math
import random

x=[random.random()*1000,random.random()*1000] #initial random point

#we can input any 2 variable function,gradient,direction and point as long they match this format:
f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
d0=[-g0[0],-g0[1]] #d0 must be descent to guarantee that an alpha step will be found


class armijoClass:
    def __init__(self,f,x,d,g0):#takes as input a 2 variables function f(x[0],x[1]) ,x s.t. x=[x[0],x[1]]
                                # d s.t. d=[d[0],d[1]],g0 s.t. g0=[g0[0],g0[1]]
        self.f=f
        self.x=x
        self.d=d
        self.g0=g0
    def armijoRule(self): 
        f=self.f
        x=self.x
        d=self.d
        g0=self.g0
        if(numpy.dot(g0,d)>=numpy.dot(g0,x)):
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
armijoC=armijoClass(f0, x, d0, g0) 
alfa=armijoClass.armijoRule(armijoC)
