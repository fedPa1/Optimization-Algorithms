# @Author: Federico Pantini 

# email:fedpa35@gmail.com

# The script computes an alpha step that minimizes the given 2 variables function using the
# armijo's rule

import numpy 
import math
import random

def armijo(f,x,d,g): 
        f0=f
        x0=x
        g0=g
        d0=d
        if(numpy.dot(g0,d0)>=0):
            print("initial direction not descent, please choose a valid direction")
            return
        alpha=1 #alfa0 step >0
        delta=0.5 # 0<delta<1
        gamma=0.5 # 0<gamma<=0.5
        j=0
        while(j<=10**4):
            #fAlpha=  your f(x+alpha*d) function
            if (fAlpha<=f0+gamma*alpha*numpy.dot(g0,d0)):
                print("after",j,"iterations found","alpha=",alpha)
                return alpha
                break
            alpha*=delta
            j+=1
