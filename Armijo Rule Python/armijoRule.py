# @Author: Federico Pantini 

# email:fedpa35@gmail.com

# The script computes an alpha step that minimizes the given 2 variables function using the
# armijo's rule

import numpy 
import math
import random

def armijo(f,x,d,g): 
        if(numpy.dot(g,d)>=0):
            print("initial direction not descent, please choose a valid direction")
            return
        alpha=1 #alfa0 step >0
        delta=0.5 # 0<delta<1
        gamma=0.5 # 0<gamma<=0.5
        j=0
        while(j<=10**4):
            fAlpha=  #your f(x+alpha*d) function
            if (fAlpha<=f+gamma*alpha*numpy.dot(g,d)):
                print("after",j,"iterations found","alpha=",alpha)
                return alpha
                break
            alpha*=delta
            j+=1
