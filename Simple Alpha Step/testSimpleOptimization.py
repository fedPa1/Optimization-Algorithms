#@Author: Federico Pantini
#Contact: fedpa35@gmail.com
#
#The script computes a random alpha step given a descent direction in a 2 variable function
#

import numpy
import math
import random

def simpleAlphaStep(f,g,x,d):
    alfa=random.random()*100 #alfa0 step >0
    #fAlfa=#your f(x+alpha*d) function
    j=0
    if(numpy.dot(g,d)>=0):
            print("initial direction not descent, please choose a valid direction")
            return
    while(j<=10**5):
        fAlfa=#your f(x+alpha*d) function
        if (fAlfa<=f):
            print("after",j,"iterations found",alfa)
            return alfa
        alfa*=0.5
        j+=1
    return

x=[random.random()*1000,random.random()*1000]
f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1]  #your function , x point etc..
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
d0=[-g0[0],-g0[1]]
simpleAlphaStep(f0,g0,x,d0)
