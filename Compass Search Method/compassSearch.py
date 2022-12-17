#@Author: Federico Pantini
#contact: fedpa35@gmail.com

# The script computes the global minimum point of a 2 variable function (that admits a global minimum) using the compass search Method

import numpy 
import random

x=[random.random()*1000,random.random()*1000]

f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]  #your input function,g0,d0
B=[[1,0],[-1,0],[0,1],[0,-1]] #R^2 basis composed of +-e1,+-e2 directions that will be used
                              #to compare the f(xk+1) values and choose the one that minimise the
                              #function the most

delta=1                     #delta step to take in Bi direction

k=0
while(k<=10**5 and delta>=1.e-20):
    f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    if(numpy.linalg.norm(g0)<=1.e-7):
        print("global minimum at",x,"after",k,"iterations")
        break
    options=[]
    for d in B:
        fValue=5*(x[0]+delta*d[0])**2-2*(x[0]+delta*d[0])+5*(x[1]+delta*d[1])+3*(x[1]+delta*d[1])**2-2*(x[0]+delta*d[0])*(x[1]+delta*d[1])
        options.append(fValue)
    best=f0
    index=0
    for i in range(len(options)):
        if(options[i]<best):
            best=options[i]
            index=i
    if(best<f0):
        x=[x[0]+delta*B[index][0],x[1]+delta*B[index][1]]
        delta=1
    else:
        delta/=2
    k+=1
