# @Author: Federico Pantini 
# email:fedpa35@gmail.com

# The script computes the global minimum point of a 2 variable function (that admits a global minimum) using the quasi Newton Method

# more info about the method: https://en.wikipedia.org/wiki/Quasi-Newton_method

import numpy
import random

def addMatrix(M1,M2,m):
    if(m=="m"):
        for i in range(len(M1)):  
            for j in range(len(M1[0])):
                M1[i][j]+=M2[i][j]
        return M1
    else:
        for i in range(len(M1)):  
            for j in range(len(M1[0])):
                M1[i][j]+=M2
        return M1

def subtractVector(V1,V2):
    for i in range(len(V1)):
        V1[i]-=V2[i]
    return V1
#initial random point x in R^2 

x=[random.random()*1000,random.random()*1000]

#function (f0 function,g0 gradient,h0 hessian matrix), we can input any 2 variables function continously differentiable

f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1] 
g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]  #your input function,h0,g0,d0
d0=[-g0[0],-g0[1]]
B0=h0=[[10,-2],[-2,6]]

#iterative scheme xk+1=xk-gk/Bk , Bk corresponds to the approximate hessian matrix 
# we can choose to approximate B0 with different formulas , we choose Broyden's formula in this script

k=0

print("Initial point",x)

while k<=10**5:
    f0=5*x[0]**2-2*x[0]+5*x[1]+3*x[1]**2-2*x[0]*x[1]
    g0=[10*x[0]-2*x[1]-2,6*x[1]-2*x[0]+5]
    g0Old=g0.copy()
    xOld=x.copy()
    if (numpy.linalg.norm(g0)<=1.e-7):
        print("global minimum at",x,"after",k,"iterations")
        break

    x=[x[0]-numpy.dot(numpy.linalg.inv(B0),g0)[0],x[1]-numpy.dot(numpy.linalg.inv(B0),g0)[1]]    
    
    xDelta=[x[0]-xOld[0],x[1]-xOld[1]]
    gradDelta=[g0[0]-g0Old[0],g0[1]-g0Old[1]]

    sk=numpy.linalg.norm(xDelta)**2
    vk=xDelta
    wk=subtractVector(gradDelta,numpy.dot(B0,xDelta)) #Broyden's formula, defining sk,wk,vk for improving readability

    bDelta=(numpy.linalg.norm(vk)**2/sk) #Broyden formula
    B0=addMatrix(B0, bDelta,"s") #new B0 using Broyden's formula, note that a non-symmetric matrix can be produced even if the old one was
    
    
    #wk=subtractVector(gradDelta,numpy.dot(B0,xDelta)) #Davidon–Fletcher–Powell formula, defining sk,wk,vk for improving readability
    #vk=wk.copy()
    #sk=numpy.dot(vk,xDelta)

    #bDelta=(numpy.dot(wk,vk)/sk) #Davidon–Fletcher–Powell formula
    #B0=addMatrix(B0, bDelta,"s") #new B0 using Davidon–Fletcher–Powell formula, in this case will always be produced a symmetric matrix
    #print(bDelta)

    k+=1
    print("point x",k,x)
