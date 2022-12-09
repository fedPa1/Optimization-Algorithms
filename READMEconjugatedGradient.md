# Optimization-Algorithms
Implementation of linear optimization algorithms in python

#we can extend this algorithm to n dimensions but alfak step will need to be calculated with armijo rule
#and we can use  Fletcher-Reeves formula to calculate bk+1= norm(gk+1)/norm(gk)

#the conjugated gradient method differs from the conjugated directions in the choiche of d0 , and dk+1
# to calculate dk+1 we'll have to introduce a new variable bk+1=gk+1^t*Q*dk/dk^t*Q*dk , gk is the gradient
#calculate at xk, the next direction dk+1 is a perturbation of the antigradient in xk+1 obtained in the following way
# dk+1=-gk+1+bk+1*dk


#Notation:
#^t indicates the transposed operator in arrays
#* indicates multiplication sign except in x* indicating a point in R^2
#sum|i=c1-->c2 of f(i) indicates summation from c1 to c2 in the input function f(i)
