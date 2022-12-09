# Optimization-Algorithms
Implementation of linear optimization algorithms in python


#Demonstration:

#let B a basis in R^n s.t. B={d0,d1,...,dn-1}, in this case we define a basis in R^2 B={d0,d1}
#we also know that exist n coefficients bi in R s.t. x*=b0*d0+b1*d1+...+bn-1*dn-1 so for each iteration k
#we can multiply by di^t*Q each member and obtain di^t*Q*x*=bi*di^t*Q*di-->bi=(di^t*Q*x*)/di^tQdi
#this because by hypothesis we are assuming di*Q*dj=0 for each i!=j

#we also know that x* in a quadratic form is x*=-Q^-1c so definitely we have:
#bi=-(di^t*c )/di^tQdi

#So we can write in the general form x*=x0+sum|i=0-->n-1 of -di^t*c*di/di^t*Q*di , x0 in R^2 initial starting point
#since our problem x*=alfa0d0+alfa1d1+...+alfan-1dn-1
#we can define (after eliminating i!=j components premultiplicating for di^t*Q in the k iteration) alfak=di^t*Q*x*/di^t*Q*di
#since x*=-Q^-1*c--> alfak=dk^t*c/dk^t*Q*dk with alfak being the optimal step from the point xk towards the direction dk

#We can now redefine our problem x*=sum|i=0-->d-1 of dk^t*c*dk/dk^t*Q*dk and xk+1=xk+alfak*dk

#dk+1 will be calculated s.t. dk^t*Q*dk+1=0 taken d0 as initial direction


#Notation:
#^t indicates the transposed operator in arrays
#* indicates multiplication sign except in x* indicating a point in R^2
#sum|i=c1-->c2 of f(i) indicates summation from c1 to c2 in the input function f(i)


