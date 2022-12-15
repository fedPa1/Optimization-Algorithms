#@Author: Federico Pantini
#Contact: fedpa35@gmail.com

#The script is an algorithm that works in 2 variables quadratic functions to find the global minimum
#using a variation created by me of conjugated directions method, the algorithm seems to be working , however it's not broadly tested
#or demonstrated and it was purely done for fun, the main difference is that we can calculate a direction d1 given as initial direction the antigradient
#iteratively avoiding to initially set a Basis in R^n of Q-conjugated directions.

#References: https://www.stat.cmu.edu/~ryantibs/convexopt-F13/scribes/lec10.pdf

#Notation:
#^t indicates the transposed operator in arrays
#* indicates multiplication sign except in x* indicating a point in R^2
#sum|i=c1-->c2 of f(i) indicates summation from c1 to c2 in the input function f(i)

#Why use conjugated directions method?
#(1)Accelarate the convergence rate of steepest decent.
#(2)Avoiding the high computational cost Newton’s method
