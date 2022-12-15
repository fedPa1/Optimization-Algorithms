#@Author: Federico Pantini
#Contact: fedpa35@gmail.com

#The script is an algorithm that works in 2 variables quadratic functions to find the global minimum using the conjugated gradient method

#References: https://www.stat.cmu.edu/~ryantibs/convexopt-F13/scribes/lec10.pdf

#the conjugated gradient method differs from the conjugated directions in the choiche of d0 , and dk+1
#to calculate dk+1 we'll have to introduce a new variable bk+1=gk+1^tQdk/dk^t*Q*dk , the next direction dk+1
#is a perturbation of the antigradient in xk+1 obtained in the following way dk+1=-gk+1+bk+1*dk

