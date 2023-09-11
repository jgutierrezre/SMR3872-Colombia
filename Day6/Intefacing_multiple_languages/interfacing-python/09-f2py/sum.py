#!/usr/bin/env python

from math import *
from sum import *

# fill data structure
num = 10
ilist = [i-1 for i in range(1,num)]

print ("Calling fortran with ints")
isum = sum_of_ints(ilist,num)
print ("Integer sum is: ", isum)

print ("Calling fortran double sum with integers")
dsum = sum_of_doubles(ilist,num)
print ("Double sum of ints is: ", dsum)

# fill data structure
dlist = [sqrt(float(i)) for i in range(1,num)]

print ("Calling fortran with doubles")
dsum = sum_of_doubles(dlist,num)
print ("Double sum is: ", dsum)
