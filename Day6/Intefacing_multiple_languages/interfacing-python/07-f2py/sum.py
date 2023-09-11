#!/usr/bin/env python

from sum import *

print("Calling fortran with ints")
isum = sum_of_int(1,2)
print("Integer sum is: ", isum)

dsum = sum_of_double(1,2)
print("Double sum of ints is: ", dsum)

dsum = sum_of_double(0.5,2.5)
print("Double sum of doubles is: ", dsum)

print("Done")
