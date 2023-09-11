#!/usr/bin/env python

from ctypes import *

dso = CDLL("./sum.so")
num = 10
dlist = (c_double * num)() # constructor: (primitive * length)()
psum = 0.0
for i in range(num):
    dlist[i] = 1.0/3.0*(i*0.5)
    psum += dlist[i]

dso.sum_of_doubles.argtypes = [POINTER(c_double), c_int]
dso.sum_of_doubles.restype = c_double

csum = dso.sum_of_doubles(dlist, num)
print("Python sum: ", psum, " C sum: ", csum)
