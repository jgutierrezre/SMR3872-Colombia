#!/usr/bin/env python

from ctypes import *
dso = CDLL("./sum.so")

isum = dso.sum_of_int(1,2)
print("Integer sum w/o prototypes is: ", isum)

#dsum = dso.sum_of_double(0.5, -2.5)
#print("Double sum w/o prototypes is: ", dsum)


dso.sum_of_int.argtypes = [c_int, c_int]
dso.sum_of_int.restype = c_int

dso.sum_of_double.argtypes = [c_double, c_double]
dso.sum_of_double.restype = c_double

isum = dso.sum_of_int(1,2)
print("Integer sum /w prototypes is: ", isum)

dsum = dso.sum_of_double(0.5, -2.5)
print("Double sum /w prototypes is: ", dsum)

