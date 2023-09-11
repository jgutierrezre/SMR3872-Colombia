#!/usr/bin/env python

from ctypes import *

class parm(Structure):
    _fields_ = [ ("type", c_int), ("label", c_char_p),
                 ("epsilon",c_double), ("sigma",c_double) ]

# import DSO
dso = CDLL("./data.so")

print ("Calling DSO")
print ("Global inum is: ", c_int.in_dll(dso,"inum").value)
print ("Global dnum is: ", c_double.in_dll(dso,"dnum").value)

p = parm(type=1, label=b"LJ-12-6", epsilon=0.1, sigma=3.4)
dso.pass_by_value(p)

p = parm(type=2, label=b"LJ-10-4", epsilon=0.5, sigma=4)
dso.pass_by_reference(byref(p))

