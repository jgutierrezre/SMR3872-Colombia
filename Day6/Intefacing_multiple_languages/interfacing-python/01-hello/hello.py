#!/usr/bin/env python

from ctypes import CDLL

# import shared object and assign to handle
dso = CDLL("./hello.so")

# call symbol in shared object. Function name becomes argument to dlsym()

dso.hello()

