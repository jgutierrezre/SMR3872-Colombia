#!/usr/bin/env python

from ctypes import *

hello = CDLL("./hello.so")

hello.hello.argtypes = [c_char_p]

# passing string as "byte sequence" to C
hello.hello(b"World")
hello.hello(b"Axel")
# regular strings must be encoded
hello.hello(str(666).encode('UTF-8'))

makeupper = CDLL("./makeupper.so")
makeupper.makeupper.argtypes = [c_char_p]
makeupper.makeupper.restype = c_char_p

# create string buffer to avoid segfault
buf = create_string_buffer(b"Axel")
# result must be decoded from a byte sequence into a string
upper = makeupper.makeupper(buf).decode('UTF-8')
print("Hello, %s!" % upper)

