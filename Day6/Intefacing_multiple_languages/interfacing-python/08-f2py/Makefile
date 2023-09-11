SHELL=/bin/sh

hello.so: hello.f90
	f2py -c $< -m hello --quiet

clean:
	rm -f *~ *.so *.o
