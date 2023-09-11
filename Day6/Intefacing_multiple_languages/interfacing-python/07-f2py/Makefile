SHELL=/bin/sh

sum.so: sum.f90
	f2py -c $< -m sum --quiet

clean:
	rm -f *~ *.so *.o
