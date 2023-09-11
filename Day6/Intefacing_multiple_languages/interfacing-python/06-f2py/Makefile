SHELL=/bin/sh

hello.so: hello.f90
	python -m numpy.f2py -c $< -m hello --quiet

clean:
	rm -f *~ *.so *.o
