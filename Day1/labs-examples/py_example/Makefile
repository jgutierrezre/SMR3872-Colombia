TIMEIT = python3 -m timeit -s 'import sum_abs_timing as t'

FUNCTIONS = \
ranged_forloop \
functional_numpy

all: $(addsuffix .time, $(FUNCTIONS))

%.time: 
	$(TIMEIT) "t.$*(t.xlst)"
	$(TIMEIT) "t.$*(t.xarr)"
