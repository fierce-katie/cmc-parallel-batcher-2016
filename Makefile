OBJMODULES = tools.o point.o
CXXFLAGS = -g -Wall
MPICXX = mpicxx
MPIXLCXX = mpixlcxx_r

%.o: %.cpp %.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.o bsort qsort

