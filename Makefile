OBJMODULES = tools.o point.o dhsort.o
CXXFLAGS = -g -Wall
MPICXX = mpic++

%.o: %.cpp %.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	$(MPICXX) $(CXXFLAGS) $^ -lpthread -o $@

clean:
	rm -f *.o bsort

