OBJMODULES = tools.o point.o
CXXFLAGS = -g -Wall
MPICXX = mpicxx
MPIXLCXX = mpixlcxx_r

%.o: %.cpp %.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	$(MPICXX) $^ -o $@

clean:
	rm -f *.o bsort qsort hsort dhsort

qsort: qsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

hsort: hsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

dhsort: dhsort.cpp
	$(CXX) $(CXXFLAGS) $^ -lpthread -o $@
