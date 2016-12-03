OBJMODULES = tools.o point.o
CXXFLAGS = -g -Wall -O3
MPICXX = mpicxx
MPIXLCXX = mpixlcxx_r
SEQ = qsort dsort hsort dhsort

%.o: %.cpp %.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	$(MPICXX) $^ -lpthread  -o $@

clean:
	rm -f *.o bsort $(SEQ)

qsort: qsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

hsort: hsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

dhsort: dhsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

dsort: dsort.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

seq: $(SEQ)

all: bsort seq
