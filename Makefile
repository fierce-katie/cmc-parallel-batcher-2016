CXXFLAGS = -g -Wall
MPICXX = mpicxx

bsort: bsort.cpp
	$(MPICXX) $(CXXFLAGS) -O3 $^ -lpthread -o $@

deompose_seq:
	$(CXX) $(CXXFLAGS) $^ -o $@

all: bsort decompose_seq

clean:
	rm -f *.o bsort decompose_seq

