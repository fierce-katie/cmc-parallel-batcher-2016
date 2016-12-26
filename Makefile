CXXFLAGS = -g -Wall -O3
MPICXX = mpicxx

bsort: bsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lpthread -o $@

deompose_seq:
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.o bsort decompose_seq

