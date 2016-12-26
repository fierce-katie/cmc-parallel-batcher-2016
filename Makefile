CXXFLAGS = -g -Wall -O3
MPICXX = mpicxx

bsort: bsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lpthread -o $@

clean:
	rm -f *.o bsort

