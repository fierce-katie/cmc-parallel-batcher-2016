CXXFLAGS = -g -Wall -O3
MPICXX = mpicxx

bisect: bisect.cpp
	$(MPICXX) $(CXXFLAGS)  $^ -o $@

bisect_seq: bisect_seq.cpp
	$(MPICXX) $(CXXFLAGS)  $^ -o $@

all: bisect bisect_seq

clean:
	rm -f *.o bisect bisect_seq

