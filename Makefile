CXXFLAGS = -g -Wall
MPICXX = mpicxx

bisect: bisect.cpp
	$(MPICXX) $(CXXFLAGS) -O3 $^ -lpthread -o $@

bisect_seq: bisect_seq.cpp
	$(MPICXX) $(CXXFLAGS) -o3 $^ -o $@

all: bisect bisect_seq

clean:
	rm -f *.o bisect bisect_seq

