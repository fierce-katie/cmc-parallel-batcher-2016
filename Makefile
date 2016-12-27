CXXFLAGS = -g -Wall
MPICXX = mpicxx

bisect: bisect.cpp
	$(MPICXX) $(CXXFLAGS) -O3 $^ -lpthread -o $@

bisect_seq: bisect_seq.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

all: bisect bisect_seq

clean:
	rm -f *.o bisect bisect_seq

