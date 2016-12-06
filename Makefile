CXXFLAGS = -g -Wall -O3
MPICXX = mpicxx
SEQ = qsort dsort hsort dhsort

bsort: bsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lpthread -o $@

clean:
	rm -f *.o bsort $(SEQ)

qsort: qsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -o $@

hsort: hsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -o $@

dhsort: dhsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -o $@

dsort: dsort.cpp
	$(MPICXX) $(CXXFLAGS) $^ -o $@

seq: $(SEQ)

all: bsort seq
