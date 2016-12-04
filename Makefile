OBJMODULES = tools.o point.o
CXXFLAGS = -g -Wall -O3 -fomit-frame-pointer -fno-zero-initialized-in-bss -funsafe-loop-optimizations -ffast-math -funroll-all-loops
MPICXX = mpicxx
MPIXLCXX = mpixlcxx_r
SEQ = qsort dsort hsort dhsort

%.o: %.cpp %.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
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
