OBJMODULES = tools.o test.o point.o
CXXFLAGS = -g -Wall

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	mpic++ $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.o bsort
