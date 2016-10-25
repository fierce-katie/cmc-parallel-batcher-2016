OBJMODULES = tools.o test.o
CXXFLAGS = -g -Wall

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

bsort: bsort.cpp $(OBJMODULES)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.o bsort
