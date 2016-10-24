#ifndef BSORT_TOOLS
#define BSORT_TOOLS

#include <vector>

typedef std::pair<int, int> comparator;

void print_vector(std::vector<int> &v, int n);
void print_comparators(std::vector<comparator> &cmp);
bool check_args(int argc, char **argv, int &n0, int &n1);

#endif
