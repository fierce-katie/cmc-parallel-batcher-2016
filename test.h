#ifndef BSORT_TEST
#define BSORT_TEST

#include <vector>
#include "tools.h"

bool is_sorted_binary(std::vector<int> &v, int n);
void to_binary_vector(int n, int total, std::vector<int> &res);
void generate_tests(int n, std::vector<std::vector<int> > &tests);
void run_tests(int n, std::vector<comparator> &cmp);

#endif
