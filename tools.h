// Catherine Galkina, group 524, year 2016
// File tools.h
#ifndef BSORT_TOOLS
#define BSORT_TOOLS

#include <vector>
#include "point.h"

typedef std::pair<int, int> comparator;

void swap(comparator cmp, std::vector<int> &v);
void print_vector(std::vector<int> &v, int n);
void print_comparators(std::vector<comparator> &cmp);
void print_points(std::vector<Point> p, int n);
bool check_args(int argc, char **argv, int &n0, int &n1);
int count_tacts(int n, std::vector<comparator> &cmp);

#endif
