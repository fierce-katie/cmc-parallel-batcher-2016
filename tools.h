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
void print_points(Point* p, int n, int rank, const char *comment);
bool check_args(int argc, char **argv, int &n0, int &n1);
int count_tacts(int n, std::vector<comparator> &cmp);
void swap_ptr(void *ptr1_ptr, void *ptr2_ptr);
void write_output(Point *a, int n, char *file, int nx, int ny, int rank);

#endif
