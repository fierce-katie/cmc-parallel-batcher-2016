// Catherine Galkina, group 524, year 2016
// File point.h
#ifndef BSORT_POINT
#define BSORT_POINT

#include <mpi.h>

struct Point {
    float coord[2];
    int index;
};

Point* init_points(int n, int ny, int procs, int proc_elems, int rank);
float x(int i, int j);
float y(int i, int j);
int compare_points(const void *a, const void *b);
MPI_Datatype pointType();

#endif
