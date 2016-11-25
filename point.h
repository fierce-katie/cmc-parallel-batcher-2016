// Catherine Galkina, group 524, year 2016
// File point.h
#ifndef BSORT_POINT
#define BSORT_POINT

#include <mpi.h>

class Point {
    float coord[2];
    int index;
public:
    Point() { coord[0] = coord[1] = index = -1; }
    Point(float x, float y, int idx) : index(idx)
        { coord[0] = x; coord[1] = y; }
    float GetX() const { return coord[0]; }
    float GetY() const { return coord[1]; }
    int GetIndex() const { return index; }
    MPI_Datatype getType();
    bool operator<(const Point &other) const
        { return coord[0] < other.coord[0]; }
    bool operator>(const Point &other) const
        { return coord[0] > other.coord[0]; }
    bool operator<=(const Point &other) const
        { return coord[0] <= other.coord[0]; }
    bool operator>=(const Point &other) const
        { return coord[0] >= other.coord[0]; }
};

Point* init_points(int n, int ny, int procs, int proc_elems, int rank);
float x(int i, int j);
float y(int i, int j);
int compare_points(const void *a, const void *b);

#endif
