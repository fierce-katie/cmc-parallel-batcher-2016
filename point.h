// Catherine Galkina, group 524, year 2016
// File point.h
#ifndef BSORT_POINT
#define BSORT_POINT

#include <vector>

class Point {
    float coord[2];
    int index;
public:
    Point() { coord[0] = coord[1] = index = 0; }
    Point(float x, float y, int idx) : index(idx)
        { coord[0] = x; coord[1] = y; }
    float GetX() const { return coord[0]; }
    float GetY() const { return coord[1]; }
    int GetIndex() const { return index; }
};

std::vector<Point> init_points(int n1, int n2, int fake);
float x(int i, int j);
float y(int i, int j);

#endif
