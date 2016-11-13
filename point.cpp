// Catherine Galkina, group 524, year 2016
// File point.cpp

#include "point.h"

using namespace std;

Point* init_points(int n1, int n2)
{
    Point *res = new Point[n1 * n2];
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            int idx = i*n2 + j;
            res[idx] = Point(0, 0, idx);
        }
    }
    return res;
}
