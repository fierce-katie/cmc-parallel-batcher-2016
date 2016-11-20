// Catherine Galkina, group 524, year 2016
// File point.cpp

#include <stdlib.h>

#include "point.h"

using namespace std;

vector<Point> init_points(int n1, int n2, int fake)
{
    int n = n1 * n2 + fake;
    vector<Point> res(n);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            int idx = i*n2 + j;
            res[idx] = Point(x(i, j), y(i, j), idx);
        }
    }
    for (int i = 0; i < fake; i++)
        res[n1*n2 + i] = Point(0, 0, -i-1);
    return res;
}

float x(int i, int j)
{
    return (float)rand()/(float)(RAND_MAX/(i*j+1));
}

float y(int i, int j)
{
    return (float)rand()/(float)(RAND_MAX/(i*j+1));
}
