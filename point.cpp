// Catherine Galkina, group 524, year 2016
// File point.cpp

#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include "point.h"

Point* init_points(int nx, int ny, int rank, int n, bool fake)
{
    srand(time(NULL) + rank);
    Point *res = new Point[n];
    //int i0, i1, j0, j1; //FIXME
    int k = 0;
    for (int i = 0, j = 0; i < n; i++)
        res[k++] = Point(x(i, j), y(i, j), i*ny + j);
    if (fake)
        res[n] = Point();
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

int compare_points(const void *a, const void *b)
{
  float ax = ((Point * const)a)->GetX();
  float bx = ((Point * const)b)->GetX();
  if (ax < bx)
      return -1;
  else
      return 1;
}
