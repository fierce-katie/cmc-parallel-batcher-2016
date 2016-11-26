// Catherine Galkina, group 524, year 2016
// File qsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "point.h"
#include "tools.h"

int main(int argc, char **argv)
{
    int nx, ny;

    // Parsing command line arguments
    if (!check_args(argc, argv, nx, ny))
        return 1;

    // Initializing points
    clock_t clocks = clock();
    Point *points = new Point[nx*ny];
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int idx = i*ny + j;
            points[idx] = Point(x(i, j), y(i, j), idx);
        }
    }
    qsort(points, nx*ny, sizeof(*points), compare_points);
    clocks = clock() - clocks;
    double total_time = (double)clocks/CLOCKS_PER_SEC;
    printf("Elems: %d\nTotal time: %f sec.\n", nx*ny, total_time);
    return 0;
}

