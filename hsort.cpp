// Catherine Galkina, group 524, year 2016
// File hsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct Point {
    float coord[2];
    int index;
};

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
  Point *aptr = (Point * const)a;
  Point *bptr = (Point * const)b;
  float ax = aptr->coord[0];
  float bx = bptr->coord[0];

  if (ax == bx)
      return 0;
  else if (ax > bx)
      return 1;
  return -1;
}

void heapify(Point* a, int i, int n)
{
    int imax, l, r;
    Point tmp;
    while (i < n) {
        imax = i;
        l = 2*i + 1;
        r = l + 1;
        if (l < n && (compare_points(&a[l], &a[imax]) > 0))
            imax = l;
        if (r < n && (compare_points(&a[r], &a[imax]) > 0))
            imax = r;
        if (imax == i)
            return;
        tmp = a[i];
        a[i] = a[imax];
        a[imax] = tmp;
        i = imax;
    }
}

void make_heap(Point* a, int n)
{
    for (int i = n/2 - 1; i >= 0; i--)
        heapify(a, i, n);
}

void hsort(Point *a, int n)
{
    make_heap(a, n);
    Point tmp;
    for (int i = n - 1; i > 0; i--) {
        tmp = a[0];
        a[0] = a[i];
        a[i] = tmp;
        heapify(a, 0 ,i);
    }
}

bool check_args(int argc, char **argv, int &nx, int &ny)
{
    if (argc < 3) {
        printf("Wrong arguments. Usage: bsort nx ny\n");
        return false;
    }
    int check = sscanf(argv[1], "%d", &nx);
    if (!check) {
        printf("nx must be int: %s\n", argv[1]);
        return false;
    }
    check = sscanf(argv[2], "%d", &ny);
    if (!check) {
        printf("ny must be int: %s\n", argv[2]);
        return false;
    }
    if (!((nx >= 1) && (ny >= 1))) {
        printf("Wrong n1 or n2\n");
        return false;
    }
    return true;
}

int main(int argc, char **argv)
{
    int nx, ny;

    // Parsing command line arguments
    if (!check_args(argc, argv, nx, ny))
        return 1;
    int n = nx * ny;

    // Initializing points
    srand(time(NULL));
    Point *points = new Point[n];
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Point p;
            p.coord[0] = x(i, j);
            p.coord[1] = y(i, j);
            int idx = i*ny + j;
            p.index = idx;
            points[idx] = p;
        }
    }

    // Sorting
    clock_t sort_time = clock();
    hsort(points, n);
    sort_time = clock() - sort_time;

#if 0
    for (int i = 0; i < n; i++)
        printf("%f\n", points[i].coord[0]);
#endif

    printf("Sort time (hsort, %d x %d): %lf\n", nx , ny,
            (double)sort_time/CLOCKS_PER_SEC);

    delete [] points;
    return 0;
}
