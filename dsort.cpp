// Catherine Galkina, group 524, year 2016
// File dsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

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

void dsort(Point *array, int n)
{
    Point *a = array;
    Point *b = new Point[n];
    Point *c;

    for (int i = 1; i < n ; i *= 2) {
        for (int j = 0; j < n; j = j + 2*i) {
            int r = j + i;

            int n1 = (i < n - j) ? i : n - j;
            int n2 = (i < n - r) ? i : n - r;
            n1 = (n1 < 0) ? 0 : n1;
            n2 = (n2 < 0) ? 0 : n2;

            for (int ia = 0, ib = 0, k = 0; k < n1 + n2; k++) {
                if (ia >= n1)
                    b[j+k] = a[r+ib++];
                else if (ib >= n2)
                    b[j+k]=a[j+ia++];
                else if (compare_points(&a[j+ia], &a[r+ib]) < 0)
                    b[j+k]=a[j+ia++];
                else
                    b[j+k]=a[r+ib++];
            }
        }
        c = a;
        a = b;
        b = c;
    }

    c = a;
    a = b;
    b = c;
    if (b != array) {
        memcpy(array, b, n*sizeof(*array));
        delete [] b;
    } else {
        delete [] a;
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
    dsort(points, n);
    sort_time = clock() - sort_time;

#if 0
    for (int i = 0; i < n; i++)
        printf("%f\n", points[i].coord[0]);
#endif

    printf("Sort time (dsort, %d x %d): %lf\n", nx, ny,
           (double)sort_time/CLOCKS_PER_SEC);

    delete [] points;
    return 0;
}
