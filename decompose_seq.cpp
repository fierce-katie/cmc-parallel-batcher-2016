// Catherine Galkina, group 524, year 2016
// File decompose_seq.cpp

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

bool axis = true;

struct Point {
    float coord[2];
    int index;
    int domain;
};

float x(int i, int j)
{
    return (float)rand()/(float)(RAND_MAX/(i*j+1));
}

float y(int i, int j)
{
    return (float)rand()/(float)(RAND_MAX/(i*j+1));
}

bool check_args(int argc, char **argv, int &nx, int &ny, int &k)
{
    if (argc < 4) {
        printf("Wrong arguments. Usage: decompose nx ny k\n");
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
    check = sscanf(argv[3], "%d", &k);
    if (!check) {
        printf("k must be int: %s\n", argv[3]);
        return false;
    }
    if (!((nx >= 1) && (ny >= 1) && (k >= 1))) {
        printf("Wrong nx or nx or k\n");
        return false;
    }
    return true;
}

Point* init_points(int n, int ny, int procs, int proc_elems, int rank)
{
    Point *res = new Point[proc_elems];
    int tmp = n/procs;
    int not_fake = n % procs;
    int real_elems = rank < not_fake ? tmp + 1 : tmp;
    int delta;
    if (rank < not_fake)
        delta = rank*proc_elems;
    else
        delta = not_fake*proc_elems + tmp*(rank - not_fake);
    for (int k = 0; k < real_elems; k++) {
        int i = (k + delta) / ny, j = (k + delta) % ny;
        Point p;
        p.coord[0] = x(i, j);
        p.coord[1] = y(i, j);
        p.index = i*ny + j;
        p.domain = -1;
        res[k] = p;
    }
    return res;
}

inline int compare_points(const void *a, const void *b)
{
  Point *aptr = (Point * const)a;
  Point *bptr = (Point * const)b;
  int c = axis ? 0 : 1;
  float ax = aptr->coord[c];
  float bx = bptr->coord[c];

  if (ax == bx)
      return 0;
  else if (ax > bx)
      return 1;
  return -1;
}

void dsort (Point *array, int n, int sorted)
{
    Point *a = array;
    Point *b = new Point[n];
    Point *c;

    for (int i = sorted; i < n ; i *= 2) {
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

void print_domains(Point *points, int nx, int ny)
{
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      Point cur = points[i*ny + j];
      printf("%d %d %f %f %d\n", i, j, cur.coord[0], cur.coord[1], cur.domain);
    }
  }
}

int main(int argc, char **argv)
{
  int nx, ny, k;

  // Parsing command line arguments
  if (!check_args(argc, argv, nx, ny, k))
    return 1;

  int n = nx*ny;
  srand(time(NULL));
  Point *points = init_points(n, ny, 1, n, 0);

  clock_t decompose_time = clock();

  decompose_time = clock() - decompose_time;
  printf("Decomposition time: %lf\n", (double)decompose_time/CLOCKS_PER_SEC);
  print_domains(points, nx, ny);

}
