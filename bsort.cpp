// Catherine Galkina, group 524, year 2016
// File bsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cstddef>


#define MIN_THREADS_NUM 4
#define MAX_HSORT_ELEMS 100000

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
        res[k] = p;
    }
    return res;
}

MPI_Datatype pointType()
{
    MPI_Datatype point;
    MPI_Datatype types[2] = { MPI_FLOAT, MPI_INT };
    int blocks[2] = { 2, 1 };
    MPI_Aint disps[2] = { offsetof(Point, coord), offsetof(Point, index) };
    MPI_Type_create_struct(2, blocks, disps, types, &point);
    MPI_Type_commit(&point);
    return point;
}

typedef std::pair<int, int> comparator;

void swap(comparator cmp, std::vector<int> &v)
{
    int fst = cmp.first;
    int snd = cmp.second;
    if (v[fst] > v[snd]) {
        int tmp = v[fst];
        v[fst] = v[snd];
        v[snd] = tmp;
    }
}

void print_vector(std::vector<int> &v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

void print_comparators(std::vector<comparator> &cmp)
{
    std::vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++)
        printf("%d %d\n", it->first, it->second);
    printf("%lu\n", cmp.size());
}

void print_points(Point* p, int n, int rank, const char *comment)
{
    for (int i = 0; i < n; i++) {
        Point point = p[i];
        if (point.index >= 0)
            printf("%d %d: %f %s\n", rank, point.index, point.coord[0],
                   comment);
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

int count_tacts(int n, std::vector<comparator> &cmp)
{
    std::vector<int> v(n);
    std::vector<comparator>::iterator it;
    int max;
    for (it = cmp.begin(); it != cmp.end(); it++) {
        int fst = it->first;
        int snd = it->second;
        max = v[fst] > v[snd] ? v[fst] : v[snd];
        v[fst] = max + 1;
        v[snd] = max + 1;
    }
    max = 0;
    for (int i = 0; i < n; i++)
        if (v[i] > max)
            max = v[i];
    return max;
}

void swap_ptr(void *ptr1_ptr, void *ptr2_ptr)
{
    void **ptr1 = (void **)ptr1_ptr;
    void **ptr2 = (void **)ptr2_ptr;

    void *tmp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = tmp;
}

void write_output(Point *a, int n, char *file, int nx, int ny, int rank)
{
    float *a_out = new float[n];
    int out_cnt = 0;
    for (int i = 0; i < n; i++) {
        if (a[i].index != -1)
            a_out[out_cnt++] = a[i].coord[0];
    }

    MPI_Status status;
    MPI_File output;
    int res = MPI_File_open(MPI_COMM_WORLD, file,
                            MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                            &output);
    if (res != MPI_SUCCESS) {
        if (rank == 0) {
            printf("Cannot open %s\n", file);
        }
        MPI_Finalize();
        exit(1);
    }
    MPI_File_set_size(output, 0);
    MPI_File_write_ordered(output, a_out, out_cnt, MPI_FLOAT, &status);
    MPI_File_close(&output);
    delete [] a_out;

    return;
}
inline int compare_points(const void *a, const void *b)
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

struct PthreadArgs {
    pthread_t tid;
    Point *a;
    int n;
};

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

void* hsort_start(void *arg)
{
    PthreadArgs *args = (PthreadArgs*)arg;
    hsort(args->a, args->n);
    return NULL;
}


void hsort_threads(Point* a, int n, int nthreads)
{
    PthreadArgs *pthread_args = new PthreadArgs[nthreads];
    int tmp = ceil(n / (double)nthreads);
    int elems;
    int offset = 0;
    for (int th = 0; th < nthreads; th++) {
        if (n - offset >= tmp)
            elems = tmp;
        else
            elems = (n - offset > 0) ? n - offset : 0;

        pthread_args[th].a = a + offset;
        pthread_args[th].n = elems;
        pthread_create(&pthread_args[th].tid, NULL, hsort_start,
                       &pthread_args[th]);
        offset += elems;
    }
    for (int th = 0; th < nthreads; th++)
        pthread_join(pthread_args[th].tid, NULL);
    delete [] pthread_args;
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

void dhsort(Point *a, int n)
{
    int nthreads = ceil(n / (double)MAX_HSORT_ELEMS);
    nthreads = nthreads > MIN_THREADS_NUM ? nthreads : MIN_THREADS_NUM;
    hsort_threads(a, n, nthreads);
    dsort(a, n, ceil(n / (double)nthreads));

    return;
}

void join(std::vector<int> idx_up, int n0, std::vector<int> idx_down, int n1,
          std::vector<comparator> &cmp)
{
    int n = n0 + n1;
    if (n == 1)
        return;

    if (n0 == 1 && n1 == 1) {
        cmp.push_back(comparator(idx_up[0], idx_down[0]));
        return;
    }

    int n0_even = n0/2;
    int n0_odd = n0 - n0_even;
    std::vector<int> idx_up_even(n0_even);
    std::vector<int> idx_up_odd(n0_odd);

    int n1_even = n1/2;
    int n1_odd = n1 - n1_even;
    std::vector<int> idx_down_even(n1_even);
    std::vector<int> idx_down_odd(n1_odd);

    std::vector<int> idx_result;

    int i, i0 = 0, i1 = 0;
    for (i = 0; i < n0; i++)
        if (i%2) {
            idx_up_even[i0] = idx_up[i];
            i0++;
        } else {
            idx_up_odd[i1] = idx_up[i];
            i1++;
        }
    i0 = i1 = 0;
    for (i = 0; i < n1; i++)
        if (i%2) {
            idx_down_even[i0] = idx_down[i];
            i0++;
        } else {
            idx_down_odd[i1] = idx_down[i];
            i1++;
        }

    join(idx_up_odd, n0_odd, idx_down_odd, n1_odd, cmp);
    join(idx_up_even, n0_even, idx_down_even, n1_even, cmp);

    for (i = 0; i < n0; i++)
        idx_result.push_back(idx_up[i]);
    for (i = 0; i < n1; i++)
        idx_result.push_back(idx_down[i]);

    for (int i = 1; i < n - 1; i += 2)
        cmp.push_back(comparator(idx_result[i], idx_result[i + 1]));
}

void sort(std::vector<int>idx, int n, std::vector<comparator> &cmp)
{
    if (n == 1) {
        return;
    }

    int n0 = n/2;
    int n1 = n - n0;

    std::vector<int> idx_up;
    std::vector<int> idx_down;

    int i;
    for (i = 0; i < n0; i++)
        idx_up.push_back(idx[i]);
    for (i = n0; i < n; i++)
        idx_down.push_back(idx[i]);

    sort(idx_up, n0, cmp);
    sort(idx_down, n1, cmp);
    join(idx_up, n0, idx_down, n1, cmp);
}

void make_comparators(int procs, std::vector<comparator> &cmp)
{
    std::vector<int> idx;
    for (int i = 0; i < procs; i++)
        idx.push_back(i);
    sort(idx, procs, cmp);
    return;
}

int main(int argc, char **argv)
{
    int nx, ny;

    // Parsing command line arguments
    if (!check_args(argc, argv, nx, ny))
        return 1;

    MPI_Init(&argc, &argv);

    int rank, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    // Comparators network
    std::vector<comparator> cmp;
    make_comparators(procs, cmp);

    // Calculating elems per processor
    int n = nx*ny;
    int fake = n % procs ? (procs - n % procs) : 0;
    int elems = n + fake;
    int proc_elems = elems / procs;

    // Initializing points
    srand(time(NULL) + rank);
    Point *proc_points =
        init_points(n, ny, procs, proc_elems, rank);

    // Sorting
    double start_time = MPI_Wtime();
    dhsort(proc_points, proc_elems);

    // Exchanging elements
    Point *tmp_points = new Point[proc_elems];
    Point *other_points = new Point[proc_elems];
    MPI_Status status;
    MPI_Datatype MPI_POINT = pointType();
    std::vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++) {
        if (rank == it->first) {
            MPI_Send(proc_points, proc_elems, MPI_POINT,
                     it->second, 0, MPI_COMM_WORLD);
            MPI_Recv(other_points, proc_elems, MPI_POINT,
                     it->second, 0, MPI_COMM_WORLD, &status);
            int idx = 0;
            int other_idx = 0;
            for (int tmp_idx = 0; tmp_idx < proc_elems; tmp_idx++) {
                Point my = proc_points[idx];
                Point other = other_points[other_idx];
                if (my.coord[0] < other.coord[0]) {
                    tmp_points[tmp_idx] = my;
                    idx++;
                } else {
                    tmp_points[tmp_idx] = other;
                    other_idx++;
                }
            }
            swap_ptr(&proc_points, &tmp_points);
        }

        if (rank == it->second) {
            MPI_Recv(other_points, proc_elems, MPI_POINT,
                     it->first, 0, MPI_COMM_WORLD, &status);
            MPI_Send(proc_points, proc_elems, MPI_POINT,
                     it->first, 0, MPI_COMM_WORLD);
            int idx = proc_elems - 1;
            int other_idx = proc_elems - 1;
            for (int tmp_idx = proc_elems - 1; tmp_idx >= 0; tmp_idx--) {
                Point my = proc_points[idx];
                Point other = other_points[other_idx];
                if (my.coord[0] > other.coord[0]) {
                    tmp_points[tmp_idx] = my;
                    idx--;
                } else {
                    tmp_points[tmp_idx] = other;
                    other_idx--;
                }
            }
            swap_ptr(&proc_points, &tmp_points);
        }
    }
    double end_time = MPI_Wtime();
    double time = end_time - start_time;
    double sort_time = 0;
    MPI_Reduce(&time, &sort_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!rank) {
        printf("Elems: %d\nProcs: %d\n", n, procs);
        printf("Sort time: %f sec.\n", sort_time);
    }

    if (argc > 3) {
        write_output(proc_points, proc_elems, argv[3], nx, ny, rank);
    }

    delete [] proc_points;
    MPI_Finalize();
    return 0;
}

