// Catherine Galkina, group 524, year 2016
// File bsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <pthread.h>

#include "tools.h"
#include "point.h"

#define MIN_THREADS_NUM 4
#define MAX_HSORT_ELEMS 100000

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

