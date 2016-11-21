// Catherine Galkina, group 524, year 2016
// File bsort.cpp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#include "tools.h"
#include "point.h"

using namespace std;

void join(vector<int> idx_up, int n0, vector<int> idx_down, int n1,
          vector<comparator> &cmp)
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
    vector<int> idx_up_even(n0_even);
    vector<int> idx_up_odd(n0_odd);

    int n1_even = n1/2;
    int n1_odd = n1 - n1_even;
    vector<int> idx_down_even(n1_even);
    vector<int> idx_down_odd(n1_odd);

    vector<int> idx_result;

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

void sort(vector<int>idx, int n, vector<comparator> &cmp)
{
    if (n == 1) {
        return;
    }

    int n0 = n/2;
    int n1 = n - n0;

    vector<int> idx_up;
    vector<int> idx_down;

    int i;
    for (i = 0; i < n0; i++)
        idx_up.push_back(idx[i]);
    for (i = n0; i < n; i++)
        idx_down.push_back(idx[i]);

    sort(idx_up, n0, cmp);
    sort(idx_down, n1, cmp);
    join(idx_up, n0, idx_down, n1, cmp);
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
    vector<comparator> cmp;
    vector<int> idx;
    for (int i = 0; i < procs; i++)
        idx.push_back(i);
    sort(idx, procs, cmp);

    // Calculating elems per processor
    int n = nx*ny;
    int fake = n % procs ? (procs - n % procs) : 0;
    int elems = n + fake;
    int proc_elems = elems / procs;
    // Initializing points
    srand(time(NULL));
    vector<Point> points = init_points(nx, ny, fake);
    if (!rank) {
        printf("Procs: %d\nElems per proc: %d\n", procs, proc_elems);
        print_points(points, elems);
    }

    // Printing result
    if (!rank) {
        for (int r = 0; r < procs; r++) {
            printf("\nI am number %d and my elems are:\n", r);
            for (int i = r; i < elems; i += procs)
                printf("%d: (%f, %f)\n", points[i].GetIndex(), points[i].GetX(),
                        points[i].GetY());
        }
    }

    MPI_Finalize();
    return 0;
}
#if 0
int main(int argc, char **argv)
{
    vector<comparator> cmp;
    int n0, n1;

    // Parsing command line arguments
    if (!check_args(argc, argv, n0, n1))
        return 1;

    // Sorting or joining
    printf("%d %d %d\n", n0, n1, 0);
    if (n1) {
        vector<int> idx_up;
        vector<int> idx_down;
        int i;
        for (i = 0; i < n0; i++)
            idx_up.push_back(i);
        for (i = n0; i < n0 + n1; i++)
            idx_down.push_back(i);
        join(idx_up, n0, idx_down, n1, cmp);
    } else {
        vector<int> idx;
        for (int i = 0; i < n0; i++)
            idx.push_back(i);
        sort(idx, n0, cmp);
    }

    // Printing result
    print_comparators(cmp);
    printf("%d\n", count_tacts(n0 + n1, cmp));

    // Run tests for sorting <= 15 elems
    if (!n1 && n0 <= 15)
        run_tests(n0, cmp);
    return 0;
}
#endif
