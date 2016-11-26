// Catherine Galkina, group 524, year 2016
// File dhsort.cpp

#include "dhsort.h"
#include "tools.h"

#include <pthread.h>

#define MIN_THREADS_NUM 4
#define MAX_HSORT_ELEMS 10000

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
        if (l < n && a[l] > a[imax])
            imax = l;
        if (r < n && a[r] > a[imax])
            imax = r;
        if (imax == i)
            return;
        tmp = a[i];
        a[i] = a[imax];
        a[imax] = tmp;
        i = imax;
    }
    return;
}

void make_heap(Point* a, int n)
{
    for (int i = n/2 - 1; i >= 0; i--)
        heapify(a, i, n);
    return;
}

void *hsort(void *arg)
{
    PthreadArgs *points = static_cast<PthreadArgs *>(arg);
    Point *a = points->a;
    int n = points->n;
    make_heap(a, n);
    Point tmp;
    for (int i = n - 1; i > 0; i--) {
        tmp = a[0];
        a[0] = a[i];
        a[i] = tmp;
        heapify(a, 0 ,i);
    }
    return NULL;
}

void hsort_threads(Point* a, int n, int nthreads)
{
    PthreadArgs *pthread_args = new PthreadArgs[nthreads];
    int tmp = n/nthreads;
    int full = n % nthreads;
    for (int th = 0; th < nthreads; th++) {
        int elems = th < full ? tmp + 1 : tmp;
        int offset;
        if (th < full)
            offset = th*(tmp + 1);
        else
            offset = full*(tmp + 1) + (th - full)*tmp;
        pthread_args[th].a = a + offset;
        pthread_args[th].n = elems;
        pthread_create(&pthread_args[th].tid, NULL, hsort, &pthread_args[th]);
    }
    for (int th = 0; th < nthreads; th++)
        pthread_join(pthread_args[th].tid, NULL);
}

Point* dsort(Point *a, int n, int sorted_size)
{
    Point *res = new Point[n];
    for (int i = 0; i < n; i++)
        res[i] = a[i];
    return res;
}

void dhsort(Point *a, int n)
{
    int nthreads = n / MAX_HSORT_ELEMS;
    nthreads = nthreads ? nthreads : MIN_THREADS_NUM;
    hsort_threads(a, n, nthreads);

    /*
    int full = n % nthreads;
    int a1_sorted_size = n/nthreads + 1;
    int a1_elems = full*a1_sorted_size;
    int a2_sorted_size = a1_sorted_size - 1;
    int a2_elems = (nthreads - full)*a2_sorted_size;
    Point *a1 = dsort(a, a1_elems, a1_sorted_size);
    Point *a2 = dsort(a + a1_elems, a2_elems, a2_sorted_size);
    int i1 = 0, i2 = 0;
    for (int i = 0; i < n; i++) {
        if (a1[i1] < a2[i2])
            a[i] = a1[i1++];
        else
            a[i] = a2[i2++];
    }*/
    return;
}
