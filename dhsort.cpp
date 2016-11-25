// Catherine Galkina, group 524, year 2016
// File dhsort.cpp

#include "dhsort.h"
#include "tools.h"

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

void hsort(Point* a, int n)
{
    make_heap(a, n);
    Point tmp;
    for (int i = n - 1; i > 0; i--) {
        tmp = a[0];
        a[0] = a[i];
        a[i] = tmp;
        heapify(a, 0 ,i);
    }
    return;
}

void dhsort(Point* a, int n)
{
    int nthreads = 4;
    int tmp = n/nthreads;
    int full = n % nthreads;
    for (int th = 0; th < nthreads; th++) {
        int elems = th < full ? tmp + 1 : tmp;
        int offset;
        if (th < full)
            offset = th*(tmp + 1);
        else
            offset = full*(tmp + 1) + (th - full)*tmp;
        printf("th = %d, elems = %d, offset = %d\n", th, elems, offset);
        print_points(a + offset, elems, 0, "dh-initial");
        hsort(a + offset, elems);
        print_points(a + offset, elems, 0, "dh-sorted");
    }
}

