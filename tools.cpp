// Catherine Galkina, group 524, year 2016
// File tools.cpp
#include <stdio.h>
#include "tools.h"

using namespace std;


void swap(comparator cmp, vector<int> &v)
{
    int fst = cmp.first;
    int snd = cmp.second;
    if (v[fst] > v[snd]) {
        int tmp = v[fst];
        v[fst] = v[snd];
        v[snd] = tmp;
    }
}

void print_vector(vector<int> &v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

void print_comparators(vector<comparator> &cmp)
{
    vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++)
        printf("%d %d\n", it->first, it->second);
    printf("%lu\n", cmp.size());
}

void print_points(Point* p, int n, int rank, const char *comment)
{
    for (int i = 0; i < n; i++) {
        Point point = p[i];
        if (point.GetIndex() >= 0)
            printf("%d %d: %f %s\n", rank, point.GetIndex(), point.GetX(),
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

#if 0
bool check_args(int argc, char **argv, int &n0, int &n1)
{
    if (argc < 2) {
        printf("Wrong arguments. Usage: bsort n0 [n1]\n");
        return false;
    }
    int check = sscanf(argv[1], "%d", &n0);
    if (!check) {
        printf("n0 must be int: %s\n", argv[1]);
        return false;
    }
    if (argc >= 3) {
        check = sscanf(argv[2], "%d", &n1);
        if (!check) {
            printf("n1 must be int: %s\n", argv[2]);
            return false;
        }
    } else
        n1 = 0;
    if (!((n0 >= 1) && (n1 >= 0))) {
        printf("Wrong n0 or n1\n");
        return false;
    }
    return true;
}
#endif

int count_tacts(int n, std::vector<comparator> &cmp)
{
    vector<int> v(n);
    vector<comparator>::iterator it;
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

