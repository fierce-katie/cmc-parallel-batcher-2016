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

