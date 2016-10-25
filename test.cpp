// Catherine Galkina, group 524, year 2016
// File test.cpp
#include <vector>
#include <cmath>
#include <stdio.h>

#include "test.h"

using namespace std;

bool is_sorted_binary(vector<int> &v, int n)
{
    for (int i = 0; i < n; i++)
        if (v[i]) {
            for (int j = i; j < n; j++)
                if (!v[j]) return false;
            break;
        }
    return true;
}

void to_binary_vector(int n, int total, vector<int> &res)
{
    int last = total - 1;
    while (n) {
        res[last--] = n%2;
        n /= 2;
    }
}

void generate_tests(int n, vector<vector<int> > &tests)
{
    int test_num = pow(2, n);
    for (int i = 0; i < test_num; i++) {
        vector<int> test(n);
        to_binary_vector(i, n, test);
        tests.push_back(test);
    }
}

void do_sort(vector<int> &v, vector<comparator> &cmp)
{
    vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++)
        swap(*it, v);
    return;
}

void run_tests(int n, vector<comparator> &cmp) {
    printf("\nRunning tests...\n");
    vector<vector<int> > tests;
    generate_tests(n, tests);
    vector<vector<int> >::iterator it;
    for (it = tests.begin(); it != tests.end(); it++) {
        do_sort(*it, cmp);
        if (!is_sorted_binary(*it, n)) {
            printf("Test failed!\n");
            return;
        }
    }
    printf("All tests have passed.\n");
    return;
}
