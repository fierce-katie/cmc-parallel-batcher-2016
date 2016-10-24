#include <stdio.h>
#include <vector>

using namespace std;

typedef pair<int, int> comparator;

void print_vector(vector<int> v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

void print_comparators(vector<comparator> cmp)
{
    vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++)
        printf("%d %d\n", it->first, it->second);
    printf("%lu\n", cmp.size());
}

void join(vector<int> idx_up, int n0, vector<int> idx_down, int n1,
          vector<comparator> &cmp)
{
    int n = n0 + n1;
    if (n == 1)
        return;

    if (n == 2) {
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
    //printf("Sorting idx:\n");
    //print_vector(idx, n);
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
    print_comparators(cmp);
    return 0;
}


