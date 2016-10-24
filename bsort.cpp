#include <stdio.h>
#include <string.h>

void join(int *idx_up, int n0, int *idx_down, int n1, int &cmp)
{
    int n = n0 + n1;
    if (n == 1)
        return;

    if (n == 2) {
        printf("%d %d\n", idx_up[0], idx_down[0]);
        cmp++;
        return;
    }

    int n0_even = n0/2;
    int n0_odd = n0 - n0_even;
    int* idx_up_even = new int[n0_even];
    int* idx_up_odd = new int[n0_odd];

    int n1_even = n1/2;
    int n1_odd = n1 - n1_even;
    int* idx_down_even = new int[n1_even];
    int* idx_down_odd = new int[n1_odd];

    int* idx_result = new int[n];

    int i0 = 0, i1 = 0;
    for (int i = 0; i < n0; i++)
        if (i%2) {
            idx_up_even[i0] = idx_up[i];
            i0++;
        } else {
            idx_up_odd[i1] = idx_up[i];
            i1++;
        }
    i0 = i1 = 0;
    for (int i = 0; i < n1; i++)
        if (i%2) {
            idx_down_even[i0] = idx_down[i];
            i0++;
        } else {
            idx_down_odd[i1] = idx_down[i];
            i1++;
        }

    join(idx_up_odd, n0_odd, idx_down_odd, n1_odd, cmp);
    join(idx_up_even, n0_even, idx_down_even, n1_even, cmp);

    memmove(idx_result, idx_up, n0*sizeof(int));
    memmove(idx_result + n0, idx_down, n1*sizeof(int));

    for (int i = 1; i < n - 1; i += 2) {
        printf("%d %d\n", idx_result[i], idx_result[i + 1]);
        cmp++;
    }

    /*delete [] idx_up_even;
    delete [] idx_up_odd;
    delete [] idx_down_even;
    delete [] idx_down_odd;
    delete [] idx_result;
    */
}

void sort(int *idx, int n, int &cmp)
{
    //printf("Sorting idx:\n");
    //print_array(idx, n);
    if (n == 1) {
        return;
    }

    int n0 = n/2;
    int *idx_up = new int[n0];
    int n1 = n - n0;
    int *idx_down = new int(n1);

    memmove(idx_up, idx, n0*sizeof(int));
    memmove(idx_down, idx + n0, n1*sizeof(int));

    sort(idx_up, n0, cmp);
    sort(idx_down, n1, cmp);
    join(idx_up, n0, idx_down, n1, cmp);

    //delete [] idx_up;
    //delete [] idx_down;
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
    int cmp = 0;
    int n0, n1;

    // Parsing command line arguments
    if (!check_args(argc, argv, n0, n1))
        return 1;

    // Sorting or joining
    printf("%d %d %d\n", n0, n1, 0);
    if (n1) {
        int *idx_up = new int[n0];
        int *idx_down = new int[n1];
        int i;
        for (i = 0; i < n0; i++)
            idx_up[i] = i;
        for (i = n0; i < n0 + n1; i++)
            idx_down[i] = i;
        join(idx_up, n0, idx_down, n1, cmp);
    } else {
        int *idx = new int[n0];
        for (int i = 0; i < n0; i++)
            idx[i] = i;
        sort(idx, n0, cmp);
    }
    printf("%d\n", cmp);
    return 0;
}


