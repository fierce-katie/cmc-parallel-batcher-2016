#include <stdio.h>

void join(int start0, int n0, int start1, int n1, int step,
          int &cmp, int &t)
{
    int n = n0 + n1;
    if (n > 2) {
        join(start0, (n0 + 1)/2,
             start1, (n1 + 1)/2,
             step*2, cmp, t);
        join(start0 + step, n0 - (n0 + 1)/2,
             start1 + step, n1 - (n1 + 1)/2,
             step*2, cmp, t);
    }
    if (n == 2) {
        printf("%d %d\n", start0, start1);
        cmp++;
    }
    if (n > 2)
        for (int i = start0 + step; i + step < start1 + n1*step; i += step*2) {
            printf("%d %d\n", i, i + step);
            cmp++;
        }
    return;
}

void sort(int start, int step, int n,
          int &cmp, int &t)
{
    int n0, n1;
    n0 = n/2;
    n1 = n - n0;
    if (n >= 2) {
        sort(start, step, n0, cmp, t);
        sort(start + n0*step, step, n1, cmp, t);
        join(start, n0, start + n0, n1, step, cmp, t);
    }
    return;
}

int main(int argc, char **argv)
{
    int cmp = 0, t = 0;
    int n0, n1;

    // Parsing command line arguments
    if (argc < 2) {
        printf("Wrong arguments. Usage: bsort n0 [n1]\n");
        return 1;
    }
    int check = sscanf(argv[1], "%d", &n0);
    if (!check) {
        printf("n0 must be int: %s\n", argv[1]);
        return 1;
    }
    if (argc >= 3) {
        check = sscanf(argv[2], "%d", &n1);
        if (!check) {
            printf("n1 must be int: %s\n", argv[2]);
            return 1;
        }
    } else
        n1 = 0;
    if (!((n0 >= 1) && (n1 >= 0))) {
        printf("Wrong n0 or n1\n");
        return 1;
    }

    // Sorting or joining
    printf("%d %d %d\n", n0, n1, 0);
    if (n1)
        join(0, n0, n0, n1, 1, cmp, t);
    else
        sort(0, 1, n0, cmp, t);
    printf("%d\n", cmp);
    return 0;
}


