#include <stdio.h>
#include <cmath>
#include <vector>

using namespace std;

void print_vector(vector<int> &v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

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

int main(int argc, char **argv)
{
    int n;
    sscanf(argv[1], "%d", &n);
    vector<vector<int> > tests;
    generate_tests(n, tests);
    vector<vector<int> >::iterator it;
    for (it = tests.begin(); it != tests.end(); it++) {
        print_vector(*it, n);
        if (is_sorted_binary(*it, n))
            putchar('+');
        else
            putchar('-');
        putchar('\n');
    }
    return 0;
}
