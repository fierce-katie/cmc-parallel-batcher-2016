#include <vector>
#include <cmath>

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

