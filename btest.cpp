#include "tools.h"
#include "test.h"

using namespace std;

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
