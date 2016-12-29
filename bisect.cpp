// Catherine Galkina, group 524, year 2016
// File bisect.cpp
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cstddef>
#include <float.h>

#define MIN_THREADS_NUM 4
#define MAX_HSORT_ELEMS 100000

bool axis = true;
int *domain_array;
MPI_Datatype MPI_POINT;

struct Domain {
    float coord[2];
    int index;
    int domain;
};

struct Point {
    float coord[2];
    int index;
};

float x(int i, int j)
{
    return i;
    //return 10*(float)rand()/(float)(RAND_MAX/(i*j+1));
}

float y(int i, int j)
{
    return j;
    //return 10*(float)rand()/(float)(RAND_MAX/(i*j+1));
}

Point* init_points(int n, int ny, int procs, int proc_elems, int rank)
{
    Point *res = new Point[proc_elems];
    int tmp = n/procs;
    int not_fake = n % procs;
    int real_elems = rank < not_fake ? tmp + 1 : tmp;
    int delta;
    if (rank < not_fake)
        delta = rank*proc_elems;
    else
        delta = not_fake*proc_elems + tmp*(rank - not_fake);
    for (int k = 0; k < real_elems; k++) {
        int i = (k + delta) / ny, j = (k + delta) % ny;
        Point p;
        p.coord[0] = x(i, j);
        p.coord[1] = y(i, j);
        p.index = i*ny + j;
        res[k] = p;
    }
    for (int k = real_elems; k < proc_elems; k++) {
        Point p;
        p.coord[0] = FLT_MAX;
        p.coord[1] = FLT_MAX;
        p.index = -1;
        res[k] = p;
    }

    return res;
}

MPI_Datatype pointType()
{
    MPI_Datatype point;
    MPI_Datatype types[2] = { MPI_FLOAT, MPI_INT };
    int blocks[2] = { 2, 1 };
    MPI_Aint disps[2] = { offsetof(Point, coord), offsetof(Point, index) };
    MPI_Type_create_struct(2, blocks, disps, types, &point);
    MPI_Type_commit(&point);
    return point;
}

typedef std::pair<int, int> comparator;

void swap(comparator cmp, std::vector<int> &v)
{
    int fst = cmp.first;
    int snd = cmp.second;
    if (v[fst] > v[snd]) {
        int tmp = v[fst];
        v[fst] = v[snd];
        v[snd] = tmp;
    }
}

void print_vector(std::vector<int> &v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

void print_comparators(std::vector<comparator> &cmp)
{
    std::vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++)
        printf("%d %d\n", it->first, it->second);
    printf("%lu\n", cmp.size());
}

void print_points(Point* p, int n, int rank, const char *comment)
{
    for (int i = 0; i < n; i++) {
        Point point = p[i];
        if (point.index >= 0)
            printf("%d %d: %f %f %s\n", rank, point.index, point.coord[0],
                   point.coord[1], comment);
    }
}

bool check_args(int argc, char **argv, int &nx, int &ny, int &k)
{
    if (argc < 4) {
        printf("Wrong arguments. Usage: bisect nx ny k\n");
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
    check = sscanf(argv[3], "%d", &k);
    if (!check) {
        printf("k must be int: %s\n", argv[3]);
        return false;
    }
    if (!((nx >= 1) && (ny >= 1) && (k >= 1))) {
        printf("Wrong nx or nx or k\n");
        return false;
    }
    return true;
}

int count_tacts(int n, std::vector<comparator> &cmp)
{
    std::vector<int> v(n);
    std::vector<comparator>::iterator it;
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

inline int compare_points(const void *a, const void *b)
{
    Point *aptr = (Point * const)a;
    Point *bptr = (Point * const)b;
    int c = axis ? 0 : 1;
    float ax = aptr->coord[c];
    float bx = bptr->coord[c];

    if (ax == bx)
        return 0;
    else if (ax > bx)
        return 1;
    return -1;
}

struct PthreadArgs {
    pthread_t tid;
    Point *a;
    int n;
};

void heapify(Point *a, int i, int n)
{
    int imax, l, r;
    Point tmp;
    while (i < n) {
        imax = i;
        l = 2*i + 1;
        r = l + 1;
        if (l < n && (compare_points(&a[l], &a[imax]) > 0))
            imax = l;
        if (r < n && (compare_points(&a[r], &a[imax]) > 0))
            imax = r;
        if (imax == i)
            return;
        tmp = a[i];
        a[i] = a[imax];
        a[imax] = tmp;
        i = imax;
    }
}

void make_heap(Point *a, int n)
{
    for (int i = n/2 - 1; i >= 0; i--)
        heapify(a, i, n);
}

void hsort(Point *a, int n)
{
    make_heap(a, n);
    Point tmp;
    for (int i = n - 1; i > 0; i--) {
        tmp = a[0];
        a[0] = a[i];
        a[i] = tmp;
        heapify(a, 0 ,i);
    }
}

void* hsort_start(void *arg)
{
    PthreadArgs *args = (PthreadArgs*)arg;
    hsort(args->a, args->n);
    return NULL;
}


void hsort_threads(Point *a, int n, int nthreads)
{
    PthreadArgs *pthread_args = new PthreadArgs[nthreads];
    int tmp = ceil(n / (double)nthreads);
    int elems;
    int offset = 0;
    for (int th = 0; th < nthreads; th++) {
        if (n - offset >= tmp)
            elems = tmp;
        else
            elems = (n - offset > 0) ? n - offset : 0;

        pthread_args[th].a = a + offset;
        pthread_args[th].n = elems;
        pthread_create(&pthread_args[th].tid, NULL, hsort_start,
                       &pthread_args[th]);
        offset += elems;
    }
    for (int th = 0; th < nthreads; th++)
        pthread_join(pthread_args[th].tid, NULL);
    delete [] pthread_args;
}

void dsort(Point *array, int n, int sorted)
{
    Point *a = array;
    Point *b = new Point[n];
    Point *c;

    for (int i = sorted; i < n ; i *= 2) {
        for (int j = 0; j < n; j = j + 2*i) {
            int r = j + i;

            int n1 = (i < n - j) ? i : n - j;
            int n2 = (i < n - r) ? i : n - r;
            n1 = (n1 < 0) ? 0 : n1;
            n2 = (n2 < 0) ? 0 : n2;

            for (int ia = 0, ib = 0, k = 0; k < n1 + n2; k++) {
                if (ia >= n1)
                    b[j+k] = a[r+ib++];
                else if (ib >= n2)
                    b[j+k]=a[j+ia++];
                else if (compare_points(&a[j+ia], &a[r+ib]) < 0)
                    b[j+k]=a[j+ia++];
                else
                    b[j+k]=a[r+ib++];
            }
        }
        c = a;
        a = b;
        b = c;
    }

    c = a;
    a = b;
    b = c;
    if (b != array) {
        memcpy(array, b, n*sizeof(*array));
        delete [] b;
    } else {
        delete [] a;
    }
}

void dhsort(Point *a, int n)
{
    int nthreads = ceil(n / (double)MAX_HSORT_ELEMS);
    nthreads = nthreads > MIN_THREADS_NUM ? nthreads : MIN_THREADS_NUM;
    hsort_threads(a, n, nthreads);
    dsort(a, n, ceil(n / (double)nthreads));

    return;
}

void join(std::vector<int> idx_up, int n0, std::vector<int> idx_down, int n1,
          std::vector<comparator> &cmp)
{
    int n = n0 + n1;
    if (n == 1)
        return;

    if (n0 == 1 && n1 == 1) {
        cmp.push_back(comparator(idx_up[0], idx_down[0]));
        return;
    }

    int n0_even = n0/2;
    int n0_odd = n0 - n0_even;
    std::vector<int> idx_up_even(n0_even);
    std::vector<int> idx_up_odd(n0_odd);

    int n1_even = n1/2;
    int n1_odd = n1 - n1_even;
    std::vector<int> idx_down_even(n1_even);
    std::vector<int> idx_down_odd(n1_odd);

    std::vector<int> idx_result;

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

void sort(std::vector<int>idx, int n, std::vector<comparator> &cmp)
{
    if (n == 1) {
        return;
    }

    int n0 = n/2;
    int n1 = n - n0;

    std::vector<int> idx_up;
    std::vector<int> idx_down;

    int i;
    for (i = 0; i < n0; i++)
        idx_up.push_back(idx[i]);
    for (i = n0; i < n; i++)
        idx_down.push_back(idx[i]);

    sort(idx_up, n0, cmp);
    sort(idx_down, n1, cmp);
    join(idx_up, n0, idx_down, n1, cmp);
}

void make_comparators(int procs, std::vector<comparator> &cmp)
{
    std::vector<int> idx;
    for (int i = 0; i < procs; i++)
        idx.push_back(i);
    sort(idx, procs, cmp);
    return;
}

void batcher(Point* &proc_points, int proc_elems, std::vector<comparator> cmp,
             MPI_Comm comm)
{
    dhsort(proc_points, proc_elems);

    // Exchanging elements
    Point *tmp_points = new Point[proc_elems];
    Point *other_points = new Point[proc_elems];
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;
    std::vector<comparator>::iterator it;
    for (it = cmp.begin(); it != cmp.end(); it++) {
        if (rank == it->first) {
            MPI_Send(proc_points, proc_elems, MPI_POINT,
                     it->second, 0, comm);
            MPI_Recv(other_points, proc_elems, MPI_POINT,
                     it->second, 0, comm, &status);
            int idx = 0;
            int other_idx = 0;
            for (int tmp_idx = 0; tmp_idx < proc_elems; tmp_idx++) {
                Point my = proc_points[idx];
                Point other = other_points[other_idx];
                if (my.coord[0] < other.coord[0]) {
                    tmp_points[tmp_idx] = my;
                    idx++;
                } else {
                    tmp_points[tmp_idx] = other;
                    other_idx++;
                }
            }
            swap_ptr(&proc_points, &tmp_points);
        }

        if (rank == it->second) {
            MPI_Recv(other_points, proc_elems, MPI_POINT,
                     it->first, 0, comm, &status);
            MPI_Send(proc_points, proc_elems, MPI_POINT,
                     it->first, 0, comm);
            int idx = proc_elems - 1;
            int other_idx = proc_elems - 1;
            for (int tmp_idx = proc_elems - 1; tmp_idx >= 0; tmp_idx--) {
                Point my = proc_points[idx];
                Point other = other_points[other_idx];
                if (my.coord[0] > other.coord[0]) {
                    tmp_points[tmp_idx] = my;
                    idx--;
                } else {
                    tmp_points[tmp_idx] = other;
                    other_idx--;
                }
            }
            swap_ptr(&proc_points, &tmp_points);
        }
    }
}

void bisect_seq(Point *points, int n0, int n, int dom0, int k)
{
    // One point
    if (n == 1) {
      domain_array[n0] = dom0;
      return;
    }

    // One domain
    if (k == 1) {
      for (int i = 0; i < n; i++)
        domain_array[n0 + i] = dom0;
      return;
    }

    // Sort and change axis
    dhsort(points + n0, n);
    axis = !axis;

    // Split ratio
    int k1 = (k + 1) / 2;
    int k2 = k - k1;
    int n1 = n*k1/(double)k;
    int n2 = n - n1;

    // Recursively bisect parts
    bisect_seq(points, n0, n1, dom0, k1);
    bisect_seq(points, n0 + n1, n2, dom0 + k1, k2);
}

int remove_fake(Point **points, int n)
{
    Point *res_points = new Point[n];
    int res_n = 0;
    for (int i = 0; i < n; i++) {
        if ((*points)[i].index != -1)
          res_points[res_n++] = (*points)[i];
    }
    delete [] (*points);
    *points = res_points;
    return res_n;
}

double min_coord(Point *points, int n, int c)
{
    double res = FLT_MAX;
    for (int i = 0; i < n; i++)
      if ((points[i].index != -1) && (points[i].coord[c] < res))
        res = points[i].coord[c];
    return res;
}

double max_coord(Point *points, int n, int c)
{
    double res = FLT_MIN;
    for (int i = 0; i < n; i++)
      if ((points[i].index != -1) && (points[i].coord[c] > res))
        res = points[i].coord[c];
    return res;
}

void change_axis(Point *points, int n, MPI_Comm comm)
{
    double local_minx = min_coord(points, n, 0);
    double local_miny = min_coord(points, n, 1);
    double local_maxx = max_coord(points, n, 0);
    double local_maxy = max_coord(points, n, 1);
    double minx = 0, miny = 0, maxx = 0, maxy = 0;
    MPI_Allreduce(&local_minx, &minx, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_miny, &miny, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_maxx, &maxx, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&local_maxy, &maxy, 1, MPI_DOUBLE, MPI_MAX, comm);
    double dx = maxx - minx, dy = maxy - miny;
    axis = dx <= dy;
}

Point *add_fake(Point *points, int n, int procs)
{
    int proc_elems = ceil(n/(double)procs);
    Point *res = new Point[procs*proc_elems];
    int idx = 0;
    for (int i = 0; i < n; i++, idx++) {
        res[idx] = points[i];
        points[i].index = -1;
        points[i].coord[0] = FLT_MAX;
        points[i].coord[1] = FLT_MAX;
    }
    for (; idx < procs*proc_elems; idx++) {
        res[idx].index = -1;
        res[idx].coord[0] = FLT_MAX;
        res[idx].coord[1] = FLT_MAX;
    }
    return res;
}

void bisect(Point **points, int &proc_elems, int n, int k, int dom0,
            MPI_Comm comm)
{
    // One domain
    if (k == 1) {
        int real_elems = remove_fake(points, proc_elems);
        domain_array = new int[real_elems];
        for (int i = 0; i < real_elems; i++)
            domain_array[i] = dom0;
        proc_elems = real_elems;
        return;
    }

    int rank, procs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    // One proc
    if (procs == 1) {
        int real_elems = remove_fake(points, proc_elems);
        domain_array = new int[real_elems];
        bisect_seq(*points, 0, real_elems, dom0, k);
        proc_elems = real_elems;
        return;
    }

    // Split ratio
    int k1 = (k + 1) / 2;
    int k2 = k - k1;
    int n1 = n*k1/(double)k;
    int n2 = n - n1;
    int middle = n1 % proc_elems;
    int procs1 = n1 / proc_elems; // on the left
    int psplit = procs1 ? (rank >= procs1 ? 0 : 1) : (rank > procs1 ? 0 : 1);

    change_axis(*points, proc_elems, comm);

    // Sort
    std::vector<comparator> cmp;
    make_comparators(procs, cmp);
    batcher(*points, proc_elems, cmp, comm);

    // Split procs
    MPI_Comm comm_split;
    MPI_Comm_split(comm, psplit, rank, &comm_split);

    // Rearrange elements on edges for recursive call
    MPI_Status s;
    if (procs1) {
        if (rank <= procs1) {
            int proc_elems1 = ceil(middle/(double)procs1);
            if (rank == procs1) { // send to the left
                Point *tmp = add_fake(*points, middle, procs1);
                for (int i = 0; i < procs1; i++)
                    MPI_Send(tmp + i*proc_elems1, proc_elems1, MPI_POINT, i, 0,
                             comm);
                delete [] tmp;
                bisect(points, proc_elems, n2, k2, dom0 + k1, comm_split);
            } else { // recieve new elems
                Point *other = new Point[proc_elems + proc_elems1];
                MPI_Recv(other + proc_elems, proc_elems1, MPI_POINT, procs1, 0,
                         comm, &s);
                memcpy(other, *points, proc_elems * sizeof(Point));
                proc_elems += proc_elems1;
                delete [] (*points);
                *points = other;
                bisect(points, proc_elems, n1, k1, dom0, comm_split);
            }
        } else { // nothing to exchange, just decompose
            bisect(points, proc_elems, n2, k2, dom0 + k1, comm_split);
        }
    } else { // we'll have to send some to the right
        int proc_elems1 = ceil((proc_elems - middle)/(double)(procs - procs1));
        if (rank) { // recieve new elems
            Point *other = new Point[proc_elems + proc_elems1];
            MPI_Recv(other + proc_elems, proc_elems1, MPI_POINT, procs1, 0,
                     comm, &s);
            memcpy(other, *points, proc_elems * sizeof(Point));
            proc_elems += proc_elems1;
            delete [] (*points);
            *points = other;
            bisect(points, proc_elems, n2, k2, dom0 + k1, comm_split);
        } else { //rank == 0
            Point *tmp = add_fake(*points + middle, proc_elems - middle,
                                  procs - procs1 - 1);
            for (int i = procs1 + 1, j = 0; i < procs; i++, j++)
                MPI_Send(tmp + j*proc_elems1, proc_elems1, MPI_POINT, i, 0,
                         comm);
            delete [] tmp;
            bisect(points, proc_elems, n1, k1, dom0, comm_split);
        }
    }
}

bool connected(Domain p1, Domain p2, int ny)
{
    int i1 = p1.index / ny;
    int j1 = p1.index % ny;
    int i2 = p2.index / ny;
    int j2 = p2.index % ny;
    return (((i1 == i2) && (abs(j1 - j2) == 1)) ||
            ((j1 == j2) && (abs(i1 - i2) == 1)));
}

int edges(Domain *p, int n, int ny)
{
    int res = 0;
    for (int i = p[0].domain, j = 0; j < n; i++)
        for (; p[j].domain == i && j < n; j++) {
            Domain p1 = p[j];
            for (int k = j + 1; p[k].domain == i && k < n; k++) {
                Domain p2 = p[k];
                if (connected(p1, p2, ny))
                    res++;
            }
        }
    return res;
}
int main(int argc, char **argv)
{
    int nx, ny, k;

    // Parsing command line arguments
    if (!check_args(argc, argv, nx, ny, k))
        return 1;

    MPI_Init(&argc, &argv);

    int rank, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    // Calculating elems per processor
    int n = nx*ny;
    int fake = n % procs ? (procs - n % procs) : 0;
    int elems = n + fake;
    int proc_elems = elems / procs;

    // Initializing points
    srand(time(NULL) + rank);
    Point *proc_points =
        init_points(n, ny, procs, proc_elems, rank);

    print_points(proc_points, proc_elems, rank, "initial");

    // Decomposition
    MPI_POINT = pointType();
    double start_time = MPI_Wtime();
    bisect(&proc_points, proc_elems, n, k, 0, MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double time = end_time - start_time;
    double max_time = 0;
    MPI_Reduce(&time, &max_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!rank)
        printf("Decomposition time: %f sec.\n", max_time);

    /*
    for (int i = 0; i < proc_elems; i++) {
        Point p = proc_points[i];
        printf("%d %d %d %f %f %d\n", p.index, p.index / ny, p.index % ny, p.coord[0], p.coord[1], domain_array[i]);
    }*/

    // Remove fake
    Domain *proc_domains = new Domain[proc_elems];
    int new_proc_elems = 0;
    for (int i = 0; i < proc_elems; i++) {
        if (proc_points[i].index != -1) {
            proc_domains[new_proc_elems].index = proc_points[i].index;
            proc_domains[new_proc_elems].coord[0] = proc_points[i].coord[0];
            proc_domains[new_proc_elems].coord[1] = proc_points[i].coord[1];
            proc_domains[new_proc_elems++].domain = domain_array[i];
        }
    }

    // Cut edges
    int total_edges = nx*(ny - 1) + ny*(nx - 1);
    int local_edges = edges(proc_domains, new_proc_elems, ny);
    int sum_edges = 0;
    MPI_Reduce(&local_edges, &sum_edges, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    int cut = total_edges - sum_edges;

    if (!rank)
        printf("Edges cut: %d of %d\n", cut, total_edges);

    delete [] proc_points;
    delete [] proc_domains;
    delete [] domain_array;
    MPI_Finalize();
    return 0;
}

