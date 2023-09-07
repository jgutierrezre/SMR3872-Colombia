#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

static double get_wtime()
{
#if !defined(_OPENMP) || (defined(_OPENMP) && (_OPENMP <= 201511))
    double time;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    time = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
    return time;
#else
    return omp_get_wtime();
#endif
}

int main(int argc, char *argv[])
{
    int n, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;
    double start;

#if defined(_OPENMP)
    numprocs = omp_get_max_threads();
#else
    numprocs = 1;
#endif

    if (argc > 1) n = atoi(argv[1]);
    else n = 0;

    if (n == 0) {
        printf("\nUsage: %s <num intervals>\n\n",argv[0]);
        return 1;
    }

    start = get_wtime();
    h   = 1.0 / (double) n;
    sum = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for private(i,x) shared(h,n) reduction(+:sum)
#endif
    for (i = 1; i <= n; ++i) {
        x = h * ((double)i - 0.5);
        sum += (4.0 / (1.0 + x*x));
    }
    pi = h * sum;
    printf("\npi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
    printf("Time: %.16f seconds on %d CPUs\n\n",get_wtime()-start,numprocs);
    return 0;
}
