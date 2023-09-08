/*
  TASK: Calculate a * X + Y using OpenMP offloading
*/
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static double wallclock() {
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

double* saxpy(int n, double a, double* x, double* y) {
#pragma acc kernels
    for (int i = 0; i < n; i += 1) {
        y[i] = a * x[i] + y[i];
    }
    return y;
}

void printvec(int n, double* vec) {
    for (int i = 0; i < n; i += 1) {
        printf("%f\n", vec[i]);
    }
}

int main(int argc, char* argv[]) {
    // Your code goes here!
    int N = 1 << 30;
    double *x, *y;
    double t_start;

    x = (double*)malloc(N * sizeof(double));
    y = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i += 1) {
        x[i] = (double)i;
        y[i] = 1.0;
    }

    printf("Start SAXPY of length %d\n", N);

    t_start = wallclock();

    saxpy(N, 2.0, x, y);
    // printvec(N, saxpy(N, 2.0, x, y));

    printf("End SAXPY\nTime to execute: %10.3fs\n", wallclock() - t_start);

    // printvec(N, saxpy(N, 2.0, x, y));
    return (0);
}