#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define CHUNK_SIZE (1 << 20)

static double wallclock() {
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

double* saxpy(int n, double a, double* x, double* y) {
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
    int N = 1 << 30;
    double *x, *y;
    double t_start;

    x = (double*)malloc(CHUNK_SIZE * sizeof(double));
    y = (double*)malloc(CHUNK_SIZE * sizeof(double));

    printf("Start SAXPY of length %d\n", N);

    t_start = wallclock();

    for (int i = 0; i < N; i += CHUNK_SIZE) {
        int chunk_size = CHUNK_SIZE;
        if (i + chunk_size > N) {
            chunk_size = N - i;
        }
        for (int j = 0; j < chunk_size; j += 1) {
            x[j] = (double)(i + j);
            y[j] = 1.0;
        }
        saxpy(chunk_size, 2.0, x, y);
    }

    printf("End SAXPY\nTime to execute: %10.3fs\n", wallclock() - t_start);

    return (0);
}
