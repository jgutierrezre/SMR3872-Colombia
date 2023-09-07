#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int n, myid, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double mypi, pi, h, sum, x;
    double start;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    if (argc > 1) n = atoi(argv[1]);
    else n = 0;

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n == 0) {
        if (myid == 0) {
            printf("\nUsage: %s <num intervals>\n\n",argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    start = MPI_Wtime();
    h   = 1.0 / (double) n;
    sum = 0.0;
    for (i = myid + 1; i <= n; i += numprocs) {
        x = h * ((double)i - 0.5);
        sum += (4.0 / (1.0 + x*x));
    }
    mypi = h * sum;
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0) {
         printf("\npi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
         printf("Time: %.16f seconds on %d CPUs\n\n",MPI_Wtime()-start,numprocs);
    }
    MPI_Finalize();
    return 0;
}
