
#include <stdio.h>

int sum_of_int(int a, int b)
{
    int c = a + b;
    printf("sum of %d and %d is %d\n", a, b, c);
    return c;
}

double sum_of_double(double a, double b)
{
    double c = a + b;
    printf("sum of %g and %g is %g\n", a, b, c);
    return c;
}
