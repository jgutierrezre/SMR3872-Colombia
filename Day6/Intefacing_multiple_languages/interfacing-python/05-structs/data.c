
#include<stdio.h>

/* global variables */
int    inum = -1;
double dnum = 1.1;

struct parm {
   int  type;
   char *label;
   double epsilon, sigma;
};

void pass_by_value(struct parm p)
{
    printf("type=%d  label=%s epsilon=%g  sigma=%g\n",
           p.type, p.label, p.epsilon, p.sigma);
}

void pass_by_reference(struct parm *p)
{
    printf("type=%d  label=%s epsilon=%g  sigma=%g\n",
           p->type, p->label, p->epsilon, p->sigma);
}
