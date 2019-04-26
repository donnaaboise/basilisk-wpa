//#include "grid/multigrid1D.h"
#include "grid/bitree.h"
//#include "utils.h"   /* Needed for DT */

#include "waveprop.h"

#include "rp1_swe.h"

scalar h[];
vector hu[];
scalar *scalars = {h};
vector *vectors = {hu};

int mwaves = 2;
int limiters[2] = {1,1}, *mthlim = limiters;

bool dt_fixed = false;
double dt_initial = 1e-3;


vector *aux = NULL;
int maux = 0;

#if TREE
int matlab_out = false;
#else
int matlab_out = false;
#endif


#define MAXLEVEL 15
#define MINLEVEL 15

double a = 0.1, b = 12.2;

double sum0;

int main() 
{
    origin(0);
    L0 = 4000;

    CFL = 0.9;

    N = 1 << MINLEVEL;
    run();
}

event defaults(i=0)
{
}


event init (i = 0) 
{
    wpa_rp1 = rp1_swe;   
    foreach() 
    {
        h[] = 1.0 + a*(fabs(x) < b);
        hu.x[] = 0;
    }
    dt_initial = 0.1047505286036417;   /* to get CFL= 0.9 on first step */
}

#if 0
event verbosity(i++)
{
    printf("Step %6d at t = %20.16f\n",i,t);
}
#endif

event plot (t = {0, 35, 90, 700, 2000})
{
    if (!matlab_out)
    {
        char name[80];
        sprintf (name, "t-%d", Frame);
        FILE * fp = fopen (name, "w");
        fprintf(fp,"%20.16f %20.16f %20.16f\n",t,0.,0.);
        foreach()
        {
            fprintf (fp, "%20.16f %24.16e %24.16e\n", x, h[], hu.x[]);
        }
        fclose (fp);
    }
}


event conservation_check(t = {0,2000})
{
    double sum = 0;
    foreach()
    {
        sum += Delta*h[];
    }
    if (i == 0)
        sum0 = sum;

    printf("t = %20.16f; mass = %24.16e; diff = %24.16e\n",t, sum, fabs(sum0-sum));
    fflush(stdout);
}



/* 1e-3 is too large */
#if TREE
event adapt (i++) {
  adapt_wavelet ({h}, (double []){1e-5}, maxlevel = MAXLEVEL);
}
#endif







