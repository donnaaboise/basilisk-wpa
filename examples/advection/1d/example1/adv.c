//#include "grid/multigrid1D.h"
#include "grid/bitree.h"

#include "waveprop.h"
#include "rp1_adv.h"

extern double ubar;


scalar q[];
scalar *scalars = {q};
vector *vectors = NULL;

double dt_fixed;

int mwaves = 1;

int limiters[1] = {3}, *mthlim = limiters;

vector *aux = NULL;
int maux = 0;

#if TREE
int matlab_out = false;
#else
int matlab_out = false;
#endif

#define MAXLEVEL 15

int main() 
{
    origin (-1.0);
    L0 = 2;

    CFL = 0.9;

    periodic (right);        

    N = 1 << MAXLEVEL;
    run();
}

event defaults(i=0)
{
    dt = 1e-4;
    dt_fixed = DT;
    claw_rp1 = rp1_adv;
}


event init (i = 0) 
{
    foreach() 
    {
        //q[] = exp(-100.0*sq(x));     
        //q[] = sin(4*pi*x);
        q[] = fabs(x) <= 0.25;
    }
}

event plot (t+=0.5; t<=4)
{
    if (!matlab_out)
    {        
        char name[80];
        sprintf (name, "t-%d", Frame);
        FILE * fp = fopen (name, "w");
        fprintf(fp,"%20.16f %20.16f\n",t,0.);
        foreach()
        {
            fprintf (fp, "%20.16f %24.16e\n", x, q[]);
        }
        fclose (fp);
    }
}

/* Threshold > 1e-3 */
#if TREE
event adapt (i++) {
  adapt_wavelet ({q}, (double []){1e-4}, maxlevel = MAXLEVEL);
}
#endif









