/** 
St. Venant using Wave Propagation Algorithm 
*/

#include "grid/quadtree.h"

#include "waveprop.h"
#include "rpn2_swe.h"

scalar h[];
vector hu[];
scalar *scalars = {h,hu};


int mwaves = 3;
int limiters[3] = {4,4,4}, *mthlim = limiters;    

/* Set these in 'init' */
bool dt_fixed;
double dt_initial;

bool conservation_law;
int order;


vector *aux = NULL;
int maux = 0;

int matlab_out = true;

#if TREE
#define ADAPT 1
#else
#define ADAPT 0
#endif


#define MAXLEVEL 8
#define MINLEVEL 5

static double sum0[3];
static double threshold = 4e-3;

int main() 
{
    origin(-0.5,-0.5);

    CFL = 0.9;

    periodic(right);
    periodic(top);

    N = 1 << MINLEVEL;
    run();
}

event init (i = 0) 
{
    wpa_rpsolver = rpn2_swe;

    dt_fixed = true;
    dt_initial = 2.5e-3;   /* to get CFL= 0.9 on first step */

    conservation_law = true;
    order = 2;

#if ADAPT
    astats s;
    int k = 0;
    do {
        foreach()       
        {
            h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
        }
        boundary (all);      
        s = adapt_wavelet ({h}, (double []){threshold}, maxlevel = MAXLEVEL);
        fprintf(stdout,"Initial conditions : k = %d\n",k);
        if (k > MAXLEVEL-MINLEVEL)
        {
            fprintf(stdout,"Initial conditions : Maximum refinement exceeded\n");
            break;
        }
        k++;
    } while (s.nf > 0);
#else
    foreach()       
    {
        h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
    }
#endif    
    boundary(scalars);
}

event logfile (i++)
{
  fprintf (stderr, "%g %g\n", t, normf(hu.x).max);
}

event movie (i += 10)
{
  output_ppm (h, linear = true, file = "h.mp4", n = 512);
  scalar l[];
  foreach()
      l[] = level;
  output_ppm (l, file = "level.mp4", min = MINLEVEL, max = MAXLEVEL, n = 512);
}

event conservation_check (t += 1)
{
    double sum[3] = {0,0,0};
    foreach()
    {
        int m = 0;
        for(scalar s in scalars)
        {
            sum[m] += Delta*Delta*s[];
            m++;
        }
    }

    /* Save initial mass for comparison */
    if (i == 0)
    {
        for(int m = 0; m < 3; m++)
            sum0[m] = sum[m];
    }

    fprintf(stdout,"Conservation check: Step %3d : t = %16.8e\n",i,t);
    int m = 0;
    for(scalar s in scalars)
    {
        fprintf(stdout,"sum[%d] = %24.16f %24.16e\n",m, sum[m], fabs(sum0[m]-sum[m]));
        m++;
    }
    fprintf(stdout,"\n");
}


event end (i = 1600)
{
  printf ("i = %d t = %g\n", i, t);
}


#if ADAPT
event adapt (i++) 
{
  adapt_wavelet ({h}, (double []){threshold}, maxlevel = MAXLEVEL);
}
#endif

/**
Here are the results.

![Animation of free surface height](swe/h.mp4)

![Animation of the level of refinement](swe/level.mp4)

~~~gnuplot Maximum velocity as a function of time
set xlabel 'Time'
set ylabel 'Maximum velocity'
unset key
plot 'log' w l
~~~
*/




