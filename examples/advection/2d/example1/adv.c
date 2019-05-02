#include "grid/quadtree.h"  

#include "fractions.h"
#include "waveprop.h"
#include "rpn2_adv.h"

scalar q[];
scalar *scalars = {q};
vector *vectors = NULL;


bool dt_fixed;
double dt_initial;

bool conservation_law;
int order;

int mwaves = 1;

int limiters[1] = {3}, *mthlim = limiters;

int maux;
vector *aux;

int matlab_out = true;

#define MINLEVEL 5
#define MAXLEVEL 9
#define ADAPT 1

static double threshold = 1e-4;
static double sum0;


int main() 
{
    origin (-1.0,-1.0);
    L0 = 2;

    CFL = 0.9;

    periodic (left);        
    periodic (bottom);

    N = 1 << MAXLEVEL;
    run();
}




event init (i = 0) 
{
    dt_initial = 1.13e-4;
    dt_fixed = false;
    conservation_law = false;
    order = 2;

    wpa_rpsolver = rpn2_adv;

    astats s;
    do {
        vertex scalar phi[];   
        foreach_vertex()       
        {
            phi[] = 0.5 - sqrt(x*x + y*y);                
        }
        boundary ({phi});      
        fractions (phi, q);    
        s = adapt_wavelet ({q}, (double []){1e-5}, maxlevel = MAXLEVEL);
    } while (s.nf);

#if 0
    //astats s = adapt_wavelet ({q}, (double []){threshold}, maxlevel = MINLEVEL);
    while (s.nf > 0)
        s = adapt_wavelet ({q}, (double []){1e-5}, maxlevel = MINLEVEL);
#endif        
}

event conservation_check(t=0; )
{
    double sum = 0;
    foreach()
    {
        sum += Delta*Delta*q[];
    }
    if (i == 0)
        sum0 = sum;

    printf("t = %20.16f; mass = %24.16e; diff = %24.16e\n",t, sum, fabs(sum0-sum));
    fflush(stdout);
}


event images (t += 0.25; t <= 4)
{
    scalar l[];
    foreach()
        l[] = level;
    static FILE * fp = fopen ("grid.ppm", "w");
    output_ppm (l, fp, min = 0, max = MAXLEVEL);
}


event plot (t += 0.25; t <= 4)
{
#if 0    
  static FILE * fp = fopen ("out.ppm", "w");
  output_ppm (q, fp, min = 0, max = MAXLEVEL);
#endif  
}

/* Threshold > 1e-3 */

#if ADAPT //TREE
event adapt (i++) 
{
  adapt_wavelet ({q}, (double []){threshold}, maxlevel = MAXLEVEL);
}
#endif










