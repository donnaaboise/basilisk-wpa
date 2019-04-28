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

#define MINLEVEL 8
#define MAXLEVEL 10

static double threshold = 1e-4;


int main() 
{
    origin (-1.0,-1.0);
    L0 = 2;

    CFL = 0.9;

    periodic (left);        
    periodic (bottom);

    N = 1 << MINLEVEL;
    run();
}




event init (i = 0) 
{
    dt_initial = 1.13e-4;
    dt_fixed = false;
    conservation_law = false;
    order = 2;

    wpa_rpsolver = rpn2_adv;

#if 0
    vertex scalar phi[];   
    foreach_vertex()       
    {
        phi[] = 0.5 - sqrt(x*x + y*y);                
    }
    boundary ({phi});      
    fractions (phi, q);    
#else
    foreach() 
    {
        //q[] = exp(-100.0*sq(x*x + y*y));     
        //q[] = sin(4*pi*x);
        q[] = sqrt(x*x + y*y) <= 0.5;
    }
#endif    
    astats s = adapt_wavelet ({q}, (double []){threshold}, maxlevel = MAXLEVEL);
    while (s.nf > 0)
        s = adapt_wavelet ({q}, (double []){threshold}, maxlevel = MAXLEVEL);
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
#if TREE
event adapt (i++) 
{
  adapt_wavelet ({q}, (double []){threshold}, maxlevel = MAXLEVEL);
}
#endif










