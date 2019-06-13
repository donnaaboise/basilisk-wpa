/**
 This solver solves the conservative transport problem given by

 $$
  \partial_t \mathbf q + \nabla(\mathbf u \cdot \mathbf q) = -
$$
where $\mathbf u = (u(x,y),v(x,y))$ is a prescribed velocity field and 
$\mathbf q(x,y,t)$ is a vector of tracer quantities. 

The solver is based on the wave propagation algorithms, first described by 
R. J. LeVeque [LeVeque, 2002](/src/refs.bib#le:2002)
*/

#include "waveprop.h"

scalar u[];
scalar v[];

extern scalar* tracers;

#define MWAVES 1
int limiters[MWAVES], *mthlim = limiters;  

/* Defaults for this Riemann solver */
event defaults(i=0)
{
    conservation_law = true;
    use_fwaves = true;

    mwaves = MWAVES;
    /* The user should set definitions of limiters in main */
    // limiters[0] = 6;  /* This would overwrite user defs in main */

    statevars = list_copy(tracers);
    auxvars = list_concat({u},{v});
}

void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double* flux) 
{
    if (mwaves != 1) {
        fprintf(stdout,"Error  (wpa_transport.h : mwaves != 1\n");
        exit(0);
    }

    double ur = auxr[dir];  
    double ul = auxl[dir];

    double uhat = (ur + ul)/2.;

    for(int m = 0; m < meqn; m++) {
        waves[m] = ur*qr[m] - ul*ql[m];            
        if (uhat >= 0) {
            amdq[m] = 0;
            apdq[m] = waves[m];    
        }
        else {
            amdq[m] = waves[m];
            apdq[m] = 0;    
        }        
        flux[m] = ur*qr[m]; 
    }
    speeds[0] = uhat;
}
