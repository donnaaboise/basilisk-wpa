/**
 Compute waves, speeds and fluctuations at each interface for the 
   scalar advection problem. 
*/

#include "waveprop.h"

scalar u[];
scalar v[];
scalar *aux = {u,v};  /* So we can create a list of scalars */


/* Defaults for this Riemann solver */
event defaults(i=0)
{
    conservation_law = true;
    use_fwaves = true;
    mwaves = 1;
}

void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double* flux)
{
    if (meqn != 1) 
    {
        printf("Error  (rpn2_adv_fwaves : meqn != 1\n");
        exit(0);
    }

    double ur = auxr[dir];  
    double ul = auxl[dir];

    double uhat = (ur + ul)/2.;

    waves[0] = ur*qr[0] - ul*ql[0];
    if (uhat >= 0)
    {
        amdq[0] = 0;
        apdq[0] = waves[0];    
    }
    else
    {
        amdq[0] = waves[0];
        apdq[0] = 0;            
    }
    speeds[0] = uhat;

    flux[0] = ur*qr[0]; 
}
