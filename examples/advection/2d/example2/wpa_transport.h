/**
 Compute waves, speeds and fluctuations at each interface for the 
   scalar advection problem. 
*/

#include "waveprop.h"

scalar u[];
scalar v[];

#define MWAVES 1

extern scalar* tracers;
int limiters[MWAVES], *mthlim = limiters;  

/* Defaults for this Riemann solver */
event defaults(i=0)
{
    /* Set up variables recognized by waveprop algorithm, but specific
       to this application. */
    conservation_law = true;
    use_fwaves = true;

    mwaves = MWAVES;
    /* don't set defaults here because they will get overwritten by definitions 
       in main */
    // limiters[0] = 6;  

    statevars = list_copy(tracers);  
    auxvars = list_concat({u},{v});
}

void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double* flux)
{
    if (mwaves != 1) 
    {
        printf("Error  (wpa_transport.h : mwaves != 1\n");
        exit(0);
    }

    double ur = auxr[dir];  
    double ul = auxl[dir];

    double uhat = (ur + ul)/2.;

    for(int m = 0; m < meqn; m++)
    {
        waves[m] = ur*qr[m] - ul*ql[m];            
        if (uhat >= 0)
        {
            amdq[m] = 0;
            apdq[m] = waves[m];    
        }
        else
        {
            amdq[m] = waves[m];
            apdq[m] = 0;    
        }        
        flux[m] = ur*qr[m]; 
    }
    speeds[0] = uhat;
}
