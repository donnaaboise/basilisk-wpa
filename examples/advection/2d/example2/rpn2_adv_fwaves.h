/**
 Compute waves, speeds and fluctuations at each interface for the 
   scalar advection problem. 
*/

void rpn2_adv_fwaves(int dir, int meqn, int mwaves, 
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
