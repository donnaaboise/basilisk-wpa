/* Compute waves, speeds and fluctuations at each interface for the 
   scalar advection problem. */


static double ubar[2] = {1,0.5};  

void rpn2_adv(int dir, int meqn, int mwaves, double *ql, double *qr, double *waves, 
              double *speeds, double *amdq, double *apdq, double* flux)
{
    if (meqn != 1) 
    {
        printf("Error : meqn != 1\n");
        exit(0);
    }

    waves[0] = qr[0] - ql[0];  

    speeds[0] = ubar[dir];  /* Assume same direction in both x and y for now */

    amdq[0] = 0;
    apdq[0] = 0;
    for(int m = 0; m < meqn; m++)
        if (speeds[0] < 0) 
            amdq[0] = speeds[0]*waves[0];
        else
            apdq[0] += speeds[0]*waves[0];    

    flux[0] = 0;  /* Not used */
}
