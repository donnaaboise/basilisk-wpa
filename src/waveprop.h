/**
# Wave propagation Algorithm
*/

#include "utils.h"   /* Needed for DT, CFL, among other things */


/* ----------------------------- Defined by the user ---------------------------------- */
/* Lists of state variables */
extern scalar* scalars;   

/* Auxiliary info such as prescribed velocity fields, and so on */
extern scalar *aux;

/* Number of waves in system;  usually == number of equations, but not always */

/* Limiters to use for each wave. Should have length mwaves */
int mwaves;
int *mthlim;

/* Time stepping */
bool dt_fixed = false;      /* if true, use dt_initial for all time steps */
double dt_initial;

int order = 2;               /* 1 or 2 */
bool conservation_law = true;   /* true or false */
bool use_fwaves = false;

/* ------------------ Variables defined by WPA and referenced elsewhere ----------------- */
void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double *flux);

int Frame = 0; /* Used by plotting */

/* ----------------------- Static values used internally ------------------------------ */

scalar *statevars = NULL;
scalar *auxvars = NULL;


double dt;

/* ------------------------- Used by Basilisk and set here ---------------------------- */

/* Needed to get two layers of ghost cells at physical boundary.  Used here so we can do 
   limiting without storing wave fields first. */
#define BGHOSTS 2

/* ---------------------------------- Events  ----------------------------------------- */

event defaults (i = 0)
{
    CFL = 0.9;
    dt_initial = DT;

#if TREE
    for (scalar q in statevars) 
    {
        q.refine = q.prolongation = refine_linear;
        q.restriction = restriction_volume_average;
    }
#endif
}

event init (i = 0)
{
    /* User defined-init is called first */

    /* Call user defined values that should be set before anything else is checked */
#if 0    
    if (!(dt_initial < DT))
    {        
        printf("dt_initial is not set (event:user_parameters).\n");
        exit(0);
    }
#endif    
    dt = dt_initial;
    boundary (statevars);
    Frame = 0;
}

event plot(t >=0)
{
}

event cleanup(i=end, last)
{
    free(statevars);
}

/* ------------------------ Wave Propagation Algorithm -------------------------------- */

trace
void wpa_initialize(vector **wpa_fm, vector **wpa_fp, vector **wpa_flux)  
{
    *wpa_fm = NULL;
    *wpa_fp = NULL;
    *wpa_flux = NULL;

    int meqn = list_len(statevars);

    if (conservation_law)
        for(int m = 0; m < meqn; m++) {
            vector f = new face vector;
            *wpa_flux = vectors_append(*wpa_flux,f);            
        }
    else
        for(int m = 0; m < meqn; m++) {
            vector fmv = new face vector;
            *wpa_fm = vectors_append(*wpa_fm,fmv);

            vector fpv = new face vector;
            *wpa_fp = vectors_append(*wpa_fp,fpv);
        }
}

void wpa_cleanup(vector** wpa_fm, vector** wpa_fp, vector** wpa_flux)
{
    if (conservation_law) {
        delete((scalar*) *wpa_flux);
        free(*wpa_flux);   
    }
    else {
        delete((scalar*) *wpa_fp);
        free(*wpa_fp);

        delete((scalar*) *wpa_fm);
        free(*wpa_fm);        
    }                
}

attribute
{
    int limiter;
}

double wpa_limiter(double a, double b,int mlim)
{
    double wlimiter, r;
    r = b/a;

    switch (mlim)
    {
        double c,th;
        case 0:
            wlimiter = 1;
            break;

        /* Minmod */
        case 1:
            wlimiter = max(0., min(1., r));
            break;

        /* Superbee */
        case 2:
            wlimiter = max(0., max(min(1., 2.*r), min(2., r)));
            break;

        /* Van Leer */
        case 3:

            wlimiter = (r + fabs(r)) / (1. + fabs(r));
            break;

        /* Monotonized - Centered */
        case 4:
            c = (1. + r)/2.;
            wlimiter = max(0., min(c, min(2., 2.*r)));
            break;

        /* Sweby */
        case 5:
            wlimiter = max(0., max(min(r,1.5),min(1.5*r,1.)));
            break;

        /* Generalized minmod */
        case 6:
            th = 1.3;
            wlimiter = max(0., min(th*r,min((1+r)/2.,th)));
            break;

        default:
            printf("Invalid limiter method.\n");
            exit(0);

    }
    return wlimiter; 
}

double wpa_advance(double dt, double *cflmax)                   
{

    double dtnew = DT;
    double cflmax_local = 0;


    boundary(statevars);

    vector *wpa_fp, *wpa_fm, *wpa_flux;
    wpa_initialize(&wpa_fm,&wpa_fp, &wpa_flux);

    int meqn = list_len(statevars);
    int maux = list_len(auxvars); 

    int dir = 0;
    foreach_dimension() {
        /* Sweep over x, y, z dimensions using dim-split algorithm */
        foreach_face(x, reduction(min:dtnew) reduction(max:cflmax_local)) {
            double amdq[meqn];
            double apdq[meqn];

            double speeds[mwaves];

            double flux[meqn];
            double waves[meqn*mwaves];    /* waves at interface I */

            double ql[meqn];
            double qr[meqn];

            double auxl[maux];
            double auxr[maux];

            int m = 0;        
            for (scalar q in statevars) {
                qr[m] = q[0];
                ql[m] = q[-1];
                m++;
            }

            m = 0;
            for (scalar a in auxvars) {
                auxr[m] = a[0];
                auxl[m] = a[-1];                    
                m++;
            }

            /* Check s.v.x.i == s.i   --> is x component */

            /* In case the user doesn't define a flux */
            for(int m = 0; m < meqn; m++)
                flux[m] = 0;

            wpa_rpn2(dir,meqn, mwaves, ql, qr, auxl, auxr, waves, speeds, amdq, apdq, flux); 

            /* Get new time step for next step */
            for(int mw = 0; mw < mwaves; mw++) {
                /* Compute next step */
                double s = fabs(speeds[mw]);
                //fprintf(stdout,"speeds[%d]= %24.16f\n",mw,speeds[mw]);

                if (s == 0)
                    continue;  /* speed of 0 does not constrain the time step */

                double dt_local = CFL*Delta/s;
                if (dt_local < dtnew)
                    dtnew = dt_local;  
                                  
                /* Test current CFL */
                double cfl = dt*s/Delta;
                if (cfl > cflmax_local)
                    cflmax_local = cfl;

#if 0
                if (cflmax_local > 1)
                {
                    fprintf(stdout,"CFL > 1 (stop); cflmax = %6.2f; dt = %g; s = %g\n",
                            *cflmax, dt, s);
                }
#endif                
            } 

            /* First order update */            
            if (conservation_law) {
                m = 0;
                for (vector f in wpa_flux) 
                {
                    f.x[] = flux[m] - apdq[m];   
                    m++;
                }
            }
            else {
                vector fm;
                vector fp;
                m = 0;
                for (fp, fm in wpa_fp, wpa_fm) 
                {
                    fm.x[] = amdq[m];                
                    fp.x[] = -apdq[m]; 
                    m++;
                }                
            }

            /* --------------------- Second order corrections ----------------------------- */
            if (order == 2) {
                /* None of these are saved */
                double wavesl[meqn*mwaves];   
                double wavesr[meqn*mwaves];               
                double s[mwaves]; 
                double ap[meqn];
                double am[meqn];
                double flx[meqn];

                double qm2[meqn];
                double qm1[meqn];
                double q0[meqn];
                double qp1[meqn];

                double auxm2[maux];
                double auxm1[maux];
                double aux0[maux];
                double auxp1[maux];

                /* Right and left waves */
                int m = 0;
                for (scalar q in statevars) {
                    qm2[m] = q[-2];
                    qm1[m] = q[-1];
                    q0[m]  = q[0];
                    qp1[m] = q[1];

                    m++;
                }

                m = 0;
                for (scalar a in auxvars) {
                    auxm2[m] = a[-2];
                    auxm1[m] = a[-1];
                    aux0[m]  = a[0];
                    auxp1[m] = a[1];
                    m++;                    
                }

                bool use_limiting = false;
                for(int mw = 0; mw < mwaves; mw++)
                    use_limiting = use_limiting || mthlim[mw] > 0;

                double wlimiter[mwaves];
                if (!use_limiting) {
                    for(int mw = 0; mw < mwaves; mw++)
                        wlimiter[mw] = 1;                                    
                }
                else {
                    /* Get left and right waves so we can do limiting in place */
                    wpa_rpn2(dir,meqn, mwaves, qm2, qm1, auxm2, auxm1, wavesl, s, am, ap,flx);                           
                    wpa_rpn2(dir,meqn, mwaves, q0,  qp1, aux0, auxp1, wavesr, s, am, ap,flx);                
                    /* -------------- Compute limited second order corrections ---------------- */
                    for(int mw = 0; mw < mwaves; mw++) {
                        double w2 = 0;
                        double dr = 0;  
                        double dl = 0;
                        for (int m = 0; m < meqn; m++) {
                            double w = waves[m + mw*meqn];
                            double wr = wavesr[m + mw*meqn];
                            double wl = wavesl[m + mw*meqn];
                            w2 += w*w;
                            dr += w*wr;
                            dl += w*wl;
                        }


                        if (w2 == 0 || mthlim[mw] == 0)
                            wlimiter[mw] = 1;
                        else {                
                            if (speeds[mw] >= 0)
                                wlimiter[mw] = wpa_limiter(w2,dl,mthlim[mw]);
                            else if (speeds[mw] < 0)
                                wlimiter[mw] = wpa_limiter(w2,dr,mthlim[mw]);
                        }
                    }
                }

                /* Use dt passed in as argument;  dx is mesh width in cell x[0] */
                double dtdx = dt/Delta;
                double cqxx[meqn];
                for(int m = 0; m < meqn; m++)
                    cqxx[m] = 0;

                for(int mw = 0; mw < mwaves; mw++) {
                    double abs_sign;
                    if (use_fwaves)
                        abs_sign = sign(speeds[mw]);
                    else
                        abs_sign = fabs(speeds[mw]);

                    double cq = abs_sign*(1 - fabs(speeds[mw])*dtdx)*wlimiter[mw];
                    for(int m = 0; m < meqn; m++)
                        cqxx[m] += cq*waves[m+mw*meqn];
                }

                if (conservation_law) {
                    m = 0;
                    for (vector f in wpa_flux) {
                        f.x[]  += 0.5*cqxx[m];  
                        m++;
                    }  
                }
                else {
                    vector fm;
                    vector fp;
                    m = 0;
                    for (fp,fm in wpa_fp, wpa_fm) {
                        fm.x[] += 0.5*cqxx[m];
                        fp.x[] += 0.5*cqxx[m];   
                        m++;
                    }  
                }
            }   /* End of second order corrections */
        }

        /* Replace coarse grid fluxes with fine grid fluxes */
        if (conservation_law)
            boundary_flux(wpa_flux);
        else {
            boundary_flux(wpa_fp);
            boundary_flux(wpa_fm);            
        }

        /* ------------------------------ Update solution --------------------------------- */
        foreach() {
            double dtdx = dt/Delta;

            if (conservation_law) {            
                vector f;
                scalar q;
                for (f,q in wpa_flux, statevars) 
                    q[] += -dtdx*(f.x[1] - f.x[]);
            }
            else {
                vector fp;
                vector fm;
                scalar q;
                for (fp, fm, q in wpa_fp, wpa_fm, statevars) 
                    q[] += -dtdx*(fm.x[1] - fp.x[]);  
            }
        }
        dir++;
    } /* End of dimension loop */

    wpa_cleanup(&wpa_fm, &wpa_fp, &wpa_flux);

    *cflmax = cflmax_local;
    return dtnew;
}

void run()
{
    t = 0.;
    iter = 0;
    init_grid (N);  

    double cflmax = 0;
    double dtnew;
    //double dt_curr;

    /* Events are processed first, followed by statements in the while loop */
    perf.nc = perf.tnc = 0;
    perf.gt = timer_start();
    while (events (true)) {
        dtnew = wpa_advance(dt, &cflmax);
        if (cflmax > 1) 
        {
            /* CFL using time step dt */
            printf("CFL exceeds 1; cflmax = %g\n",cflmax);
            //exit(0);
        }    
#if 0        
        dt_curr = dt;
        t += dt_curr; /* Use step just taken */        

        if (cflmax > 1) 
        {
            /* CFL using time step dt */
            printf("CFL exceeds 1; cflmax = %g\n",cflmax);
            exit(0);
        }    
        if (dt_fixed)
        {
            dt = dt_initial;  
        }    
        else
        {
            dt = dtnext(dtnew);            
            //t = tnext;  This is t += dtnew;  we want t += dt (above). 
        }
#endif
        dt = dtnext(dtnew);
        t = tnext;
        iter = inext;
        fprintf(stdout,"Step %6d : dt = %8.4e; cflmax = %12.6f; Time = %12.4f\n",
                iter, dt, cflmax,t);
        fflush(stdout);
    }
    timer_print (perf.gt, iter, perf.tnc);

    free_grid();
}





