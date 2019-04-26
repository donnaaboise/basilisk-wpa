#include "utils.h"   /* Needed for DT, CFL, among other things */

/* ----------------------------- Defined by the user ---------------------------------- */
extern scalar* scalars;
extern vector* vectors;
scalar *statevars = NULL;

/* Number of waves in system;  usually == number of equations, but not always */
extern int mwaves;

/* Limiters to use for each wave. Should have length mwaves */
extern int *mthlim;

/* Time stepping */
extern bool dt_fixed;      /* if true, use dt_initial for all time steps */
extern double dt_initial;

/* Plotting */
extern int matlab_out;

/* Not sure what to do with these yet */
extern int maux;
extern vector *aux;

/* -------------------------- DEFS defined by WPA and used here ----------------------- */

/* These here are a hack and should be set so the user can change them */
#define SECOND_ORDER 1
#define CONSERVATION_LAW 0

/* ------------------ Variables defined by WPA and referenced elsewhere ----------------- */
int Frame = 0; /* Used by plotting */

typedef void (*wpa_rp1_t)(int meqn, int mwaves, 
                           double *ql, double *qr, 
                           double *waves, double *speeds, 
                           double *amdq, double *apdq, double *flux);

wpa_rp1_t wpa_rp1;   /* riemann solver */   


/* ----------------------- Static values used internally ------------------------------ */

//static int meqn;


/* ------------------------- Used by Basilisk and set here ---------------------------- */
/* Needed to get two layers of ghost cells at physical boundary */
#define BGHOSTS 2

static double dt;

/* ---------------------------------- Events  ----------------------------------------- */

event defaults (i = 0)
{
    dt_initial = DT;
#if TREE
    for (scalar q in statevars) 
    {
        q.refine = q.prolongation = refine_linear;
        q.restriction = restriction_volume_average;
    }
#endif
}

event setaux(i++)
{
    /* Set auxiliary variables such as velocity fields */
}

event init (i = 0)
{
    /* User defined-init is called first */
    if (!(dt_initial < DT))
    {        
        printf("dt_initial is not set.\n");
        exit(0);
    }

    dt = dt_initial;

    boundary (statevars);
}


event plot (i >= 0)
{
    if (matlab_out)
    {        
        char name[80];

        /* Write header file */
        sprintf(name,"fort.t%04d",Frame);
        FILE *fp = fopen(name,"w");
        fprintf(fp,"%20.16f %20s\n",t,"time");

        int meqn = list_len(statevars);
        fprintf(fp,"%10d %20s\n",meqn,"meqn");
        fprintf(fp,"%10d %20s\n",1,"ngrids");
        fprintf(fp,"%10d %20s\n",maux,"maux");
        fprintf(fp,"%10d %20s\n",dimension,"dim");
        fclose(fp);

        /* Write data file */
        sprintf(name,"fort.q%04d",Frame);
        fp = fopen(name,"w");
        fprintf(fp,"%10d %20s\n",1,"grid_number");
        fprintf(fp,"%10d %20s\n",1,"AMR_level");
        fprintf(fp,"%10d %20s\n",N,"mx");
        fprintf(fp,"%24.16f %20s\n",X0,"xlow");  /* xlow */

        /* Assume uniform Cartesian grid */
        double dx = L0/N;
        fprintf(fp,"%24.16f %20s\n",dx,"dx");
        fprintf(fp,"\n");
        foreach()
        {
            for(scalar s in statevars)
                fprintf(fp,"%24.16e",s[]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    printf("Output Frame %d at time t = %12.4e\n",Frame, t);
    printf("\n");
    Frame++;
}


event cleanup(i=end, last)
{
    free(statevars);
}

/* ------------------------ Wave Propagation Algorithm -------------------------------- */

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
            th = 2;
            wlimiter = max(0., min(th*r,min((1+r)/2.,th)));
            break;

        default:
            printf("Invalid limiter method.\n");
            exit(0);

    }
    return wlimiter; 
}


void wpa_initialize(vector **wpa_fm, vector **wpa_fp, vector **wpa_flux)  
{
    statevars = list_concat (scalars, (scalar *) vectors);       

    int meqn = list_len(statevars);

    *wpa_fm = NULL;
    *wpa_fp = NULL;
    *wpa_flux = NULL;
    for(int m = 0; m < meqn; m++)
    {
        vector fmv = new face vector;
        *wpa_fm = vectors_append(*wpa_fm,fmv);

        vector fpv = new face vector;
        *wpa_fp = vectors_append(*wpa_fp,fpv);

        vector f = new face vector;
        *wpa_flux = vectors_append(*wpa_flux,f);
    }
}

double wpa_advance(double dt, double* cflmax, vector* wpa_fm,
                    vector* wpa_fp, vector* wpa_flux)
{
    int meqn = list_len(statevars);

    double flux[meqn];
    double waves[meqn*mwaves];    /* waves at interface I */
    double amdq[meqn];
    double apdq[meqn];
    double speeds[mwaves];

    double ql[meqn];
    double qr[meqn];

    double dtnew = DT;
    *cflmax = 0;

    boundary(statevars);

    foreach_face()
    {
        int m;   

        /* ------------------ CENTERED waves ------------------ */
        /* These waves are used for the first order update */
        m = 0;        
        for (scalar q in scalars) 
        {
            qr[m] = q[0];
            ql[m] = q[-1];
            m++;
        }

        for (vector w in vectors) 
        {
            qr[m] = w.x[0];
            ql[m] = w.x[-1];
            m++;
        }

        wpa_rp1(meqn, mwaves, ql, qr, waves, speeds, amdq, apdq, flux); 

        /* Get new time step for next step */
        for(int mw = 0; mw < mwaves; mw++)
        {
            /* Compute next step */
            double dt_local = CFL*Delta/fabs(speeds[mw]);
            if (dt_local < dtnew)
                dtnew = dt_local;

            /* Test current CFL */
            double cfl = dt*fabs(speeds[mw])/Delta;
            if (cfl > *cflmax)
                *cflmax = cfl;
        } 

        /* First order update */
        vector fm;
        vector fp;
        vector f;

        m = 0;
        for (fm,fp,f in wpa_fm,wpa_fp,wpa_flux) 
        {
            fm.x[] = amdq[m];                
            fp.x[] = -apdq[m]; 

            f.x[] = flux[m] + fp.x[];   
            m++;
        }

/* --------------------------- SECOND ORDER CORRECTIONS ------------------------------- */
#if SECOND_ORDER            
    /* Compare waves at face I to waves at face I+1 and I-1.  Store norm of 
       wave(I), and projection of wave(I+1) onto wave(I).  These stored values
       will be used below in wave limiting 

       We take advantage of the fact that the Basilisk stencil has width 5 (2 ghost
       cells in each direction). 
       */
        
        {
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


            /* -------------------------- RIGHT AND LEFT waves -------------------------*/
            m = 0;
            for (scalar q in scalars) 
            {
                qm2[m] = q[-2];
                qm1[m] = q[-1];
                q0[m]  = q[0];
                qp1[m] = q[1];
                m++;
            }
            for (vector w in vectors) 
            {
                qm2[m] = w.x[-2];
                qm1[m] = w.x[-1];
                q0[m]  = w.x[0];
                qp1[m] = w.x[1];
                m++;
            }

            /* Get left and right waves so we can do limiting "in place" */
            wpa_rp1(meqn, mwaves, qm2, qm1, wavesl, s, am, ap,flx);                            
            wpa_rp1(meqn, mwaves, q0,  qp1, wavesr, s, am, ap,flx);                


            /* -------------- COMPUTE LIMITED SECOND ORDER CORRECTIONS ---------------- */
            double wlimiter[mwaves];
            for(int mw = 0; mw < mwaves; mw++) 
            {
                double w2 = 0;
                double dr = 0;  /* These are being computed 2-3 times each (:-(( */
                double dl = 0;
                for (int m = 0; m < meqn; m++)
                {
                    double w = waves[m + mw*meqn];
                    double wr = wavesr[m + mw*meqn];
                    double wl = wavesl[m + mw*meqn];
                    w2  += w*w;
                    dr += w*wr;
                    dl += w*wl;
                }

                if (w2 == 0 || mthlim[mw] == 0)
                    wlimiter[mw] = 1;
                else
                {                
                    if (speeds[mw] >= 0)
                        wlimiter[mw] = wpa_limiter(w2,dl,mthlim[mw]);
                    else if (speeds[mw] < 0)
                        wlimiter[mw] = wpa_limiter(w2,dr,mthlim[mw]);
                }
            }

            /* Use dt passed in as argument (set by user at t=0) */
            double dtdx = dt/Delta;
            double cqxx[meqn];
            for(int m = 0; m < meqn; m++)
                cqxx[m] = 0;

            for(int mw = 0; mw < mwaves; mw++)
            {
                double cq = fabs(speeds[mw])*(1 - fabs(speeds[mw])*dtdx)*wlimiter[mw];
                for(int m = 0; m < meqn; m++)
                {
                    cqxx[m] += cq*waves[m+mw*meqn];
                }
            }

            m = 0;
            vector f;
            for (fm,fp,f in wpa_fm,wpa_fp,wpa_flux) 
            {
                fm.x[] += 0.5*cqxx[m];
                fp.x[] += 0.5*cqxx[m];   
                f.x[]  += 0.5*cqxx[m];  /* Note sign */
                m++;
            }  
        }   /* End of second order scope */
#endif   
/* ----------------------- END OF SECOND ORDER CORRECTIONS ---------------------------- */

    }  /* End of foreach () */  

    boundary_flux(wpa_fp);
    boundary_flux(wpa_fm);
    boundary_flux(wpa_flux);

    if (dt_fixed)
    {
        dtnew = dt_initial;
    }

    /* Update solution */

    foreach()
    {
        double dtdx = dt/Delta;
        scalar q;

#if CONSERVATION_LAW
        vector f;
        for (f,q in wpa_flux, statevars) 
        {
            q[] += -dtdx*(f.x[1] - f.x[]);
        }
#else    
        vector fp;
        vector fm;
        for (fp, fm, q in wpa_fp, wpa_fm, statevars) 
        {
            q[] += -dtdx*(fm.x[1] - fp.x[]);  
        }
#endif        
    }
    return dtnew;
}

void wpa_cleanup(vector **wpa_fm, vector **wpa_fp, vector **wpa_flux)  
{
    free(*wpa_fp);
    free(*wpa_fm);
    free(*wpa_flux);    
}


void run()
{
    t = 0.;
    iter = 0;
    init_grid (N);  
    double cflmax;
    double dtnew;

    /* Set up global variables */
    vector *wpa_fp;
    vector *wpa_fm;
    vector *wpa_flux;

    wpa_initialize(&wpa_fm, &wpa_fp, &wpa_flux);

    /* Events are processed first, followed by statements in the while loop */
    perf.nc = perf.tnc = 0;
    perf.gt = timer_start();
    while (events (true)) 
    {
        dtnew = dtnext(wpa_advance(dt, &cflmax, wpa_fm, wpa_fp, wpa_flux));
        if (cflmax > 1) 
        {
            printf("CLF exceeds 1; cflmax = %g\n",cflmax);
        }        
        //t += dt;
        //dt = dtnext(dtnew);  /* Updates tnext with new dt. */


        iter = inext;
        t = tnext;    
    }
    timer_print (perf.gt, iter, perf.tnc);

    wpa_cleanup(&wpa_fm, &wpa_fp, &wpa_flux);

    free_grid();
}





