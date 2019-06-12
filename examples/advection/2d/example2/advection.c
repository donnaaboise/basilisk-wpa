#include "grid/multigrid.h"

#define USE_WAVEPROP 1

#if USE_WAVEPROP
#include "wpa_transport.h"
#else
//#include "advection.h"

/* Custom version that fixes the time step */
#include "advection_basilisk.h"
#endif



#if 1
/* Single tracer */
scalar f[];
scalar *tracers = {f};
#else
/* Test multiple tracers */
scalar f1[];
scalar f2[];
scalar *tracers = {f, f1, f2};
#endif



/* Don't change these */
#define DT0     2e-2
#define N0      16
#define NOUT0   250    /* nout for N=16; dt = 2e-2; T_final = 5 */
#define NSTEP0  25

/* Change FACTOR to increase resolution */
#define FACTOR  16    /* 1,2,4,8,16,32,64, ... */


/* Computed from FACTOR.  Note : Don't re-define N here! */
#define NOUT    FACTOR*NOUT0
#define NSTEP   FACTOR*NSTEP0
#define TFINAL  DT0/FACTOR*NOUT

int main()
{
    // coordinates of lower-left corner
    origin (-0.5, -0.5);

    N = FACTOR*N0;

#if USE_WAVEPROP    
    limiters[0] = 0;  /* Number corresponds to limiter type;  0 = no limiter */
    dt_fixed = true;
    order = 2;
#else    
    gradient = NULL;
#endif    

    run();
}


#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))
#define Hsmooth(r) ((tanh((r)/0.015625) + 1)/2.0)

static
double qinit(double x, double y)
{
    double x0 = -0.2;
    double y0 = -.236338;
    double r2 = sq(x-x0) + sq(y-y0);
    //return exp(-100*r2);

    double r = sqrt(r2);
    double r0 = 0.125;
    return Hsmooth(r+r0) - Hsmooth(r-r0);
}

event init (i = 0)
{
#if USE_WAVEPROP
    dt_initial = DT0/FACTOR;
#endif    
    foreach()
    {
        for (scalar q in tracers) 
        {
            q[] = qinit(x,y);
        }
    }

    boundary((scalar*) tracers);
}


#if USE_WAVEPROP
event velocity(i++)
{
    trash ({u,v});
    if (!use_fwaves)

    {
        vertex scalar psi[];
        foreach_vertex()
        {
            psi[] = -1.5*sin(2.*pi*(t)/5.)*sin(pi*(x + 0.5))*sin(pi*(y + 0.5))/pi;        
        }
        foreach()
        {
            u[] = -(psi[0,1] - psi[])/Delta;      /* u = -dpsi/dy */
            v[] = (psi[1,0] - psi[])/Delta;       /* v = dpsi/dx  */
        }
    }
    else
    {
        foreach()
        {
            u[] =  1.5*pi*sin(2.*pi*(t)/5.)*sin(pi*(x + 0.5))*cos(pi*(y + 0.5))/pi;  
            v[] = -1.5*pi*sin(2.*pi*(t)/5.)*cos(pi*(x + 0.5))*sin(pi*(y + 0.5))/pi;  
        }
    }
    boundary ((scalar *){u, v});    
}
#else

event velocity (i++) 
{
    vertex scalar psi[];
    foreach_vertex()
    {
        psi[] = - 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;        
    }
    trash ({u});
    struct { double x, y; } f = {-1.,1.};
    foreach_face()
    {
        u.x[] = f.x*(psi[0,1] - psi[])/Delta;        
    }
    boundary ((scalar *){u});
    dt = DT0/FACTOR;
}
#endif

event logfile (i = {0,NOUT}) 
{
    stats s = statsf (f);
    fprintf (stderr, "# %8.4f %24.16e %24.16e %24.16e\n", t, s.sum, s.min, s.max);
}

event verbosity (i += NSTEP)
{
    fprintf(stdout,"Step %6d : dt = %8.4e; Time = %12.4f\n",i, dt, t);
    fflush(stdout);
}

event field (i = NOUT) 
{
    scalar e[];
    foreach()
    {
        e[] = f[] - qinit(x,y);       /* Choose 1 of f1, f2 or f3 */ 
    }
    norm n = normf (e);
    fprintf (stderr, "%-5d %24.16e %24.16e %24.16e\n", N, n.avg, n.rms, n.max);
  
#if 0
    if (N == 256)
        output_field ({e}, stdout, N);
#endif        
}

#if 1
event movie (i += 10)
{
  output_ppm (f, linear = true, file = "f.mp4", n = N, min = 0, max = 1);
#if TREE
  scalar l[];
  foreach()
      l[] = level;
  output_ppm (l, file = "level.mp4", min = 0, max = 1, n = N);
#endif
}
#endif

event end (i = NOUT)
{
  printf ("i = %d t = %g\n", i, t);
}

