#include "grid/quadtree.h"

#define USE_WAVEPROP 1

#if USE_WAVEPROP
#include "wpa_transport.h"
#else
#include "advection.h"
#endif

scalar f[];
scalar *tracers = {f};

int main() {
    origin (-0.5, -0.5);

#if USE_WAVEPROP    
    limiters[0] = 6;  /* 6 = generalized minmod (Basilisk default) */
    order = 2;
#else    
    gradient = minmod;
#endif    

  CFL = 0.8;
  for (N = 64; N <= 256; N *= 2)
    run();
}


#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

event init (i = 0) {
#if USE_WAVEPROP    
    dt_initial = 1e-5;  /* Probably too small */
#endif    
    /* Set this here so it doesn't get overwritten by defaults */
    foreach() {
        f[] = bump(x,y);    
    }
    boundary((scalar*) tracers);
}


#if USE_WAVEPROP
event velocity(i++) {
    trash ({u,v});
    /* Define a cell-centered velocity field */
    foreach() {
        u[] =  1.5*pi*sin(2.*pi*(t)/5.)*sin(pi*(x + 0.5))*cos(pi*(y + 0.5))/pi;  
        v[] = -1.5*pi*sin(2.*pi*(t)/5.)*cos(pi*(x + 0.5))*sin(pi*(y + 0.5))/pi;  
    }
    boundary ((scalar *){u, v});    
}

#else

event velocity (i++) {
    /* Define edge velocities using a stream-function*/
    vertex scalar psi[];
    foreach_vertex() {
        psi[] = - 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;        
    }
    trash ({u});
    struct { double x, y; } f = {-1.,1.};
    foreach_face() {
        u.x[] = f.x*(psi[0,1] - psi[])/Delta;        
    }
    boundary ((scalar *){u});
}
#endif


event logfile (t = {0,5})  {
    stats s = statsf (f);
    fprintf (stderr, "# %8.4f %24.16e %24.16e %24.16e\n", t, s.sum, s.min, s.max);
}

event field (t = 5) {
    scalar e[];
    foreach() {
        e[] = f[] - bump(x,y);       
    }
    norm n = normf (e);
    fprintf (stderr, "%-5d %24.16e %24.16e %24.16e\n", N, n.avg, n.rms, n.max);
    if (N == 256)
        output_field ({e}, stdout, N);      
}

event movie (t += 0.1) {
    if (N == 256) {
        output_ppm (f, linear = true, file = "f.mp4", n = N, min = 0, max = 1);
    }
}


