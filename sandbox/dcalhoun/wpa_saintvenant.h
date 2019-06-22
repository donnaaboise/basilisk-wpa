/** 

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']],
                    processEscapes: true},
            menuSettings: {
                    zoom: "Click",
               mpContext: true,
         mpMouse: true},
         errorSettings: {
             message: ["[Math Error]"] }
         });
</script>

<script type="text/javascript"
    src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# Saint Venant solver 

This script contains the Riemann solver for St. Venant equations.  The Riemann solver 
returns waves, speeds and fluctuations needed for the wave propagation algorithm.  The
equations are given by

$$
h_t + (u h)_x + (vh)_y  =  0 \\
(hu)_t + (hu^2 + \frac{1}{2}gh^2)_x  + (huv)_y   =  -ghb_x \\
(hv)_t + (hvu)_x + (hv^2 + \frac{1}{2}gh^2)_y    =  -ghb_y \\
$$

*/

#include "waveprop.h"

/**

The state variables are the conserved quantities $(h,hu,hv)$, which are set up 
as a list of scalars, 

*/
scalar h[];
vector hu[];
scalar * scalars = {h,hu};

/**
The governing equations are governed by a single parameter, the acceleration due to gravity. 
For this problem, we set this value to 1.   Users can override this value in their 
example code, e.g. [bump.c]().

*/

double grav = 1;

/**

The solution to the governing equations are updated via three waves in each direction 
- two gravity waves moving at speeds $u \pm \sqrt{gh}$ and one contact discontinuity moving at speed $u$. 

The number of waves expected by the Riemann solver, as well as the array of limiters 
used for each wave are given here. 
*/

#define MWAVES 3
int limiters[MWAVES];

event defaults(i=0) {

/**
The wave propagation algorithm in [waveprop.h]() requires several parameters 
which are specific to the problem being solved.  
In this Riemann solver, we are solving a conservation law, and will use the standard
solver (not the F-wave solver) and set we set
*/

    conservation_law = true;
    use_fwaves = false;

/**
For the "well-balanced" solver that handles bathymetry, we use the 
F-wave approach, described in the solver [wpa_stvenant_fwave.h]().
*/

/**
Set the number of waves used in this solver. 

*/

    mwaves = MWAVES;

/**

We can use use the option to set the same limiter for each wave.  In this case, we use
the generalized minmod limiter for all three waves, 
*/
    
    wave_limiter = 6;

/**
More generally, the user 
may wish to override the default behavior set above and specify limiters for each wave. 
In this case, they would set 3 entries of the array ${\tt limiters}$ (defined above)
to the desired wave limiters.  

In any case, the wave  propagation algorithm expects that ${\tt mthlim}$ is 
set and points to an integer array of length ${\tt mwaves}$.
*/

    mthlim = limiters;

/**

Set up the  list of state variables expected by the wave propagation algorithm.  This solver does not require any auxiliary variables.  On the other hand, if we were solving a problem with bathymetry, we would use the  auxiliary array to store topographic terms.
*/

    statevars = list_copy(scalars);
}

/**
This Riemann solver is called for each cell face.  The arguments are the left and  right states (stored in ${\tt ql}$ and ${\tt qr}$) and left and right auxiliary data  (stored in ${\tt auxl}$ and ${\tt auxr}$). An additional parameter ${\tt dir}$ is used to  indication the direction of the Riemann solve (x-direction or y-direction).  

The Riemann solver for the wave-propagation algorithm should then return waves, wave speeds and "fluctuations".  These will be used to compute second order corrections (if required) and update the solution.
*/


void wpa_rpsolver (int dir, int meqn, int mwaves, 
                   double *ql, double *qr, 
                   double *auxr, double *auxl,
                   double *waves, double *speeds, 
                   double *amdq, double *apdq, 
                   double *flux)
{
  if (meqn != 3) {
    fprintf (stderr, "Error : meqn != 3\n");
    exit(0);
  }

  int mu = (dir == 0) ? 1 : 2;
  int mv = (dir == 0) ? 2 : 1;
  
  double  hl = ql[0],  hr  = qr[0];
  double hul = ql[mu], hur = qr[mu];
  double hvl = ql[mv], hvr = qr[mv];
  
  double delta[3];  
  delta[0] = qr[0]   - ql[0];
  delta[1] = qr[mu]  - ql[mu];
  delta[2] = qr[mv]  - ql[mv];
  
  if (hl <= 0 || hr <= 0) {
    fprintf (stderr, "hr, hl <= 0;   %g %g\n", hl, hr);
    exit (0);
  }

/**
Roe averaged values $\overline{u}$, $\overline{v}$ and $\overline{c}$ are
computed and used in computing waves and speeds, below.
*/

  double ubar, vbar, cbar, divsqrt;
  divsqrt = sqrt(hl) + sqrt(hr);
  
  ubar = (hul/sqrt(hl) + hur/sqrt(hr))/divsqrt;
  vbar = (hvl/sqrt(hl) + hvr/sqrt(hr))/divsqrt;
  cbar = sqrt(grav*(hr + hl)/2.0);
  
  double a1, a2, a3, cbar2;
  cbar2 = 2*cbar;
  a1 = (-delta[1] + (ubar + cbar) * delta[0])/cbar2;
  a2 = -vbar*delta[0] + delta[2];
  a3 = (delta[1] - (ubar - cbar) * delta[0])/cbar2;
  

/**
### Waves and speeds

At each interface, we define waves and speeds.   To get the waves, we solve

$$
R \alpha = \delta, \qquad \mbox{where} \qquad \delta \equiv \mathbf q_r - \mathbf q_l
$$

where $R$ is a $3 \times 3$ matrix whose columns are the eigenvectors of the flux Jacobian matrix, evaluated the Roe averaged states (defined above).  

The waves are then given by

$$
\mathcal W^p = \alpha^p r^p, \qquad p = 1,2,3
$$

where $r^p$ is the $p^{th}$ eigenvector of the flux Jacobian, and the 
$p^{th}$ column of $R$. 

The speeds are the eigenvalues of the flux Jacobian matrix.
*/

/**

### Wave 1

*/
  int mw = 0;
  waves[0  + mw*meqn] = a1;
  waves[mu + mw*meqn] = a1*(ubar - cbar);
  waves[mv + mw*meqn] = a1*vbar;
  speeds[mw] = ubar - cbar;
  
/**

### Wave 2

*/

  mw = 1;
  waves[0  + mw*meqn] = 0;
  waves[mu + mw*meqn] = 0;
  waves[mv + mw*meqn] = a2;
  speeds[mw] = ubar;
  

/**

### Wave 3

*/

  mw = 2;
  waves[0  + mw*meqn] = a3;
  waves[mu + mw*meqn] = a3*(ubar + cbar);
  waves[mv + mw*meqn] = a3*vbar;
  speeds[mw] = ubar + cbar;
  

/**

### Fluctuations

Fluctuations $\mathcal A^+ \Delta Q$ and $\mathcal A^- \Delta Q$ 
are defined according to the speed of the f-wave at each interface.  
*/                

  for(int m = 0; m < meqn; m++) {
    amdq[m] = 0;
    apdq[m] = 0;
  }
  
  for(int m = 0; m < meqn; m++)
    for(int mw = 0; mw < mwaves; mw++)
      if (speeds[mw] < 0) 
    	amdq[m] += speeds[mw]*waves[m+mw*meqn];
      else
	    apdq[m] += speeds[mw]*waves[m+mw*meqn];  
  
  flux[0] = hur;
  flux[mu] = hur*hur/hr + 0.5*grav*hr*hr;  
  flux[mv] = hur*hvr/hr;                
}
