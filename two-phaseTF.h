/**
# Two-phase interfacial flows

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in drop is $f1=1$ and $f2=0$. In the thin film, it is $f2=1$ and $f1=0$. Air (fluid 3) is $f1 = f2 = 0$. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*, *mu1*, *rho4*, *mu2*, respectively.
**Note:** The drop and the film are defined by different VoF fields, but have same properties (density and viscosity).
*/

#include "vof.h"
/**
Instead of one VoF tracer, we define three, f1, f2, and f3.
*/
scalar f1[], f2[], f3[], *interfaces = {f1, f2, f3};
double rho1 = 1., mu1 = 0., rho4 = 1., mu2 = 0., mu3 = 0., mu4 = 0.;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */
face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  if (mu1 || mu2 || mu3)
    mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). The difference comes in how we call these averages.
$$
\hat{A} = (f_1+f_2) + (1-f_1-f_2)\frac{A_g}{A_l}\,\,\,\forall\,\,\,A \in \{\mu,\rho\}
$$
*/

#ifndef rho
#define rho(f) (clamp(f,0.,1.)*(rho1 - rho4) + rho4)
#endif
#ifndef mu
#define mu(f1, f2, f3)  (clamp(f1,0.,1.)*mu1 + clamp(f2,0.,1.)*mu2 + clamp(f3,0.,1.)*mu3 + clamp(1.-f1-f2-f3,0.,1.)*mu4)
#endif


/**
We have the option of using some "smearing" of the density/viscosity
jump. It is modified to take into account that there are two VoF tracers. */

#ifdef FILTERED
scalar sf1[], sf2[], sf3[], *smearInterfaces = {sf1, sf2, sf3};
#else
#define sf1 f1
#define sf2 f2
#define sf3 f3
scalar *smearInterfaces = {sf1, sf2, sf3};
#endif

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. Introduce for loops to ensure that smearing is done properly. */

  #ifdef FILTERED
    int counter1 = 0;
    for (scalar sf in smearInterfaces){
      counter1++;
      int counter2 = 0;
      for (scalar f in interfaces){
        counter2++;
        if (counter1 == counter2){
          // fprintf(ferr, "%s %s\n", sf.name, f.name);
        #if dimension <= 2
            foreach(){
              sf[] = (4.*f[] +
          	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
          	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
            }
        #else // dimension == 3
            foreach(){
              sf[] = (8.*f[] +
          	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
          	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
          		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
          		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
          	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
          	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
            }
        #endif
        }
      }
    }
    #endif
  #if TREE
    for (scalar sf in smearInterfaces){
      sf.prolongation = refine_bilinear;
      boundary ({sf});
    }
  #endif



  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    double ff3 = (sf3[] + sf3[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff1+ff2+ff3);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1, ff2, ff3);
  }
  foreach()
    rhov[] = cm[]*rho(sf1[]+sf2[]+sf3[]);

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    boundary ({sf});
  }
#endif
}
