/* Title: Bouncing Droplet on thin/deep films!
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

// 1 is drop

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseTF.h"
#include "tension.h"
#include "adapt_wavelet_limited.h"

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in vorticity
#define DissErr (1e-3)

// air-water
#define Rho21 (0.001)

// Calculations!
#define Xdist (1.04)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))
// domain
#define Ldomain 8                                // Dimension of the domain

// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = 0.0;
f2[left] = 1.0;

int MAXlevel;
double tmax, We, Ohd, Bo, Ohf, hf;
#define MINlevel 3                                            // maximum level
#define tsnap (0.01)

int main(int argc, char const *argv[]) {
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  if (argc < 8){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
    return 1;
  }

  MAXlevel = atoi(argv[1]);
  tmax = atof(argv[2]);
  We = atof(argv[3]); // We is 1 for 0.167 m/s <816*0.167^2*0.00075/0.017>
  Ohd = atof(argv[4]); // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  Bo = atof(argv[5]); // <816*10*0.00075^2/0.017 = 0.27>
  Ohf = atof(argv[6]);
  hf = atof(argv[7]);

  if (hf == 0){
    fprintf(ferr, "We have a problem. Wrong code. Change code or film height.\n");
    return 1;
  }
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %g, Bo %g, Ohf %3.2e, hf %3.2f\n", MAXlevel, tmax, We, Ohd, Bo, Ohf, hf);

  L0=Ldomain;
  X0=-hf; Y0=0.;
  init_grid (1 << (MINlevel));

  rho1 = 1.0; mu1 = Ohd; mu2 = Ohf; mu3 = Ohf;
  rho4 = Rho21; mu4 = 1e-5;

  f1.sigma = 1.0; f2.sigma = 1.0; f3.sigma = 1.0;

  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y) < 1.44) && (level < MAXlevel));
    fraction (f1, 1.0 - R2Drop(x,y));
    fraction (f2, -x);
    fraction (f3, sq(0.5) - R2Drop(x,y));   
    foreach () {
      if (f1[] + f3[] > 1.0) {
       f1[]=1. - f3[];
    }
    if (f1[]>0.05 || f3[]>0.05) u.x[]=-sqrt(We);
    else u.x[]=0.0;
      u.y[] = 0.0;
    }
    boundary((scalar *){f1, f2, f3, u.x, u.y});
  }
}

// Gravity
event acceleration(i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] -= Bo;
  }
}

int refRegion(double x, double y, double z){
  return (x > -hf && x < 0.1 && y < 2.0 ? MAXlevel+1 : MAXlevel);
}

event adapt(i++){
  scalar KAPPA1[], KAPPA2[], KAPPA3[], omega[];
  vorticity (u, omega);
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  curvature(f3, KAPPA3);
  scalar D2c[];
  foreach(){
    omega[] *= (f1[]+f2[]+f3[]);
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = (f1[]*Ohd+f2[]*Ohf+f3[]*Ohf)*D2;
  }
  boundary((scalar *){D2c, omega, KAPPA1, KAPPA2});
  adapt_wavelet_limited ((scalar *){f1, f2, f3, KAPPA1, KAPPA2, KAPPA3, u.x, u.y, omega, D2c},
     (double[]){fErr, fErr, fErr, KErr, KErr, KErr, VelErr, VelErr, OmegaErr, DissErr},
      refRegion, MINlevel);
}
// Outputs
// static
event writingFiles (i = 0, t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i+=100) {
  double ke = 0., zcm = 0., vcm = 0., wt = 0.;
  double zcm1 =0., vcm1 = 0., wt1 = 0.;
  foreach (reduction(+:ke), reduction(+:zcm), reduction(+:vcm), reduction(+:wt), reduction(+:zcm1), reduction(+:vcm1), reduction(+:wt1)){
    ke += 2*pi*y*(0.5*rho(f1[]+f2[]+f3[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    zcm += 2*pi*y*((f1[]+f3[])*x)*sq(Delta);
    vcm += 2*pi*y*((f1[]+f3[])*u.x[])*sq(Delta);
    wt += 2*pi*y*(f1[]+f3[])*sq(Delta);
    zcm1 += 2*pi*y*(f3[]*x)*sq(Delta);
    vcm1 += 2*pi*y*(f3[]*u.x[])*sq(Delta);
    wt1 += 2*pi*y*f3[]*sq(Delta);
  }
  zcm /= wt;
  vcm /= wt;
  zcm1 /= wt1;
  vcm1 /= wt1;
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke zcm vcm zcm1 vcm1\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke zcm vcm zcm1 vcm1\n");
    fprintf (fp, "%d %g %g %g %g %g %g %g\n", i, dt, t, ke, zcm, vcm, zcm1, vcm1);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g %g %g %g %g\n", i, dt, t, ke, zcm, vcm, zcm1, vcm1);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g %g %g %g %g\n", i, dt, t, ke, zcm, vcm, zcm1, vcm1);
}
