/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f3[];
char filename[80];

int main(int a, char const *arguments[]){
  f3[left] = dirichlet(1.);

  sprintf(filename, "%s", arguments[1]);
  restore (file = filename);
  f3.prolongation = fraction_refine;
  boundary((scalar *){f3});

  FILE * fp = ferr;
  output_facets(f3, fp);
  fflush (fp);
  fclose (fp);
}
