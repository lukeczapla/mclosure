#include "closure/closure.h"
#include <stdio.h>

int main() {
  dna *d = new dna("seq");
  closure *T = new closure();
  T->initialize_params();
  T->initialize_bdna(*d);
  T->set_beta(1.0);

  matrix W = identity(4);

  T->stress_free_state(1);

  printf("twist %lf writhe %lf\n", T->calculate_twist_print(), T->calculate_writhe());

  T->print3dna("out.dat");

  for (int i = 0; i < T->nsteps; i++) {
    printf("theta_0(%d) = %lf %lf %lf %lf %lf %lf\n", i, T->v[i](1,1), T->v[i](2,1), T->v[i](3,1), T->v[i](4,1), T->v[i](5,1), T->v[i](6,1));
    W = W * calculateW(T->v[i]);
    writematrix(stdout, W);
    printf("\n\n");
  }

  delete d;
  delete T;
 
}



