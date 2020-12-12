#include "defines.h"
#include "closure/closure-proteins.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>

int main(int argc, char **argv) {

  long int TIME = time(NULL) % 30000;

  dna *d = new dna("seq");
  unsigned short k[] = {0,0,TIME};
  closure *T = new closure(k);
  T->set_beta(1.0);

  T->initialize_bdna(*d);
  T->initialize_params();
  T->init_proteins(); 

  T->printeigenvalues();

  matrix W;

  matrix WW = identity(4);
  matrix Wsum;
  T->stress_free_state(1);
//  T->printparams();
  FLOAT rmstwist = 0.0;
  for (int i = 0; i < 20000; i++) {
    W = identity(4);
    T->generate_config1();
    T->generate_config2();
    T->fast_translate(0, T->nsteps);
//    if (i == 0) for (int k = 0; k < T->nsteps; k++) writematrix(stdout, T->v[k]);
    for (int j = 0; j < T->nsteps; j++) {
      W = W*calculateW(T->v[j]);
      rmstwist += (T->v[j](3,1)-34.28)*(T->v[j](3,1)-34.28)/T->nsteps;
    }
    if (i == 0) Wsum = W;
    else Wsum = Wsum + W;
  }
  rmstwist /= 20000;
  double q = 1.0/20000.0;
  W = q*Wsum;
  writematrix(stdout, W);

  matrix M = W;
  matrix WQ = identity(4);

  for (int i = 0; i < 50000; i++)
    WQ = WQ * M;

  writematrix(stdout, WQ);

  printf("result = %lf, rmstwist = %lf\n", WQ(3,4), sqrt(rmstwist));

  return 0;

}
