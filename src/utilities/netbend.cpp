#include "defines.h"
#include "closure/closure-proteins.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char **argv) {
  dna *d = new dna("seq");
  unsigned short k[] = {0,0,1};
  closure *T = new closure(k);
  T->set_beta(1.0);

  T->initialize_bdna(*d);
  T->initialize_params();
  T->init_proteins(); 


  matrix tp;

  matrix W = identity(4);

  char datafile[100];

//  T->resize(numsteps);
  for (int j = 1; j < argc; j++) {
  T->readconfig(argv[j]);

  printf("running %s\n", argv[j]);

  FLOAT tw = 0.0;

  W = identity(4);
  for (int i = 0; i < T->nsteps; i++) {
    W = W * calculateW(T->v[i]);
    tw += T->v[i](3,1);
  }

  tp = calculatetp(W);

  printf("\n\nBend angle = %lf, Avg twist angle = %lf\n\n", sqrt(tp(1,1)*tp(1,1)+tp(2,1)*tp(2,1)), tw/T->nsteps);

  writematrix(stdout, calculatetp(W));

  }

  return 1;

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
