#include "closure/closure-HU.h"
#include <stdio.h>

int main(int argc, char **argv) {
  dna *d = new dna("seq");
  closure T;
  T.initialize_params();
  T.initialize_bdna(*d);
  T.set_beta(1.0);
  T.init_HU();

  matrix tp(6,1);
  tp.setv(1,1,0.0);
  tp.setv(2,1,0.0);
  tp.setv(3,1,34.28);
  tp.setv(4,1,0.0);
  tp.setv(5,1,0.0);
  tp.setv(6,1,3.4);


  double PHU = 0.0129;
  if (argc != 2) printf("Syntax: dynamic P_HU\n\n");
  else sscanf(argv[1], "%lf", &PHU);
  matrix W;

  matrix WPH[16];
  matrix WW = identity(4);

  matrix Wtotal = identity(4);

  matrix Wref = identity(4);

  for (int i = 0; i < 16; i++) {
    Wref = Wref * calculateW(T.get_HU(1).v[i]);
  }

  printf("Bend = \n");
  writematrix(stdout, Wref);
  writematrix(stdout, calculatetp(Wref));

  return 0;

}
