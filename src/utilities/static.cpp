#include "closure/closure-HU.h"
#include <stdio.h>

int main(int argc, char **argv) {
  dna *d = new dna("seq");
  closure T;
  T.initialize_params();
  T.initialize_bdna(*d);
  T.set_beta(1.0);
  T.init_HU();
 

  double PHU = 0.0129;

  if (argc != 2) printf("Syntax: program P_HU\n\n");
  else
    sscanf(argv[1], "%lf", &PHU);
 
  matrix W,WW;
  double q;

for (int i = 0; i < 4; i++) {

  matrix tW;
  WW = identity(4);
  for (int j = 0; j < 14; j++) {

  tW = calculateW(T.get_HU(i).v[j]);
  WW = WW*tW;
 
  }
  q = PHU/8.0;
  if (i == 0) W = q*WW;
  else W = W + q*WW;
  WW = identity(4);
  matrix tp;
  for (int j = 0; j < 14; j++) {
    tp = T.get_HU(i).v[14-j-1];
    tp.setv(1,1,-tp(1,1));
    tp.setv(4,1,-tp(4,1));
    WW = WW*calculateW(tp);
  }
  W = W + q*WW;
}
  writematrix(stdout, W);

  matrix tp(6,1);
  tp.setv(1,1,0.0);
  tp.setv(2,1,0.0);
  tp.setv(3,1,34.28);
  tp.setv(4,1,0.0);
  tp.setv(5,1,0.0);
  tp.setv(6,1,3.4);

  matrix C = calculateW(tp);
  writematrix(stdout, C);
  q = 1.0-PHU;
  matrix M = q*C;
  writematrix(stdout, M);

  WW = M + W;
  writematrix(stdout, WW);
  matrix WQ = identity(4);

  for (int i = 0; i < 50000; i++)
    WQ = WQ * WW;

  writematrix(stdout, WQ);

  printf("result = %lf\n", WQ(3,4));

  delete d;
  return 0;

}
