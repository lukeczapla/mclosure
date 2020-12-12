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

  matrix WPH[14];
  matrix WW;

  matrix Wtotal = identity(4);

  matrix Wref = identity(4);

  for (int i = 0; i < 14; i++) {
    Wref = Wref * calculateW(T.get_HU(0).v[i]);
  }

  writematrix(stdout, Wref);

  for (int i = 0; i < 4; i++) {
  WW = identity(4);
  for (int j = 0; j < 14; j++) {

    WPH[j] = calculateW(T.get_HU(i).v[j]);
    WW = WW*WPH[j];
 
  }

  double pp = 1.0/4.0;
  WW = pp * WW;

  WW = WW*invert(Wref);
  for (int j = 0; j < 14; j++) WW = WW * calculateW(tp);
  if (i == 0) Wtotal = WW;
  else Wtotal = Wtotal + WW;
  }
  
  writematrix(stdout, Wtotal);

  double q = PHU;


  matrix WQ = identity(4);

  for (int i = 0; i < 50000; i++)
    WQ = WQ * Wtotal;

  writematrix(stdout, WQ);

  printf("result %lf\n", PHU*WQ(3,4)+(1-PHU)*480.0);

  return 0;

}
