#include "closure/closure-HU.h"
#include <stdio.h>

int main(int argc, char **argv) {
  dna *d = new dna("seq");
  closure T;
  T.initialize_params();
  T.initialize_bdna(*d);
  T.set_beta(1.0);
  T.init_HU();

  bdna HMG(8);
  HMG.readconfig("HMG-2.dat"); 

  double PHU = 0.0129;

  if (argc != 2) printf("Syntax: overall P_HU\n\n");
  else
    sscanf(argv[1], "%lf", &PHU);

  PHU = 8.0/((1.0/PHU)+7.0);
 
  matrix W, Wp;

  matrix WW = identity(4);
  for (int j = 0; j < 8; j++) {

    WW = WW*calculateW(HMG.v[j]);
 
  }

  writematrix(stdout, WW);  
  writematrix(stdout, calculatetp(WW));

  double q = PHU/1.0;

  Wp = q*WW;
  writematrix(stdout, Wp);

  matrix Wsum;
  for (int i = 0; i < 50000; i++) {
    W = identity(4);
    T.generate_config1();

    T.fast_translate(0,8);
    for (int j = 0; j < 8; j++) {
      W = W*calculateW(T.v[j]);
    }
    if (i == 0) {
      Wsum = W;
      writematrix(stdout, calculatetp(Wsum));
    }
    else Wsum = Wsum + W;
  }

  q = (1.0-PHU)/50000.0;
  WW = q*Wsum;
  writematrix(stdout, WW);

  matrix M = WW + Wp;
  writematrix(stdout, M);
  matrix WQ = identity(4);

  for (int i = 0; i < 50000; i++)
    WQ = WQ * M;

  writematrix(stdout, WQ);

  printf("result = %lf\n", WQ(3,4));

  return 0;

}
