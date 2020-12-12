#include "bdna/bdna.h"
#include <math.h>

main () {
  bdna mydna(8);
  mydna.readconfig("HMG.dat");
  mydna.calculate_twist_open();
  matrix W = identity(4);
  for (int i = 0; i < 8; i++)
    W = W * calculateW(mydna.v[i]);
  writematrix(stdout, W);
  writematrix(stdout, calculatetp(W));
  printf("\nr = %lf\n\n", sqrt(W(1,4)*W(1,4)+W(2,4)*W(2,4)+W(3,4)*W(3,4)));
  printf("\ngamma = %lf\n\n", acos(W(3,3))*180.0/M_PI);
  printf("theta = %lf\n\n", acos((W(2,2)+W(1,1))/(1+W(3,3)))*180.0/M_PI);
  W = identity(4);
  matrix tpt(6,1);
  tpt.setv(1,1, 0.0);
  tpt.setv(2,1, 0.0);
  tpt.setv(3,1, 34.28);
  tpt.setv(4,1, 0.0);
  tpt.setv(5,1, 0.0);
  tpt.setv(6,1, 3.4);
  writematrix(stdout, calculatetp(calculateW(calculatetp(calculateW(tpt)))));
  for (int i = 0; i < 8; i++)
    W = W * calculateW(tpt);
  writematrix(stdout, calculatetp(W));
}
