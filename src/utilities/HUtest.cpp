#include "bdna/bdna.h"
#include <math.h>

main () {
  bdna mydna(15);
  mydna.readconfig("HU.dat");
  matrix W = identity(4);
  for (int i = 0; i < 15; i++)
    W = W * calculateW(mydna.v[i]);
  writematrix(stdout, W);
  writematrix(stdout, calculatetp(W));
  printf("\nr = %lf\n\n", sqrt(W(1,4)*W(1,4)+W(2,4)*W(2,4)+W(3,4)*W(3,4)));
  printf("\ngamma = %lf\n\n", acos(W(3,3))*180.0/M_PI);
  printf("theta = %lf\n\n", acos((W(2,2)+W(1,1))/(1+W(3,3)))*180.0/M_PI);
}
