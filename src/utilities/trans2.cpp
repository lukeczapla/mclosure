#include "bdna/bdna.h"
#include <stdio.h>
#include <stdlib.h>

main() {
  bdna b(15);

  b.readconfig("hybrid.dat");

  matrix W = identity(4);
  for (int i = 0; i < 7; i++) {
    W = W * calculateW(b.v[i]);
  }

  matrix Win = W;

  printf("%f %f %f\n", Win(1,4), Win(2,4), Win(3,4));
  printf("%f %f %f\n", Win(1,1), Win(2,1), Win(3,1));
  printf("%f %f %f\n", Win(1,2), Win(2,2), Win(3,2));
  printf("%f %f %f\n", Win(1,3), Win(2,3), Win(3,3));
}
