#include "bdna/bdna.h"
#include <stdio.h>
#include <stdlib.h>

main() {
  bdna b(15);
  bdna c(14);

  b.readconfig("hybrid.dat");
  c.readconfig("NHP-1.dat");
  matrix W = identity(4);
  for (int i = 0; i < 7; i++) {
    W = W * calculateW(b.v[i]);
  }

  matrix W2 = identity(4);
  matrix tp;
  for (int i = 0; i < 5; i++) {
    tp = c.v[i];
    W2 = W2 * calculateW(tp);
  }

  matrix Win = W2;

  printf("%f %f %f\n", Win(1,4), Win(2,4), Win(3,4));
  printf("%f %f %f\n", Win(1,1), Win(2,1), Win(3,1));
  printf("%f %f %f\n", Win(1,2), Win(2,2), Win(3,2));
  printf("%f %f %f\n", Win(1,3), Win(2,3), Win(3,3));

  Win = W;

  printf("%f %f %f\n", Win(1,4), Win(2,4), Win(3,4));
  printf("%f %f %f\n", Win(1,1), Win(2,1), Win(3,1));
  printf("%f %f %f\n", Win(1,2), Win(2,2), Win(3,2));
  printf("%f %f %f\n", Win(1,3), Win(2,3), Win(3,3));


}
