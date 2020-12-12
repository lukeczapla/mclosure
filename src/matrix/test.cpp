#include "matrix.h"
#include <stdio.h>

main() {

  matrix A(2,2);

  A.setv(1,1, 1.0);
  A.setv(2,1, 1.0);
  A.setv(1,2, 0.0);
  A.setv(2,2, 1.0);

  matrix B = A * A * A;

  printf("%lf %lf\n", B(1,1), B(1,2));
  printf("%lf %lf\n", B(2,1), B(2,2));
}
