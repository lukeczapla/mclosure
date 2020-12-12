#include "matrix/matrix.h"

main() {

  matrix q = identity(3);

  q.setv(1,1, 0.0427);
  q.setv(2,2, 0.0597);
  q.setv(1,2, 0.03);
  q.setv(2,1, 0.03);
  q.setv(2,3, -0.1);
  q.setv(3,2, -0.1);

  double d[3];
  matrix eigenv = jeigen(q, d);

  for (int i = 0; i < 3; i++) printf("%lf\n", d[i]);

  writematrix(stdout, eigenv);

}
