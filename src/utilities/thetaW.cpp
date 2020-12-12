#include <stdio.h>
#include "bdna/bdna.h"

main() {
  FLOAT a, b, c, d, e, f;
  scanf("%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f);
  matrix tp(6,1);
  tp.setv(1, 1, a);
  tp.setv(2, 1, b);
  tp.setv(3, 1, c);
  tp.setv(4, 1, d);
  tp.setv(5, 1, e);
  tp.setv(6, 1, f);
  matrix W = calculateW(tp);
  writematrix(stdout, W);
}
