#include <stdio.h>
#include "bdna/bdna.h"

main() {
  FLOAT a, b, c, d;
  matrix W(4,4);
  scanf("%lf %lf %lf %lf", &a, &b, &c, &d);
  W.setv(1,1, a);
  W.setv(1,2, b);
  W.setv(1,3, c);
  W.setv(1,4, d);
  scanf("%lf %lf %lf %lf", &a, &b, &c, &d);
  W.setv(2,1, a);
  W.setv(2,2, b);
  W.setv(2,3, c);
  W.setv(2,4, d);
  scanf("%lf %lf %lf %lf", &a, &b, &c, &d);
  W.setv(3,1, a);
  W.setv(3,2, b);
  W.setv(3,3, c);
  W.setv(3,4, d);
  W.setv(4,1, 0.0);
  W.setv(4,2, 0.0);
  W.setv(4,3, 0.0);
  W.setv(4,4, 1.0);

  matrix tp = calculatetp(W);
  writematrix(stdout, tp);
}
