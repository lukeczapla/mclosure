#include "bdna/bdna.h"

main() {
  bdna b(15);
  b.readconfig("HMGB1.dat");
  matrix W = identity(4);
  for (int i = 0; i < 15; i++) {
    W = W * calculateW(b.v[i]);
  }
  writematrix(stdout, calculatetp(W));
}
