#include "bdna/bdna.h"
#include <stdio.h>
#include <stdlib.h>

main() {
  bdna b(15);

  b.readconfig("hybrid2.dat");

  b.print3dna("hybrid2-DNA.inp");
}
