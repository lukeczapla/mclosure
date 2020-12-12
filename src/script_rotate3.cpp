#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "translate_tools.h"
#include "bdna/bdna.h"

#define D2R 1.745329251994e-2

int main(int argc, char **argv) {
  matrix W(4,4);
  read_reference("ref_up", &W);
  matrix R = invert(W);
  writematrix(stdout, R);
  printf("written rotation file\n");
  return 0;
}
