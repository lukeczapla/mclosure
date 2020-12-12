
#include "trajectory/trajectory.h"
#include <stdio.h>


main(int argc, char **argv) {

  if (argc != 4) {
    printf("Usage: trj2pdb <filename> <frame number> <output file>\n");
    return 0;
  }

  trajectory T(argv[1], 0, READ_TRAJECTORY);

  bdna b(T.nsteps);

  int n;

  sscanf(argv[2], "%d", &n);

  for (int i = 0; i < n; i++) {
    int r = T.read_coordinates(&b);
    if (r == 0) {
      printf("Frame does not exist!\n");
      return 0;
    }
  }

  T.read_coordinates(&b);

  b.printconfig(argv[3]);

}
