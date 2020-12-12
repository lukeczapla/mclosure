
#include "file.h"
#include "trajectory/trajectory.h"
#include <stdio.h>


main(int argc, char **argv) {

  if (argc != 5) {
    printf("Usage: trj2dat <seq filename> <trajectory filename> <frame number> <output file>\n");
    return 0;
  }

  trajectory T(argv[2], 0, READ_TRAJECTORY);

  dna mydna(argv[1]);

  bdna b;

  b.initialize_bdna(mydna);

  int n;

  sscanf(argv[3], "%d", &n);

  for (int i = 0; i < n; i++) {
    int r = T.read_coordinates(&b);
    if (r == 0) {
      printf("Frame does not exist!\n");
      return 0;
    }
  }

  T.read_coordinates(&b);

  FILE *f = openfwrite(argv[4]);

  fprintf(f, "%d\n%d\nshift slide rise roll tilt twist\n", b.nsteps+1, 0);

  int base = b.bpstep[0]/4;

  fprintf(f, "%c-%c %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf\n", btranslate(base), btranslate(3-base), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  
  for (int i = 0; i < b.nsteps; i++) {
    base = b.bpstep[i] % 4;
    fprintf(f, "%c-%c %7.3lf %7.3lf %7.3lf %7.3lf %7.3lf %7.3lf\n", btranslate(base), btranslate(3-base), b.v[i](4,1), b.v[i](5,1), b.v[i](6,1),
	    b.v[i](1,1), b.v[i](2,1), b.v[i](3,1));
  }

  fclose(f);

}
