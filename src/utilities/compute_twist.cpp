#include <stdio.h>
#include "bdna/bdna.h"

main(int argc, char **argv) {
  if (argc != 3) {
    printf("usage: compute_link seq_file out_file"); 
  }

  dna *d = new dna(argv[1]);

  bdna *b = new bdna(d->nbp-1);
  b->initialize_bdna(*d);

  b->readconfig(argv[2]);

  for (;;)
    printf("%f %f\n", b->calculate_twist(), b->calculate_writhe());
  

}
