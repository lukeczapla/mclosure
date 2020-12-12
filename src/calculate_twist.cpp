#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bdna/bdna.h"

int main(int argc, char **argv) {

  dna *mydna;          

  if (argc != 3) {
    printf("usage: run [seq] [steps]\n");
    return -1;
  }
                                                                                                                            
  printf("opening %s\n", argv[1]);

  mydna = new dna(argv[1]);
                                                                                
  bdna *b = new bdna();
  b->initialize_bdna(*mydna);

  b->readconfig(argv[2]);


  matrix W = identity(4);

  for (int i = 0; i < b->nsteps; i++) {
    W = W * calculateW(b->v[i]);
  }

  double tw = b->calculate_twist_open();

  printf("new twist is %lf or %lf degrees or average %lf per step. Bend is %lf or %lf\n", tw, 360.0*tw, 360*tw/b->nsteps, 180.0*acos(W(3,3))/M_PI, W(3,3));

}
