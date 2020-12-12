#include <stdio.h>
#include "bdna/bdna.h"

int main(int argc, char **argv) {

  dna *mydna;                                                                                                                                      
  printf("opening %s\n", argv[1]);

  mydna = new dna(argv[1]);
                                                                                
  bdna *b = new bdna();
  b->initialize_bdna(*mydna);

  b->readconfig("out");

  printf("new twist is %lf\n", b->calculate_twist());

}
