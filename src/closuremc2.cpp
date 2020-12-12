
#include "defines.h"

#include "closure/closure2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clparams.cpp"

main(int argc, char **argv) {

  printf("closuremc2: SPECIFICALLY FOR CALCULATIONS INVOLVING HU PROTEIN.\n");

  dna *mydna;

  parsecommandline(argc, argv);

  printf("opening %s\n", filename);

  mydna = new dna(filename);

  closure *T;
  if (key != 0) T = new closure(key);
    else T = new closure();

  T->set_beta(beta);
  T->initialize_params();

  T->setepsilon(epsilon);
  T->set_HU(PHU);

//  T->settubelength(t);
  
  matrix boundary(6,1);

  if (j != 0) T->setboundarynum(j);
    else T->setboundarynum(1);

  for (int i = 0; i < j; i++) {
    for (int p = 0; p < 6; p++) boundary.setv(p+1, 1, bound[i][p]);
    T->addboundary(boundary);
  }

  if (j == 0) {
    printf("no boundaries defined!  defaulting to closed\n");
    T->addboundary(ringclosureboundary());
  }

  T->initialize_bdna(*mydna);
 
  T->setparams(radi, gam, twi);

  T->solveclosureproblem(num_steps, num_divisions, reduced);

  delete T;
  delete mydna;

}

