
#define _PROTEINS

#include "defines.h"

#include "closure/closure-proteins.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clparams.cpp"

int main(int argc, char **argv) {

  dna *mydna;

  parsecommandline(argc, argv);

  printf("opening %s\n", filename);

  mydna = new dna(filename);

  closure *T; 
  if (key != 0) T = new closure(key);
    else T = new closure();

  T->set_beta(beta);
  T->initialize_params();

  T->set_pics(make_pictures);
  T->set_pics_sc(make_spictures);
  T->setepsilon(epsilon);
  T->setoverlap(overlap);
  T->setrandomsamples(osamples);

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

  T->removeoverlap(neglect);

  if (persistence) {
    T->persistence(num_steps);
    return 1;
  }

  T->solveclosureproblem(num_steps, num_divisions, reduced);

  delete T;
  delete mydna;

}

