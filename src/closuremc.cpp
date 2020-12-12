
#include "defines.h"

#include "closure/closure.h"

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

  T->setepsilon(epsilon);
  T->setmem(himem);

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
  T->setpics(make_pictures, make_spictures, sc_width); 
  T->setparams(radi, gam, twi);
  T->setoverlap(overlap);
  T->setrandomsamples(osamples);
  T->setrestrictions(restrictends);
  T->setbinsize(bin_width);

  if (persistence) {
	if (radi == 0.0) T->persistence(num_steps);
	else T->solveextension(num_steps, 0.8);//T->persistence(num_steps);
	return 0;
  }
 
  T->stress_free_state(0);
//   printf("Tw = %lf Wr = %lf\n", T->calculate_twist(), T->calculate_writhe());

  T->solveclosureproblem(num_steps, num_divisions, reduced);

//  delete T;
//  delete mydna;

}

