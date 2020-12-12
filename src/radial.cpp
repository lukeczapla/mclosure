
#include "defines.h"

#include "closure/closure-radial.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rparams.cpp"

int main(int argc, char **argv) {

  dna *mydna;

  parsecommandline(argc, argv);

  printf("opening %s\n", filename);

  mydna = new dna(filename);

  closure *T; 
  if (key != 0) T = new closure(key);
    else T = new closure();

  T->initialize_params();

  T->initialize_bdna(*mydna);

  T->printeigenvalues(); 

  T->stress_free_state(1);
//   printf("Tw = %lf Wr = %lf\n", T->calculate_twist(), T->calculate_writhe());

  T->setbinsize(bin_width);
  T->combination(num_steps, radialfilename);

  delete T;
  delete mydna;

}

