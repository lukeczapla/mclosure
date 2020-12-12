#ifndef _ENERGY_H
#define _ENERGY_H

#define WRITE_ENERGY 0
#define READ_ENERGY 1

#include <stdio.h>
#include "../defines.h"

class energy {

  FLOAT buffer[1024];
  int curr;
  int pass;
  FILE *ener;

 public:

  energy(char *filename, int mode);
  ~energy();

  void add_energy(FLOAT enr);
  int read_energy(FLOAT *enr);

};

#endif
