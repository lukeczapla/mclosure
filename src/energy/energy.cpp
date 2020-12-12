
#include "energy.h"

#include <stdio.h>
#include <stdlib.h>

energy::energy(char *filename, int mode) {
  if (mode == WRITE_ENERGY) {
    curr = 0;
    pass = 0;
    ener = fopen(filename, "w");
  }
  if (mode == READ_ENERGY) {
    ener = fopen(filename, "r");
    if (ener == NULL) {
      printf("Could not open the energy file %s\n", filename);
      exit(0);
    }
  }
}

energy::~energy() {
  fwrite(buffer, sizeof(FLOAT), curr, ener);
  fclose(ener);
}

void energy::add_energy(FLOAT enr) {
  buffer[curr++] = enr;
  if (curr == 1024) {
    fwrite(buffer, sizeof(FLOAT), 1024, ener);
    curr = 0;
    pass++;
  }
}

int energy::read_energy(FLOAT *enr) {
  return fread(enr, sizeof(FLOAT), 1, ener);
}
