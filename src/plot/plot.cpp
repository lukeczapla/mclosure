
#include <stdio.h>
#include "plot.h"
#include "../file.h"

plot::plot(char *fn) {
  f = openfwrite(fn);
  curr = 0;
  pass = 0;
}

plot::~plot() {
  for (int i = 0; i < curr; i++) fprintf(f, "%d %le\n", pass*1024+i, buffer[i]);
  fclose(f);
}

void plot::add_plot(FLOAT d) {
  buffer[curr++] = d;
  if (curr == 1024) {
    for (int i = 0; i < 1024; i++) {
      fprintf(f, "%d ", pass*1024+i);
      fprintf(f, FLOAT_ESTR, buffer[i]);
      fprintf(f, "\n");
    }
    curr = 0;
    pass++;
  }
}
