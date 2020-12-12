#ifndef _PLOT_H
#define _PLOT_H

#include "../defines.h"


class plot {

  FILE *f;
  FLOAT buffer[1024];

  int curr;
  int pass;

 public:

  plot(char *fn);
  ~plot();

  void add_plot(FLOAT d);

};

#endif
