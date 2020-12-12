
#ifndef _RADIAL_H
#define _RADIAL_H

#define Xsze 100
#define Ysze 100
#define Zsze 350
#define rsze 500
#define ritr 5.0

#include "../matrix/matrix.h"

class spatial {

  int threshold;	
  FLOAT size;

 public:

  int *density[8];

  spatial(FLOAT sz, int thres);
  ~spatial();
  void add_spatial(matrix *W);

  void printpdb(char *filename);

};

#endif
