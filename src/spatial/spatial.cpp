
#include "../defines.h"
#include "../file.h"
#include "spatial.h"
#include <math.h>
#include <stdio.h>

spatial::spatial(FLOAT sz, int thres) {
  size = sz;
  threshold = thres;
  for (int i = 0; i < 8; i++) {
    density[i] = new int[Xsze*Ysze*Zsze];
  }
  for (int i = 0; i < Xsze*Ysze*Zsze; i++) {
    for (int j = 0; j < 8; j++)
      density[j][i] = 0;
  }
}

spatial::~spatial() {
  for (int i = 0; i < 8; i++) delete [] density[i];
}

void spatial::add_spatial(matrix *W) {
  FLOAT x, y, z;
  x = (*W)(1,4);
  y = (*W)(2,4);
  z = (*W)(3,4);
  int index = 0;
  if (x >= 0.0) index += 1;
  if (y >= 0.0) index += 2;
  if (z >= 0.0) index += 4;
  int i = (int)(fabs(x)/size) + Xsze*(int)(fabs(y)/size)+Xsze*Ysze*(int)(fabs(z)/size);
  if (i < Xsze*Ysze*Zsze)
    density[index][i]++;
}

void spatial::printpdb(char *filename) {

  FILE *f = fopen(filename, "a");
  if (f == NULL) {
    printf("Spatial PDB file %s could not be opened\n", filename);
    return;
  }

  int n = 1;
  for (int i = 0; i < Xsze; i++) for (int j = 0; j < Ysze; j++) for (int k = 0; k < Zsze; k++) {
    for (int l = 0; l < 8; l++) {
      if (density[l][i+Xsze*j+Xsze*Ysze*k] > threshold) {
        FLOAT X = i*size+size/2.0;
        FLOAT Y = j*size+size/2.0;
	FLOAT Z = k*size+size/2.0;
        if ((l % 2) == 0) X = -X;
        if ((l % 4) < 2) Y = -Y;
        if ((l / 4) == 0) Z = -Z;
        char atom = 'N';
        switch ((int)log10(density[l][i+Xsze*j+Xsze*Ysze*k]-(FLOAT)threshold)) {
	case 0: { atom = 'F'; break; }
	case 1: { atom = 'S'; break; }
        case 2: { atom = 'O'; break; }
	case 3: { atom = 'O'; break; }
        default: { atom = 'O'; break; }
        }
	fprintf(f, "ATOM%7d  %c   UNK     2    %8.3lf%8.3lf%8.3lf  %d\n", n++, atom, X, Y, Z, density[l][i+Xsze*j+Xsze*Ysze*k]);
      }
    }
  }
  fclose(f);
}

