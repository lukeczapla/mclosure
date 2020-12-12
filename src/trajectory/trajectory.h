#ifndef _MCTRAJECTORY_H
#define _MCTRAJECTORY_H

#include <stdio.h>
#include "../defines.h"
#include "../bdna/bdna.h"

#define WRITE_TRAJECTORY 0
#define READ_TRAJECTORY 1

class trajectory {

  int tfile;

 public:

  char *fname;

  int nsteps;

  trajectory(char *filename, int size, int mode);
  ~trajectory();

  void add_coordinates(bdna *b);
  int read_coordinates(bdna *b);
  int read_coordinates(bdna *b, int n);
  int read_coordinates(long long pos, bdna *b);
  void seek_coordinate(int n);
};

#endif

