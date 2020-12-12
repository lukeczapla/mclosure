
#include <stdio.h>
#include <stdlib.h>

#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

#include "trajectory.h"
#include "../defines.h"
#include "../file.h"

trajectory::trajectory(char *filename, int size, int mode) {
  fname = new char[30];
  strcpy(fname, filename);
  if (mode == WRITE_TRAJECTORY) {
    tfile = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0777);
    write(tfile, &size, sizeof(size));
    nsteps = size;
  }
  if (mode == READ_TRAJECTORY) {
    tfile = open(filename, O_RDONLY);
    if (tfile == -1) {
      printf("Could not open the trajectory file %s\n", filename);
      exit(0);
    }
    int q = read(tfile, &nsteps, sizeof(nsteps));
    if (q == 0) {
      printf("Empty trajectory file\n");
      exit(0);
    }
    if ((nsteps != size) && (size != 0)) {
      printf("Warning: trajectory dna sizes do not match!\n");
    }
  }
}

trajectory::~trajectory() {
  close(tfile);
}

void trajectory::add_coordinates(bdna *b) {
  if (b->nsteps < nsteps) {
    printf("Could not add to trajectory, dna is not the correct size!");
    return;
  }
  FLOAT *d = new FLOAT[6*nsteps];
  for (int i = 0; i < nsteps; i++) {
    b->v[i].getdata(&d[6*i]);
  }
  write(tfile, d, 6*nsteps*sizeof(FLOAT));
  delete [] d;
}


int trajectory::read_coordinates(bdna *b) {
  if (nsteps == 0) return 0;
  FLOAT *d = new FLOAT[6*nsteps];
  read(tfile, d, 6*nsteps*sizeof(FLOAT));
  for (int i = 0; i < nsteps; i++) {
    /*
    for (int j = 0; j < 6; j++) {
      b->v[i].setv(j+1, 1, q[6*i+j]);
    }
    */
    b->v[i].setdata(&d[6*i]);
  }
  delete [] d;
  return 1;
}

int trajectory::read_coordinates(long long pos, bdna *b) {
  int realpos = pos - pos % 512;
  int reallength = 6*sizeof(FLOAT)*b->nsteps+pos % 512;
  char *p = (char *)mmap(0, reallength, PROT_READ, MAP_PRIVATE, tfile, realpos);

  FLOAT *d = (FLOAT *)(p + pos % 512);  

  for (int i = 0; i < b->nsteps; i++) {
    /*
    for (int j = 0; j < 6; j++) {
      b->v[i].setv(j+1, 1, q[6*i+j]);
    }
    */
    b->v[i].setdata(&d[6*i]);
  }
  munmap(p, reallength);
  return 1;
  

}

int trajectory::read_coordinates(bdna *b, int n) {
  lseek(tfile, n*6*nsteps*sizeof(FLOAT)+sizeof(int), SEEK_SET);
  return read_coordinates(b);
}

void trajectory::seek_coordinate(int n) {
  lseek(tfile, n*6*nsteps*sizeof(FLOAT)+sizeof(int), SEEK_SET);
}
