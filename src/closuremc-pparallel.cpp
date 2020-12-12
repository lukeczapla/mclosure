
#define _PROTEINS

#include "defines.h"

#include "closure/closure-pparallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#include "clparams-parallel2.cpp"

closure *T;

void *calculate_combination(void *ptr) {
	int *x = (int *)ptr;
    T->reducedcombination(*x);
};


void *calculate_closure(void *ptr) {
	int *x = (int *)ptr;
    T->solveclosureproblem(*x, (num_steps/nnodes)*(*x));
};



main(int argc, char **argv) {

  dna *mydna;

  parsecommandline(argc, argv);

  printf("opening %s\n", filename);

  mydna = new dna(filename);

  if (key != 0) T = new closure(key);
    else T = new closure();

  T->set_beta(beta);
  T->initialize_params();

  T->set_pics(make_pictures);
  T->set_pics_sc(make_spictures);
  T->setepsilon(epsilon);
  T->setoverlap(overlap);
  T->setrandomsamples(osamples);

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
 
  T->setparams(radi, gam, twi);

  T->removeoverlap(neglect);


  T->setnodes(nnodes);
  printf("%d nodes\n", nnodes);
  T->initE(num_steps);

  printf("Spawning %d threads\n", nnodes);

  pthread_t threads[10];
  int x[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  printf("nnodes = %d\n", nnodes);
  for (int i = 0; i < nnodes; i++) pthread_create(&threads[i], NULL, calculate_closure, (void *) &x[i]);

  for (int i = 0; i < nnodes; i++)
    pthread_join(threads[i], NULL);

  printf("*** Program Part 2\n*** Spawning %d new threads\n", nnodes);

  for (int i = 0; i < nnodes; i++) pthread_create(&threads[i], NULL, calculate_combination, (void *) &x[i]);

  for (int i = 0; i < nnodes; i++)
    pthread_join(threads[i], NULL);

  T->writeJ(stdout);

  char outJfile[100];
  sprintf(outJfile, "Jfactor-%dbp-ID%d-allnodes", T->nsteps, T->MCkey());

  FILE *fJ = fopen(outJfile, "w");

  T->writeJ(fJ);

  T->cleanup();

  delete T;
  delete mydna;

  pthread_exit(NULL);

}

