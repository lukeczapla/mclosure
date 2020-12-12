
#include "defines.h"

#include "closure/closure-parallel.h"


#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clparams-parallel.cpp"

closure *T; 


void *calculate_closure(void *ptr) {
	int *x = (int *)ptr;
    T->solveclosureproblem(*x, (num_steps/nnodes)*(*x));
};


void *calculate_combinations(void *ptr) {
	int *x = (int *)ptr;
    T->reducedcombination(*x);
};

/*
void *wait_completion(void *ptr) {
    int *x = (int *)ptr;
	printf("waiting thread() %d\n", *x);
    do {
      sleep(1);
    } while (!done);
    T->reducedcombination(*x);
}
*/

int main(int argc, char **argv) {

  dna *mydna;

  parsecommandline(argc, argv);

  printf("opening %s\n", filename);

  mydna = new dna(filename);

  if (key != 0) T = new closure(key);
    else T = new closure();

  T->set_beta(beta);
  T->initialize_params();

  T->setepsilon(epsilon);
  T->setmem(himem);

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
  T->setpics(make_pictures, make_spictures, sc_width); 
  T->setparams(radi, gam, twi);
  T->setoverlap(overlap);
  T->setrandomsamples(osamples);
  T->setrestrictions(restrictends);
  T->setbinsize(bin_width);
  T->setnodes(nnodes);
 
  T->stress_free_state(1);
//   printf("Tw = %lf Wr = %lf\n", T->calculate_twist(), T->calculate_writhe());

  T->initE(num_steps);


  printf("Spawning %d threads\n", nnodes);

  pthread_t threads[10];
  int x[] = {0,1,2,3,4,5,6,7,8,9};
  printf("nnodes = %d\n", nnodes);
  for (int i = 0; i < nnodes; i++) pthread_create(&threads[i], NULL, calculate_closure, (void *) &x[i]);

  for (int i = 0; i < nnodes; i++)
    pthread_join(threads[i], NULL);

  printf("\nstep 2: generating second half segments and combining\n");

  for (int i = 0; i < nnodes; i++) pthread_create(&threads[i], NULL, calculate_combinations, (void *) &x[i]);

  for (int i = 0; i < nnodes; i++)
    pthread_join(threads[i], NULL);

  T->writeJ(stdout);

  T->cleanup();

  delete T;
//  delete mydna;

}

