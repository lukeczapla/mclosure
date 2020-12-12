
#include <stdio.h>
#include <math.h>
#include "trajectory/trajectory.h"
#include "plot/plot.h"
#include "spatial/spatial.h"
#include "histogram/histogram.h"

#include "trjparams.cpp"

main(int argc, char **argv) {

  plot *etime;
  plot *gtime;
  plot *rtime;

  histogram *ehist;
  histogram *ghist;
  histogram *rhist;

  parsecommandline(argc, argv);
  
  if (igt) {
    gtime = new plot(gt);
  }

  if (iet) {
    etime = new plot(et);
  }

  if (irt) {
    rtime = new plot(rt);
  }

  if (igh) {
    ghist = new histogram(2000, 0.1);
  }

  if (ieh) {
    ehist = new histogram(5000, 1.0);
  }

  if (irh) {
    rhist = new histogram(5000, 5.0);
  }

  if (argc == 1) {
    printf("Usage: calctrj <filename> <options>\n");
    return 0;
  }

  trajectory T(source, 0, READ_TRAJECTORY);

  dna mydna(filename);

  bdna b;
  b.initialize_bdna(mydna);

  matrix W;

  int i = 0;

  FLOAT rad;
  FLOAT e;

  printf("Processing. . .\n");

  for (int j = 0; j < start; j++) {
    if (!T.read_coordinates(&b)) {
      printf("Starting frame does not exist!\n");
      exit(0);
    }
    i++;
  }

  while (T.read_coordinates(&b) && (i <= end)) {

    if (iet || ieh) e = b.calculate_energy(W);

    if ((irt || irh) && !(iet || ieh)) {
      W = identity(4);
      for (int i = 0; i < b.nsteps; i++) W = W*calculateW(b.v[i]);
    }

    if (irt || irh) rad = sqrt(W(1,4)*W(1,4)+W(2,4)*W(2,4)+W(3,4)*W(3,4));
    if (igt) gtime->add_plot(b.calculateR());
    if (iet) etime->add_plot(e);
    if (irt) rtime->add_plot(rad);
    if (igh) ghist->add_data(b.calculateR());
    if (ieh) ehist->add_data(e);
    if (irh) rhist->add_data(rad);
    i++;
  }

  if (igh) ghist->printhistogram(gh);
  if (ieh) ehist->printhistogram(eh);
  if (irh) rhist->printhistogram(rh);

  printf("%d records read\n", i-start);

  if (igt) delete gtime;
  if (iet) delete etime;
  if (irt) delete rtime;
  if (igh) delete ghist;
  if (ieh) delete ehist;
  if (irh) delete rhist;

}
