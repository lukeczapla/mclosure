#include <stdio.h>
#include <stdlib.h>

#include "bdna/bdna.h"

main(int argc, char **argv) {
  matrix pretp(6,1);
  matrix posttp(6,1);
  matrix Wbase;
  matrix W = identity(4);
  posttp.setv(1,1,  -4.8233);
  posttp.setv(2,1,  58.0204);
  posttp.setv(3,1,  64.6245);
  posttp.setv(4,1,  17.8621);
  posttp.setv(5,1,  -2.2528);
  posttp.setv(6,1,  74.5819);
  if (argc == 2) {
  printf("%s\n", argv[1]);
  FILE *f = fopen(argv[1], "r");
  float x, y, z;
  fscanf(f, "%f %f %f", &x, &y, &z);
  W.setv(1,4,x);
  W.setv(2,4,y);
  W.setv(3,4,z);
  W.setv(4,4,1.0);

  fscanf(f, "%f %f %f", &x, &y, &z);
  W.setv(1,1,x);
  W.setv(2,1,y);
  W.setv(3,1,z);
  W.setv(4,1,0.0);
  
  fscanf(f, "%f %f %f", &x, &y, &z);
  W.setv(1,2,x);
  W.setv(2,2,y);
  W.setv(3,2,z);
  W.setv(4,2,0.0);

  fscanf(f, "%f %f %f", &x, &y, &z);
  W.setv(1,3,x);
  W.setv(2,3,y);
  W.setv(3,3,z);
  W.setv(4,3,0.0);

  writematrix(stdout, W);
  writematrix(stdout, W*calculateW(posttp));
  writematrix(stdout, calculatetp(W*calculateW(posttp)));

  fclose(f);    
  return 0;
  }
  pretp.setv(1,1,  4.8158);
  pretp.setv(2,1,  58.0379);
  pretp.setv(3,1,  64.6256);
  pretp.setv(4,1,  -17.8666);
  pretp.setv(5,1,  -2.2604);
  pretp.setv(6,1,  74.5757);
  posttp.setv(1,1,  -4.8233);
  posttp.setv(2,1,  58.0204);
  posttp.setv(3,1,  64.6245);
  posttp.setv(4,1,  17.8621);
  posttp.setv(5,1,  -2.2528);
  posttp.setv(6,1,  74.5819);


  Wbase = calculateW(pretp);
  printf("%lf %lf %lf\n", Wbase(1,4), Wbase(2,4), Wbase(3,4));
  printf("%lf %lf %lf\n", Wbase(1,1), Wbase(2,1), Wbase(3,1));
  printf("%lf %lf %lf\n", Wbase(1,2), Wbase(2,2), Wbase(3,2));
  printf("%lf %lf %lf\n", Wbase(1,3), Wbase(2,3), Wbase(3,3));

  

  writematrix(stdout, Wbase);

  Wbase = invert(calculateW(posttp));
  printf("%lf %lf %lf\n", Wbase(1,4), Wbase(2,4), Wbase(3,4));
  printf("%lf %lf %lf\n", Wbase(1,1), Wbase(2,1), Wbase(3,1));
  printf("%lf %lf %lf\n", Wbase(1,2), Wbase(2,2), Wbase(3,2));
  printf("%lf %lf %lf\n", Wbase(1,3), Wbase(2,3), Wbase(3,3));

  writematrix(stdout, Wbase);
  writematrix(stdout, invert(Wbase));
  writematrix(stdout, calculatetp(Wbase));
}
