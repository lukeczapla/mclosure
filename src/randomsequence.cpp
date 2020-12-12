#include <stdio.h>
#include <stdlib.h>
#include <time.h>

main(int argc, char **argv) {
  srand48(time(NULL));
  int length = 0;
  if (argc != 2) {
    printf("usage: randomsequence length\n");
    exit(0);
  }
  sscanf(argv[1], "%d", &length);
  for (int i = 0; i < length; i++) {
    int base = (int)(4.0*drand48());
    switch (base) {
      case 0: { printf("A"); break; }
      case 1: { printf("T"); break; }
      case 2: { printf("G"); break; }
      case 3: { printf("C"); break; }
      default: { break; }
    }
  }
  printf("\n");
}

