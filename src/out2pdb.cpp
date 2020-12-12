
#include "bdna/bdna.h"

int main(int argc, char **argv) {

  if (argc != 3) {
    printf("Usage: out2pdb <out file name> <pdb file name>\n");
    return 0;
  }  

  char *p;
  char s[256];

  FILE *f = fopen(argv[1], "r");

  int ns;
  do {
    p = fgets(s, 2048, f);
    if (p == NULL) {
      printf("Bad file format in structure file!\n");
      return 0;
    }
  } while (sscanf(s, "%d", &ns) != 1);

  bdna b(ns);

  fclose(f);


  b.readconfig(argv[1]);

  b.printpdb(argv[2]);


}
