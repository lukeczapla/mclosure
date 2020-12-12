#include <stdio.h>
#include <stdlib.h>

#include "closure/closure-proteins.h"
#include "translate_tools.h"

int main(int argc, char **argv) {

  if (argc < 3) {
    printf("syntax: build <seqfile> <HU position>\n");
    return 0;
  }

  dna *b = new dna(argv[1]);

  closure *c = new closure();

  c->initialize_bdna(*b);
  c->initialize_params();
  c->init_proteins();

  c->generate_config1();
  c->generate_config2();
  c->convert_steps();

  FILE *f = fopen(argv[2], "r");

  int n;
  fscanf(f, "%d", &n);

  int pos;

  printf("place proteins\n");

  for (int i = 0; i < n; i++) {
    fscanf(f, "%d", &pos);
    c->place_protein(c, (pos > 0 ? pos : -pos), 0, (pos > 0 ? 0: 1));
    matrix W = identity(4);
    for (int i = 0; i < pos; i++) W = W * calculateW(c->v[i]);
    char s1[40];
    sprintf(s1, "proteins.pdb", i);
    rewrite_pdb("HU-1.pdb", s1, identity(4), W, (i == 0 ? 0 : 1));
  }

  c->print3dna("output.steps");
  printf("%lf\n", c->calculate_writhe());

  fclose(f);

  return 0;

}
