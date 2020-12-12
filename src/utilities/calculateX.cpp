#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("usage: ./program <radial file>\n");
    return 0;
  }

  char filename[256];

  strcpy(filename, argv[1]);

  FILE *f = fopen(filename, "r");

  double x, p;

  double integral1 = 0.0;
  double normal1 = 0.0;

  double F = 0.1;

  printf("Input force (kT/Angstrom): ");
  scanf("%lf", &F);

  int nb = 0;
  printf("Number of bins: ");
  scanf("%d", &nb);

  for (int i = 0; i < nb; i++) {
    fscanf(f, "%lf %lf", &x, &p);
    if (p != 0.0) {
      integral1 += x*p*exp(F*x)/(x*x);
      normal1 += p*exp(F*x)/(x*x);
    }
  }

  printf("Integral x = %lf Angstroms\n", integral1/normal1);

}
