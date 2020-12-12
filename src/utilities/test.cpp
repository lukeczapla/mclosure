#include <stdio.h>
#include <stdlib.h>
#include "dna/dna.h"

#define number_of_bases 4

main() {
  FILE *f = fopen("Fij_AM.dat", "r");

  char temp[10];
  char bp[16][5];
  float Fij[16][6][6];
  char s[200];
  for (int i = 0; i < 16; i++) {
    fgets(s, 199, f);
    sscanf(s, "%s", bp[i]);
    int q = translate(bp[i][0])*number_of_bases+translate(bp[i][1]);
    for (int j = 0; j < 6; j++) {
      fgets(s, 199, f);
      sscanf(s, "%s %f %f %f %f %f %f", temp, &Fij[q][j][0], &Fij[q][j][1], &Fij[q][j][2], &Fij[q][j][3], &Fij[q][j][4], &Fij[q][j][5]);
    }
  }

  printf("#realFC of complete sample - alex morozov\n   CG     CA     TA     AG     GG     AA     GA     AT     AC     GC     ZZ     ZA     ZT     ZG     ZC     AZ     TZ     GZ     CZ\n");

  for (int i = 0; i < 21; i++) {    
    for (int j = 0; j < 10; j++) {
      int q;
      switch (j) {
	case 0: { q = translate('C')*number_of_bases+translate('G'); break; }
        case 1: { q = translate('C')*number_of_bases+translate('A'); break; }
        case 2: { q = translate('T')*number_of_bases+translate('A'); break; }
        case 3: { q = translate('A')*number_of_bases+translate('G'); break; }
        case 4: { q = translate('G')*number_of_bases+translate('G'); break; }
        case 5: { q = translate('A')*number_of_bases+translate('A'); break; }
        case 6: { q = translate('G')*number_of_bases+translate('A'); break; }
        case 7: { q = translate('A')*number_of_bases+translate('T'); break; }
        case 8: { q = translate('A')*number_of_bases+translate('C'); break; }
        case 9: { q = translate('G')*number_of_bases+translate('C'); break; }

      }
      switch (i) {
	case 0: { printf(" %8.3f", Fij[q][2][2]); break; }
        case 1: { printf(" %8.3f", Fij[q][1][1]); break; }
        case 2: { printf(" %8.3f", Fij[q][0][0]); break; }
        case 3: { printf(" %8.3f", Fij[q][3][3]); break; }
        case 4: { printf(" %8.3f", Fij[q][4][4]); break; }
        case 5: { printf(" %8.3f", Fij[q][5][5]); break; }
        case 6: { printf(" %8.3f", Fij[q][2][1]); break; }
        case 7: { printf(" %8.3f", Fij[q][2][0]); break; }
        case 8: { printf(" %8.3f", Fij[q][1][0]); break; }
        case 9: { printf(" %8.3f", Fij[q][3][4]); break; }
        case 10: { printf(" %8.3f", Fij[q][3][5]); break; }
        case 11: { printf(" %8.3f", Fij[q][4][5]); break; }
        case 12: { printf(" %8.3f", Fij[q][2][3]); break; }
        case 13: { printf(" %8.3f", Fij[q][1][3]); break; }
        case 14: { printf(" %8.3f", Fij[q][0][3]); break; }
        case 15: { printf(" %8.3f", Fij[q][2][4]); break; }
        case 16: { printf(" %8.3f", Fij[q][1][4]); break; }
        case 17: { printf(" %8.3f", Fij[q][0][4]); break; }
        case 18: { printf(" %8.3f", Fij[q][2][5]); break; }
        case 19: { printf(" %8.3f", Fij[q][1][5]); break; }
        case 20: { printf(" %8.3f", Fij[q][0][5]); break; }
 
      } 
    }
    printf(" %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }  

}
