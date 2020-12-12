#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix/matrix.h"
#include "bdna/bdna.h"

void read_reference(char *fname, matrix *W) {
  FILE *f = fopen(fname, "r");
  float x, y, z;
  fscanf(f, "%f %f %f", &x, &y, &z);
  W->setv(1,4,x);
  W->setv(2,4,y);
  W->setv(3,4,z);
  W->setv(4,4,1.0);

  fscanf(f, "%f %f %f", &x, &y, &z);
  W->setv(1,1,x);
  W->setv(2,1,y);
  W->setv(3,1,z);
  W->setv(4,1,0.0);
  
  fscanf(f, "%f %f %f", &x, &y, &z);
  W->setv(1,2,x);
  W->setv(2,2,y);
  W->setv(3,2,z);
  W->setv(4,2,0.0);

  fscanf(f, "%f %f %f", &x, &y, &z);
  W->setv(1,3,x);
  W->setv(2,3,y);
  W->setv(3,3,z);
  W->setv(4,3,0.0);

  writematrix(stdout, *W);

  fclose(f);
}


float det3(matrix A) { 
  return A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(3,3)*A(2,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
}


void translate_position(matrix W, float &x, float &y, float &z) {
  float x1 = x;  float y1 = y;  float z1 = z;
//  x1 -= W(1,4);  y1 -= W(2,4);  z1 -= W(3,4);

  printf("%f\n", sqrt((x-W(1,4))*(x-W(1,4))+(y-W(2,4))*(y-W(2,4))+(z-W(3,4))*(z-W(3,4))));

  matrix Q = W;

//  Q.setv(1,4,0.0);
//  Q.setv(2,4,0.0);
//  Q.setv(3,4,0.0);

  matrix P = invert(Q);

  x = x1*P(1,1)+y1*P(1,2)+z1*P(1,3)+P(1,4);
  y = x1*P(2,1)+y1*P(2,2)+z1*P(2,3)+P(2,4);
  z = x1*P(3,1)+y1*P(3,2)+z1*P(3,3)+P(3,4);

/*  printf("%f %f %f %f %f %f\n", x, y, z, x1, y1, z1);

  printf("%f\n", sqrt(x1*x1+y1*y1+z1*z1));

  matrix W2;

  matrix Wdenom(3,3);
  for (int i = 1; i <= 3; i++)
    for (int j = 1; j <= 3; j++)
      Wdenom.setv(i, j, W(i,j));

  float denom = det3(Wdenom);

  W2 = Wdenom;
  W2.setv(1,1,x1);
  W2.setv(2,1,y1);
  W2.setv(3,1,z1);

  x = -det3(W2)/denom;

  W2 = Wdenom;
  W2.setv(1,2,x1);
  W2.setv(2,2,y1);
  W2.setv(3,2,z1);

  y = -det3(W2)/denom;

  W2 = Wdenom;
  W2.setv(1,3,x1);
  W2.setv(2,3,y1);
  W2.setv(3,3,z1);

  z = -det3(W2)/denom;

  printf("%f %f %f\n", x, y, z);

  printf("%f\n", sqrt(x*x+y*y+z*z));
*/
}

void rewrite_pdb(char *in, char *out, matrix W, matrix new1) {
    char atom[8];
    char anum[7];
    char at[7];
    char res[6];
    char chain[5]; 
    char resn[6];
    char none[6];
    float x, y, z;
    FILE *f = fopen(in, "r");
    FILE *f2 = fopen(out, "w");
    char s[257];
    char *p;
    do {
      p = fgets(s, 256, f);
      strcpy(atom, "      ");
      strcpy(chain, "  ");
      strcpy(res, "    ");
      strcpy(resn, "    ");
      strcpy(at, "     ");
      strcpy(anum, "     ");
      strcpy(none, "    ");
      sscanf(s, "%6c%5c%5c%4c%2c%4c%4c%8f%8f%8f", atom, anum, at, res, chain, resn, none, &x, &y, &z);
      if (strcmp(atom, "ATOM  ") == 0) {
      translate_position(W, x, y, z);
      float xn = new1(1,4)+x*new1(1,1)+y*new1(1,2)+z*new1(1,3);
      float yn = new1(2,4)+x*new1(2,1)+y*new1(2,2)+z*new1(2,3);
      float zn = new1(3,4)+x*new1(3,1)+y*new1(3,2)+z*new1(3,3);

     // printf("%f %f %f", x, y, z);
      fprintf(f2, "%6s%5s%5s%4s%2s%4s%4s%8.3f%8.3f%8.3f  1.00 15.00\n", atom, anum, at, res, chain, resn, none, xn, yn, zn);
      }
      
    } while ((strcmp(atom, "ATOM  ") == 0) && (p != NULL));
    fclose(f);
    fclose(f2);
}



main(int argc, char **argv) {
  if (argc != 4) {
    printf("usage: translate_protein [step.dat] [input.pdb] [output.pdb]\n");
    exit(0);
  }

  matrix W(4,4);
  read_reference(argv[1], &W);


  matrix n = identity(4);

  writematrix(stdout, W);

  rewrite_pdb(argv[2], argv[3], W, n); 

}

