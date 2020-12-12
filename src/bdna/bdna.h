#ifndef _BDNA_H
#define _BDNA_H

#include "../dna/dna.h"
#include "../matrix/matrix.h"

matrix calculateW(const matrix &tp);
matrix calculateM(const matrix &tp);
matrix calculatetp(const matrix &W);
matrix invert(matrix a);

void readFparameters(char *filename, matrix F[]);
void readIparameters(char *filename, matrix I[]);


class bdna {

public:

  matrix I[number_of_bases*number_of_bases];
  matrix F[number_of_bases*number_of_bases];
  matrix *v;
  matrix *phosphate;

  FLOAT beta;
  FLOAT ion;  

  char *seq;

  int *bpstep;
  int nsteps;

  bdna();
  bdna(const bdna &b);
  bdna(int ns);
  bdna(int ns, int q);
  ~bdna();

  int overlap();
  int overlap(int a, int b);

  void add_phosphates();
  void delete_phosphates();

  void set_beta(FLOAT be);

  void resize(int ns);

  void initialize_bdna(dna D);

  void set_straight();

  void addbase(char c);

  void stress_free_state(int write);

  FLOAT calculateR();

  void print_end();
  void print_energy();
  FLOAT calculate_energy(matrix &W);
  FLOAT calculate_elenergy();
  FLOAT calculate_elenergy(int a);
  FLOAT calculate_elenergy(int a, int b);
  FLOAT calculate_elenergy(matrix tp0, int m);
  FLOAT calculate_esenergy(matrix &W);
  FLOAT calculate_penergy();

  FLOAT calculate_link3();
  FLOAT calculate_link();
  FLOAT calculate_link_new();
  FLOAT calculate_writhe();
  FLOAT calculate_writhe_closed();
  FLOAT calculate_twist();
//  FLOAT calculate_twist2();
  FLOAT calculate_twist_open();
  FLOAT calculate_twist_closed();
  FLOAT calculate_twist_print();
  int calculate_twist0();

  void printconfig(char *filename);
  void readconfig(char *filename);

  void printpdb(char *filename);
  void print3dna(char *filename);

  bdna &operator=(const bdna &b);

};

#endif
