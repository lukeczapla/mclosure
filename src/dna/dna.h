#ifndef _DNA_H
#define _DNA_H

#include <stdio.h>

#define number_of_bases 5
#define maxseq 32768

int translate(char base);
char btranslate(int base, int rna = 0);
char complement(char base, int rna = 0);
class dna {
 
public:
  int nbp;
  char *sequence;
  dna(char *filename);
  dna();
  dna(const dna &a);
  dna(char *s, int len);
  ~dna();

  void operator--();
  void operator+=(const dna &a);
  void operator+=(const char &c);

};

#endif
