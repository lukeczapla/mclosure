
#include "dna.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int translate(char base) {
  switch (toupper(base)) {
  case 'A': return 0;
  case 'T': return 3;
  case 'U': return 3;
  case 'G': return 1;
  case 'C': return 2;
  case 'Z': return 4;
  default:  return -1;
  }
}

char complement(char base, int rna) {
  switch(toupper(base)) {
  case 'A': {
	if (!rna) return 'T';
	return 'U';
  }
  case 'T': return 'A';
  case 'U': return 'A';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'Z': return 'T';
  default: return ' ';
  }
}

char btranslate(int base, int rna) {
  switch (base) {
  case 0: return 'A';
  case 1: return 'G';
  case 2: return 'C';
  case 3: {
	if (!rna) return 'T';
	return 'U';
  }
  case 4: return 'Z';
  default: return ' ';
  }
}

void dna::operator+=(const dna &a) {
  if (nbp < maxseq-a.nbp) {
    strcat(sequence, a.sequence);
    nbp += a.nbp;
    return;
  }
  printf("Could not add any more to the sequence!\n");
  return;
}

void dna::operator--() {
  nbp--;
  sequence[nbp]='\0';
  return;
}

void dna::operator+=(const char &c) {
  if (nbp < maxseq) {
    sequence[nbp++] = c;
    sequence[nbp] = '\0';
    return;
  }
  printf("Could not add any more to the sequence!\n");
}

dna::dna(char *filename) {

  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    printf("Sequence file could not be opened!\n\n");
    exit(0);
  }

  char seq[maxseq];
  fscanf(f, "%s", seq);
  sequence = new char[maxseq];
  strcpy(sequence, seq);
  nbp = strlen(sequence);
  printf("The sequence reads: \n%s\n", seq);
  printf("Contains %d base pairs\n", nbp); 
  fclose(f);
}

dna::dna(char *s, int len) {
  sequence = new char[maxseq];
  strcpy(sequence, s);
  nbp = len;
}

dna::dna() {
  char seq[maxseq];
  printf("Enter the sequence: ");
  scanf("%s", seq);
  sequence = new char[maxseq];
  strcpy(sequence, seq);
  nbp = strlen(sequence);
}

dna::~dna() {
  delete [] sequence;
}

dna::dna(const dna &a) {
  sequence = new char[maxseq];
  strcpy(sequence, a.sequence);
  nbp = strlen(sequence);
  return;
}


