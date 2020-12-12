#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "histogram.h"
#include "../file.h"

histogram::histogram(int nboxes, double binsize) {
  nb = nboxes;
  bs = binsize;
  bins[0] = new int[nboxes];
  bins[1] = new int[nboxes];
  for (int i = 0; i < nboxes; i++) {
    bins[0][i] = 0;
    bins[1][i] = 0;
  } 
}

histogram::histogram() {
}

histogram::~histogram() {
  delete [] bins[0];
  delete [] bins[1];
}

void histogram::construct_histogram(int nboxes, double binsize) {
  nb = nboxes;
  bs = binsize;
  bins[0] = new int[nboxes];
  bins[1] = new int[nboxes];
  for (int i = 0; i < nboxes; i++) {
    bins[0][i] = 0;
    bins[1][i] = 0;
  }
}


void histogram::add_data(double value) {
  int bin = abs((int)(value/bs));
  if (bin < nb) {
    if (value >= 0.0) bins[0][bin]++;
    else bins[1][bin]++; 
  }
}

void histogram::add_data(double value, int n) {
  int bin = abs((int)(value/bs));
  if (bin < nb) {
    if (value >= 0.0) bins[0][bin]=bins[0][bin]+n;
	else bins[1][bin]=bins[1][bin]+n;
  }
}

void histogram::printhistogram_mm(char *filename) {
  FILE *f = openfwrite(filename);
  if (f == NULL) {
    printf("Could not output histogram file %s", filename);
    return;
  }
  int i = nb-1;
  int j = 0;
  for (; i >= j; i--) {
    fprintf(f, "%lf %d\n", -(double)i*bs-bs/2.0, bins[1][i]);
  }
  i = 0;
  j = nb;
  for (; i < j; i++) {
    fprintf(f, "%lf %d\n", (double)i*bs+bs/2.0, bins[0][i]);
  }
  fclose(f);

}
void histogram::printhistogram_pm(char *filename) {
  FILE *f = openfwrite(filename);
  if (f == NULL) {
    printf("Could not output histogram file %s", filename);
    return;
  }
    int i = nb-1;
    int j = 0;
//  for (; i >= j; i--) {
 //   fprintf(f, "%lf %d\n", -(double)i*bs, bins[1][i]);
//  }
  i = 0;
  j = nb;
  for (; i < j; i++) {
    fprintf(f, "%lf %d\n", (double)i*bs+bs/2.0, bins[0][i]);
  }
  fclose(f);

}

void histogram::printhistogram_pi(char *filename) {
  FILE *f = openfwrite(filename);
  if (f == NULL) {
    printf("Could not output histogram file %s", filename);
    return;
  }
  long int NT = 0;
  for (int i = 0 ; i < nb; i++) NT += bins[0][i];
  for (int i = 0; i < nb; i++) {
    fprintf(f, "%f %2.7f\n", bs*i, (float)bins[0][i]/(float)NT);
  }
  fclose(f);
}

void histogram::printhistogram(char *filename) {
  FILE *f = openfwrite(filename);
  if (f == NULL) {
    printf("Could not output histogram file %s", filename);
    return;
  }
  int max = 0;
  int i;
  for (i = 0; i < nb; i++) if (bins[0][i] != 0) max = i;
  for (i = 0; i <= max; i++) {
    fprintf(f, "%d %d\n", i, bins[0][i]);
  }
  fclose(f);
}


extended_histogram::extended_histogram(int nboxes, double binsize, int n) {
  nb = nboxes;
  bs = binsize;
  N = n;

  bins = new int*[N];
  for (int i = 0; i < N; i++) {
    bins[i] = new int[nb];
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < nb; j++) {
      bins[i][j] = 0;
    }
  }

  max = 0;
 
}

extended_histogram::~extended_histogram() {
  if (N == 0) return;
  for (int i = 0; i < N; i++) {
	delete [] bins[i];
  }
  delete [] bins;
}

void extended_histogram::add_data(int n, double value) {
  bins[n][(int)(value/bs)]++;
  if (n > max) max = n;
}

void extended_histogram::printhistogram_pi(char *filename) {
  FILE *f = openfwrite(filename);
  for (int i = 0; i < nb; i++) {
	fprintf(f, "%d", (int)(bs*i));
    for (int j = 0; j <= max; j++) {
	  fprintf(f, " %d", bins[j][i]);
    }
	fprintf(f, "\n");
  }
  fclose(f);
}


void extended_histogram::printhistogram_pi(char *filename, int b1) {
  FILE *f = openfwrite(filename);
  for (int i = 0; i < nb; i++) {
	fprintf(f, "%d", (int)(bs*i));
    for (int j = b1; j <= max; j++) {
	  fprintf(f, " %d", bins[j][i]);
    }
	fprintf(f, "\n");
  }
  fclose(f);
}


void extended_histogram::printhistogram(char *filename) {
  FILE *f = openfwrite(filename);
  int tmin = 100000;
  int tmax = 0;
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j <= max; j++) {
	  if ((tmin == 100000) && (bins[j][i] != 0)) tmin = i;
      if (bins[j][i] != 0) tmax = i;
    }
  }
  for (int i = tmin; i <= tmax; i++) {
	fprintf(f, "%d", (int)(bs*i));
    for (int j = 0; j <= max; j++) {
	  fprintf(f, " %d", bins[j][i]);
    }
	fprintf(f, "\n");
  }
  fclose(f);
}



void extended_histogram::printhistogram_pm(char *filename) {
  FILE *f = openfwrite(filename);
  for (int i = 0; i < nb; i++) {
	fprintf(f, "%f", bs*i+bs/2.0);
    for (int j = 0; j <= max; j++) {
	  fprintf(f, " %d", bins[j][i]);
    }
	fprintf(f, "\n");
  }
  fclose(f);
}

