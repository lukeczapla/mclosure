#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

class histogram {

  double bs;
  int nb;
  int *bins[2];

 public:
  histogram(int nboxes, double binsize);
  histogram();
  ~histogram();

  void construct_histogram(int nboxes, double binsize);

  void add_data(double value);
  void add_data(double value, int n);
  void printhistogram(char *filename);
  void printhistogram_pm(char *filename);
  void printhistogram_mm(char *filename);
  void printhistogram_pi(char *filename);
};


class extended_histogram {

  double bs;
  int nb;
  int **bins;
  int max;
  int N;

public:
  extended_histogram(int nboxes, double binsize, int n);
  ~extended_histogram();

  void add_data(int n, double value);

  void printhistogram_pm(char *filename);
  void printhistogram(char *filename);
  void printhistogram_pi(char *filename);
  void printhistogram_pi(char *filename, int b1);


};

#endif
