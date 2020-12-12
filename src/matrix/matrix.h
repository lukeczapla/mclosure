
#ifndef _MATRIX_H
#define _MATRIX_H

#include "../defines.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class matrix {

  FLOAT *v;

public:

  unsigned char m,n;
  
  matrix();
  matrix(int a, int b);
  matrix(const matrix &M);
  ~matrix();

  inline void setdata(FLOAT *d) {
    memcpy(v, d, m*n*sizeof(FLOAT));
  };

  inline void getdata(FLOAT *d) {
    memcpy(d, v, m*n*sizeof(FLOAT));
  };

  inline void setsize(char a, char b) {
    m = a;
    n = b;
    delete [] v;
    v = new FLOAT[a*b];
  };
 
  inline void setv(char a, char b, FLOAT d) {
    v[(a-1)*n+(b-1)] = d;
  };

  inline FLOAT operator()(char a, char b) const {
    return v[(a-1)*n+(b-1)];
  };

  inline matrix &operator=(const matrix &M) {
    if (this != &M) {
      delete [] v;
      v = new FLOAT[M.m*M.n];
      m = M.m;
      n = M.n;
      //for (int i = 0; i < m*n; i++) v[i] = M.v[i];
      memcpy(v, M.v, m*n*sizeof(FLOAT));
    }
    return *this;
  };

  friend void readmatrix(FILE *f, matrix M);
  friend void writematrix(FILE *f, matrix M);
  friend int equals(matrix a, matrix b);

  friend matrix operator*(const matrix &a, const matrix &b);
  friend matrix operator*(const FLOAT &d, const matrix &a);
  friend matrix operator+(const matrix &a, const matrix &b);
  friend matrix operator-(const matrix &a, const matrix &b);

  friend matrix identity(int m);
  friend matrix jeigen(matrix a, FLOAT *d);
};

matrix identity(int m);
matrix badmatrix();

#endif







