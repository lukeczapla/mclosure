
#define vof(A,B,N) (A-1)*N+(B-1)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix.h"
#include "../defines.h"

matrix::matrix() {
  v = 0;
  m = 0;
  n = 0;
}

matrix::matrix(int a, int b) {
  v = new FLOAT[a*b];
  m = a;
  n = b;
}

matrix::~matrix() {
  delete [] v;
}

matrix badmatrix() {
  return matrix(0,0);
}

matrix identity(int m) {
  matrix P(m,m);
  int index = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      if (i == j) P.v[index] = 1.0;
        else P.v[index] = 0.0;
      index++;
    }
  }
  return P;
}

matrix operator*(const FLOAT &d, const matrix &a) {
  matrix P(a.m, a.n);
  for (int i = 0; i < a.m*a.n; i++) P.v[i] = a.v[i]*d;
  return P;
}

matrix operator*(const matrix &a, const matrix &b) {
  if (a.n == b.m) {
    matrix P(a.m, b.n);
    int index = 0;
    int indexi, indexj;
    for (int i = 1; i <= a.m; i++) {
      for (int j = 1; j <= b.n; j++) {
	P.v[index]=0.0;
        indexi = (i-1)*a.n;
        indexj = (j-1);
	for (int k = 1; k <= a.n; k++) {
	  P.v[index] += a.v[indexi]*b.v[indexj];
	  indexi += 1;
          indexj += b.n;
	}
        index++;
      }
    }
    return P;
  } 
  printf("wrong dimensions in matrix multiplication!\n");
  return badmatrix();
}
/*
void matrix::operator*=(const matrix &b) {

  if (n == b.m) {
    matrix P;
    P = *this;
    int index = 0;
    int indexi, indexj;
    for (int i = 1; i <= m; i++) {
      for (int j = 1; j <= b.n; j++) {
	v[index]=0;
        indexi = (i-1)*n;
        indexj = (j-1);
	for (int k = 1; k <= n; k++) {
	  v[index] += P.v[indexi]*b.v[indexj];
	  indexi += 1;
          indexj += b.n;
	}
        index++;
      }
    }
    return;
  } 
  printf("wrong dimensions in matrix multiplication!\n");
  return;
}
*/

matrix operator+(const matrix &a, const matrix &b) {
  if ((a.m == b.m) && (a.n == b.n)) {
    matrix S(a.m, a.n);
    int index = 0;
    for (int i = 1; i <= a.m; i++) {
      for (int j = 1; j <= a.n; j++) {
	S.v[index] = a.v[index] + b.v[index];
        index++;
      }
    }
    return S;
  }
  printf("Bad matrix addition!\n");
  return badmatrix();
}

matrix operator-(const matrix &a, const matrix &b) {
  if ((a.m == b.m) && (a.n == b.n)) {
    matrix S(a.m, a.n);
    int index = 0;
    for (int i = 1; i <= a.m; i++) {
      for (int j = 1; j <= a.n; j++) {
	S.v[index] = a.v[index] - b.v[index];
        index++;
      }
    }
    return S;
  }
  printf("Bad matrix sub-addition!\n");
  return badmatrix();
}

int equals(matrix a, matrix b) {
  if ((a.m != b.m) || (a.n != b.n)) return 0;
  for (int i = 0; i < a.m*b.n; i++) {
    if (a.v[i] != b.v[i]) return 0;
  }
  return 1;
}

void readmatrix(FILE *f, matrix M) {
  for (int i = 0; i < M.m*M.n; i++) {
      fscanf(f, FLOAT_ESTR, &M.v[i]);
  }
}

void writematrix(FILE *f, matrix M) {
  for (int i = 1; i <= M.m; i++) {
    for (int j = 1; j <= M.n; j++) {
      fprintf(f, " ");
      fprintf(f, FLOAT_ESTRF, M.v[vof(i,j,M.n)]);
    }
    fprintf(f, "\n");
  }
}


matrix jeigen(matrix a, FLOAT *d) {

  if (a.m != a.n) {
    printf("wrong dimensions in diagonalization\n");
    exit(0);
  }

  int ip, iq, i, j;
  FLOAT tresh, theta, tau, t, sm, s, h, g, c;

  matrix x;
  matrix v = identity(a.m);

  x = a;

  FLOAT *b = new FLOAT[a.n];
  FLOAT *z = new FLOAT[a.n];

  tresh = 0.0;

  for (ip = 0; ip < a.n; ip++) {
    b[ip]=(d[ip]=x.v[ip*a.n+ip]);
    z[ip]=0.0;
  }

  for (i = 0; i < 500; i++) {
    //  x.writematrix(stdout);
    sm = 0.0;
    for(ip=0; ip < a.n-1; ip++)
      for (iq=ip+1; iq < a.n; iq++) sm += fabs(x.v[ip*a.n+iq]);
    if (sm == 0.0) {
       delete [] b;
       delete [] z;
       //v.writematrix(stdout);
       return v;
    }

    for (ip = 0; ip < a.n-1; ip++) {
      for (iq = ip+1; iq < a.n; iq++) {
        g = 100.0*fabs(x.v[ip*a.n+iq]);
        if (fabs(x.v[ip*a.n+iq]) > tresh) {

          h = d[iq]-d[ip];

          if ((FLOAT)(fabs(h)+g) == (FLOAT)fabs(h)) {
	    if (h != 0.0)
  	      t = (x.v[ip*a.n+iq])/h;
            else t = 0.0;
	  } else {
            if (x.v[ip*a.n+iq] != 0.0) 
   	      theta = 0.5*h/x.v[ip*a.n+iq];
	    else theta = 0.0;
            t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
	  }

          c = 1.0/sqrt(1.0+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*x.v[ip*a.n+iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;

          x.v[ip*a.n+iq] = 0.0;
	  for (j = 0; j <= ip-1; j++) {
            g=x.v[j*a.n+ip];
            h=x.v[j*a.n+iq];
            x.v[j*a.n+ip]=g-s*(h+g*tau);
            x.v[j*a.n+iq]=h+s*(g-h*tau);
	  }
          for (j = ip+1; j <= iq-1; j++) {
            g=x.v[ip*a.n+j];
            h=x.v[j*a.n+iq];
            x.v[ip*a.n+j]=g-s*(h+g*tau);
            x.v[j*a.n+iq]=h+s*(g-h*tau);	   
          }
          for (j = iq+1; j < a.n; j++) {
            g=x.v[ip*a.n+j];
            h=x.v[iq*a.n+j];
            x.v[ip*a.n+j]=g-s*(h+g*tau);
            x.v[iq*a.n+j]=h+s*(g-h*tau);
	  }
          for (j = 0; j < a.n; j++) {
	    g=v.v[j*a.n+ip];
            h=v.v[j*a.n+iq];
            v.v[j*a.n+ip]=g-s*(h+g*tau);
            v.v[j*a.n+iq]=h+s*(g-h*tau);
	  }

        }
      }
    }
      for (ip=0; ip < a.n; ip++) {
	b[ip] += z[ip];
        d[ip] = b[ip];
	z[ip] = 0.0;
      }
  }

  printf("no luck with the matrix routines!\n");
  writematrix(stdout, x);

  return badmatrix();

}


matrix::matrix(const matrix &M) {
  v = new FLOAT[M.m*M.n];
  m = M.m;
  n = M.n;
  memcpy(v, M.v, m*n*sizeof(FLOAT));
}

