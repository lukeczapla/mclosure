

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bdna.h"
#include "../defines.h"

#define ELECTROSTATICS
#define DNA_DNA_DIST 20.0
// #define SELF_CONTACT
  const FLOAT eps = 77.4;
  const FLOAT e = 0.48*1.6021773e-19;
  const FLOAT eps0 = 8.854188e-22;
  const FLOAT kT = 1.380658e-23*300.0;
  const FLOAT coeff = e*e/(4.0*M_PI*eps0*eps*kT);


matrix invert(matrix a) {

  matrix temp(4,4);

  FLOAT A = a(1,1);
  FLOAT B = a(1,2);
  FLOAT C = a(1,3);
  FLOAT D = a(1,4);
  FLOAT E = a(2,1);
  FLOAT F = a(2,2);
  FLOAT G = a(2,3);
  FLOAT H = a(2,4);
  FLOAT I = a(3,1);
  FLOAT J = a(3,2);
  FLOAT K = a(3,3);
  FLOAT L = a(3,4);
  FLOAT denom = -C*F*I+B*G*I+C*E*J-A*G*J-B*E*K+A*F*K;

  temp.setv(1,1, (-G*J+F*K)/denom);
  temp.setv(1,2, (C*J-B*K)/denom);
  temp.setv(1,3, (-C*F+B*G)/denom);
  temp.setv(1,4, (D*G*J-C*H*J-D*F*K+B*H*K+C*F*L-B*G*L)/denom);
  temp.setv(2,1, (G*I-E*K)/denom);

  temp.setv(2,2, (-C*I+A*K)/denom);
  temp.setv(2,3, (C*E-A*G)/denom);
  temp.setv(2,4, (-D*G*I+C*H*I+D*E*K-A*H*K-C*E*L+A*G*L)/denom);
  temp.setv(3,1, (-F*I+E*J)/denom);
  temp.setv(3,2, (B*I-A*J)/denom);
  temp.setv(3,3, (-B*E+A*F)/denom);
  temp.setv(3,4, (D*F*I-B*H*I-D*E*J+A*H*J+B*E*L-A*F*L)/denom);
  temp.setv(4,1, 0);
  temp.setv(4,2, 0);
  temp.setv(4,3, 0);
  temp.setv(4,4, 1);

  return temp;

}


void crossr(FLOAT &zx, FLOAT &zy, FLOAT &zz, FLOAT xx, FLOAT xy, FLOAT xz,
	    FLOAT yx, FLOAT yy, FLOAT yz) {
  zx = xy*yz-xz*yy;
  zy = xz*yx-xx*yz;
  zz = xx*yy-xy*yx;
}


void normalize(FLOAT &x, FLOAT &y, FLOAT &z) {
  FLOAT norm = sqrt(x*x+y*y+z*z);
  x /= norm;
  y /= norm;
  z /= norm;
}


int vequal(FLOAT x, FLOAT y, FLOAT z, FLOAT x2, FLOAT y2, FLOAT z2) {
  return ((x == x2) && (y == y2) && (z == z2));
}


FLOAT dih(FLOAT tx, FLOAT ty, FLOAT tz, FLOAT ux, FLOAT uy, FLOAT uz,
	  FLOAT vx, FLOAT vy, FLOAT vz) {
#define EPSdih 1e-14
  FLOAT tux, tuy, tuz, uvx, uvy, uvz, tvx, tvy, tvz;
  crossr(tux, tuy, tuz, tx, ty, tz, ux, uy, uz);
  crossr(uvx, uvy, uvz, ux, uy, uz, vx, vy, vz);
  crossr(tvx, tvy, tvz, tx, ty, tz, vx, vy, vz);
  if ((tux*tux+tuy*tuy+tuz*tuz < EPSdih) || (uvx*uvx+uvy*uvy+uvz*uvz < EPSdih)) return 0.0;
  return atan2(sqrt(ux*ux+uy*uy+uz*uz)*(tvx*ux+tvy*uy+tvz*uz), 
               tux*uvx+tuy*uvy+tuz*uvz);
}


FLOAT correct(FLOAT mur, FLOAT wrr) {
  if (wrr-mur > M_PI) return wrr-2*M_PI;
  if (wrr-mur <= -M_PI) return wrr+2*M_PI;
  return wrr;
}
  

matrix calculateW(const matrix &tp) {

  matrix M(4,4);

  FLOAT gamma, phi, omega,
         sp, cp, sm, cm, sg, cg,
         t1, t2, t3;

  t1 = tp(1,1)*M_PI/180.0;
  t2 = tp(2,1)*M_PI/180.0;
  t3 = tp(3,1)*M_PI/180.0;

  gamma = sqrt(t1*t1+t2*t2);
  phi = atan2(t1,t2);
  omega = t3;

  sp = sin(omega/2+phi); cp = cos(omega/2+phi); sm = sin(omega/2-phi);
  cm = cos(omega/2-phi); sg = sin(gamma); cg = cos(gamma);

  M.setv(1,1, cm*cg*cp-sm*sp);
  M.setv(1,2, -cm*cg*sp-sm*cp);
  M.setv(1,3, cm*sg);
  M.setv(2,1, sm*cg*cp+cm*sp);
  M.setv(2,2, -sm*cg*sp+cm*cp);
  M.setv(2,3, sm*sg);
  M.setv(3,1, -sg*cp);
  M.setv(3,2, sg*sp);
  M.setv(3,3, cg);
  M.setv(4,1, 0);
  M.setv(4,2, 0);
  M.setv(4,3, 0);
  M.setv(4,4, 1);

  sp = sin(phi); cp = cos(phi); sg = sin(gamma/2); cg = cos(gamma/2);

  M.setv(1,4, tp(4,1)*(cm*cg*cp-sm*sp) + tp(5,1)*(-cm*cg*sp-sm*cp) + tp(6,1)*(cm*sg));
  M.setv(2,4, tp(4,1)*(sm*cg*cp+cm*sp) + tp(5,1)*(-sm*cg*sp+cm*cp) + tp(6,1)*(sm*sg));
  M.setv(3,4, tp(4,1)*(-sg*cp) + tp(5,1)*(sg*sp) + tp(6,1)*(cg));

  return M;

}


matrix calculateM(const matrix &tp) {

  matrix M(4,4);

  FLOAT gamma, phi, omega,
    sp, cp, sm, cm, sg, cg,
         t1, t2, t3;


  t1 = tp(1,1)*M_PI/180.0;
  t2 = tp(2,1)*M_PI/180.0;
  t3 = tp(3,1)*M_PI/180.0;

  gamma = 0.5*sqrt(t1*t1+t2*t2);
  phi = atan2(t1,t2);
  omega = t3;

  sp = sin(phi); cp = cos(phi); sm = sin(omega/2.0-phi);
  cm = cos(omega/2.0-phi); sg = sin(gamma); cg = cos(gamma);

  M.setv(1,1, cm*cg*cp-sm*sp);
  M.setv(1,2, -cm*cg*sp-sm*cp);
  M.setv(1,3, cm*sg);
  M.setv(2,1, sm*cg*cp+cm*sp);
  M.setv(2,2, -sm*cg*sp+cm*cp);
  M.setv(2,3, sm*sg);
  M.setv(3,1, -sg*cp);
  M.setv(3,2, sg*sp);
  M.setv(3,3, cg);
  M.setv(4,1, 0.0);
  M.setv(4,2, 0.0);
  M.setv(4,3, 0.0);
  M.setv(4,4, 1.0);

  /*  M.setv(1,4, 0.5*(tp(4,1)*(cm*cg*cp-sm*sp) + tp(5,1)*(-cm*cg*sp-sm*cp) + tp(6,1)*(cm*sg)));
  M.setv(2,4, 0.5*(tp(4,1)*(sm*cg*cp+cm*sp) + tp(5,1)*(-sm*cg*sp+cm*cp) + tp(6,1)*(sm*sg)));
  M.setv(3,4, 0.5*(tp(4,1)*(-sg*cp) + tp(5,1)*(sg*sp) + tp(6,1)*(cg)));
  */

  M.setv(1,4, 0.0);
  M.setv(2,4, 0.0);
  M.setv(3,4, 0.0);

  return M;

}


matrix calculatetp(const matrix &W) {

  matrix M(6,1);

  FLOAT cosgamma, gamma, phi, omega, sgcp, omega2_minus_phi,
         sm, cm, sp, cp, sg, cg;

  cosgamma = W(3,3);
  if (cosgamma > 1.0) cosgamma = 1.0;
  else if (cosgamma <= -1.0) cosgamma = -1.0;

  gamma = acos(cosgamma);

  sgcp = W(2,2)*W(1,3)-W(1,2)*W(2,3);

  if (gamma == 0.0) omega = -atan2(W(1,2),W(2,2));
  else omega = atan2((W(3,2)*W(1,3)+sgcp*W(2,3)),(sgcp*W(1,3)-W(3,2)*W(2,3)));

  omega2_minus_phi = atan2(W(2,3),W(1,3));

  phi = omega/2.0 - omega2_minus_phi;

  M.setv(1,1, gamma*sin(phi)*180.0/M_PI);
  M.setv(2,1, gamma*cos(phi)*180.0/M_PI);
  M.setv(3,1, omega*180.0/M_PI);

  sm = sin(omega/2.0-phi);  
  cm = cos(omega/2.0-phi);
  sp = sin(phi);
  cp = cos(phi);
  sg = sin(gamma/2.0);
  cg = cos(gamma/2.0);

  M.setv(4,1, (cm*cg*cp-sm*sp)*W(1,4)+(sm*cg*cp+cm*sp)*W(2,4)-sg*cp*W(3,4));
  M.setv(5,1, (-cm*cg*sp-sm*cp)*W(1,4)+(-sm*cg*sp+cm*cp)*W(2,4)+sg*sp*W(3,4));
  M.setv(6,1, (cm*sg)*W(1,4)+(sm*sg)*W(2,4)+cg*W(3,4));
 
  return M;

}



int locate_in_list(int l[], int n, int idx) {
  for (int i = 0; i < idx; i++)
    if (l[i] == n) return i;
  return -1;
}



void readFparameters(char *filename, matrix *F) {

  char *s;
  char *bstep;
  FLOAT d;
  int list[number_of_bases*number_of_bases];
  char *p;
  int idx,m,n;

  char pfilename[256];

  char *path = getenv("MDNA_PATH");

  if (path != NULL) {
    strcpy(pfilename, path);
    if (pfilename[strlen(pfilename)-1] != '/') strcat(pfilename, "/");
    strcat(pfilename, filename); 
  } else {
    strcpy(pfilename, filename);
  }

  FILE *f = fopen(pfilename, "r");
  
  s = new char[2048];

  if (f == NULL) {
    printf("File not found %s\n", filename);
    exit(0);
  }

  do {
    p = fgets(s, 2048, f);
    if (p != NULL) bstep = strtok(s, " ");
  } while ((p != NULL) && ((bstep[0] == '#') || (strlen(s) == 1)));

  if (p == NULL) {
    printf("Error in file format: %s\n", filename);
    exit(0);
  }

  idx = 0;

  do {
    list[idx] = translate(bstep[0])*number_of_bases+translate(bstep[1]);
    if (list[idx] < 0) {
      printf("Error in base pair step format: file %s %c %c", filename,
         bstep[0], bstep[1]);
      exit(0);
    }
    idx++;
  } while (((bstep = strtok(NULL, " ")) != NULL) && (bstep[0] != '#'));

  int a = 0;

  do {
    char *par;
    if (((p=fgets(s, 2048, f)) != NULL) && (sscanf(s, FLOAT_ESTR, &d) == 1)) {
      int q = 0;
      par = strtok(s, " ");
      sscanf(par, FLOAT_ESTR, &d);
      do {
        switch(a) {
	case 0: m=1;n=1;break;
	case 1: m=2;n=2;break;
	case 2: m=3;n=3;break;
	case 3: m=4;n=4;break;
	case 4: m=5;n=5;break;
	case 5: m=6;n=6;break;
        case 6: m=1,n=2;break;
        case 7: m=1,n=3;break;
        case 8: m=2,n=3;break;
        case 9: m=4,n=5;break;
        case 10: m=4,n=6;break;
        case 11: m=5,n=6;break;
        case 12: m=1,n=4;break;
        case 13: m=2,n=4;break;
        case 14: m=3,n=4;break;
        case 15: m=1,n=5;break;
        case 16: m=2,n=5;break;
        case 17: m=3,n=5;break;
        case 18: m=1,n=6;break;
        case 19: m=2,n=6;break;
        case 20: m=3,n=6;break;
        default: {
	  printf("Too many parameters in force constant file!");
	  exit(0);
	}
        }
	F[list[q]].setv(m,n,d); F[list[q]].setv(n,m,d);
        q++;
        par = strtok(NULL, " ");
      } while ((q < idx) && (sscanf(par, FLOAT_ESTR,&d) == 1));
      if (q != idx) {
	printf("not enough specified parameters in %s %d\n", filename, q);
        exit(0);
      }
      a++;
    }

  } while ((a < 21) && (p != NULL));

  if (p == NULL) {
    printf("not enough specified couplings in %s\n", filename);
    exit(0);
  }

  for (int i = 0; i < number_of_bases; i++) {
    for (int j = 0; j < number_of_bases; j++) 
      if (locate_in_list(list, i*number_of_bases+j, idx) == -1)
        if (locate_in_list(list, (number_of_bases-2-j)*number_of_bases+(number_of_bases-2-i), idx) == -1) {
          printf("Incomplete list of parameters in %s %d %d\n", filename, i, j);
	  exit(0);
        } else {
          F[i*number_of_bases+j] = F[(number_of_bases-2-j)*number_of_bases+(number_of_bases-2-i)];
	  F[i*number_of_bases+j].setv(2,1,-F[i*number_of_bases+j](2,1));
	  F[i*number_of_bases+j].setv(1,2,-F[i*number_of_bases+j](1,2));
	  F[i*number_of_bases+j].setv(1,3,-F[i*number_of_bases+j](1,3));
	  F[i*number_of_bases+j].setv(3,1,-F[i*number_of_bases+j](3,1));
	  F[i*number_of_bases+j].setv(5,1,-F[i*number_of_bases+j](5,1));
	  F[i*number_of_bases+j].setv(6,1,-F[i*number_of_bases+j](6,1));
	  F[i*number_of_bases+j].setv(4,2,-F[i*number_of_bases+j](4,2));
	  F[i*number_of_bases+j].setv(4,3,-F[i*number_of_bases+j](4,3));
	  F[i*number_of_bases+j].setv(1,5,-F[i*number_of_bases+j](1,5));
	  F[i*number_of_bases+j].setv(1,6,-F[i*number_of_bases+j](1,6));
	  F[i*number_of_bases+j].setv(2,4,-F[i*number_of_bases+j](2,4));
	  F[i*number_of_bases+j].setv(3,4,-F[i*number_of_bases+j](3,4));
	  F[i*number_of_bases+j].setv(5,4,-F[i*number_of_bases+j](5,4));
	  F[i*number_of_bases+j].setv(6,4,-F[i*number_of_bases+j](6,4));
	  F[i*number_of_bases+j].setv(4,5,-F[i*number_of_bases+j](4,5));
	  F[i*number_of_bases+j].setv(4,6,-F[i*number_of_bases+j](4,6));
	
	}
  }
  fclose(f);
}



void readIparameters(char *filename, matrix I[]) {
  char *s;
  char *bstep;
  FLOAT d;
  int list[number_of_bases*number_of_bases];
  char *p;
  int idx;
  char pfilename[256];

  char *path = getenv("MDNA_PATH");

  if (path != NULL) {
    strcpy(pfilename, path);
    if (pfilename[strlen(pfilename)-1] != '/') strcat(pfilename, "/");
    strcat(pfilename, filename); 
  } else {
    strcpy(pfilename, filename);
  }

  FILE *f = fopen(pfilename, "r+");

  s = new char[2048];

  if (f == NULL) {
    printf("File not found %s\n", filename);
    exit(0);
  }

  do {
    p = fgets(s, 2048, f);
    if (p != NULL) bstep = strtok(s, " ");
  } while ((p != NULL) && ((bstep[0] == '#') || (strlen(s) == 1)));

  if (p == NULL) {
    printf("Error in file format: %s\n", filename);
    exit(0);
  }

  idx = 0;

  do {
    list[idx] = translate(bstep[0])*number_of_bases+translate(bstep[1]);
    if (list[idx] < 0) {
      printf("Error in base pair step format: file %s %c %c", filename,
         bstep[0], bstep[1]);
      exit(0);
    }
    idx++;
  } while (((bstep = strtok(NULL, " ")) != NULL) && (bstep[0] != '#'));

  int a = 0;

  do {
    char *par;
    if (((p=fgets(s, 2048, f)) != NULL) && (sscanf(s, FLOAT_ESTR, &d) == 1)) {
      int q = 0;
      par = strtok(s, " ");
      sscanf(par, FLOAT_ESTR, &d);
      do {
	if (a > 5) {
	  printf("Too many parameters in theta/rho file!\n");
	  exit(0);
	}
        I[list[q]].setv(a+1,1,d);
        q++;
        par = strtok(NULL, " ");
      } while ((q < idx) && (sscanf(par, FLOAT_ESTR, &d) == 1));
      if (q != idx) {
	printf("not enough specified parameters in %s %d\n", filename, q);
        exit(0);
      }
      a++;
    }

  } while ((a < 6) && (p != NULL));

  if (p == NULL) {
    printf("not enough specified couplings in %s", filename);
  }

  for (int i = 0; i < number_of_bases; i++) {
    for (int j = 0; j < number_of_bases; j++) {
      if (locate_in_list(list, i*number_of_bases+j, idx) == -1) {
        if (locate_in_list(list, (number_of_bases-2-j)*number_of_bases+(number_of_bases-2-i), idx) == -1) {
          printf("Incomplete list of parameters in %s %d %d", filename, i, j);
        } else {
          I[i*number_of_bases+j] = I[(number_of_bases-2-j)*number_of_bases+(number_of_bases-2-i)];
	  I[i*number_of_bases+j].setv(1,1,-I[i*number_of_bases+j](1,1));
	  I[i*number_of_bases+j].setv(4,1,-I[i*number_of_bases+j](4,1));
        }
      }
    }
  }
  fclose(f);
}






bdna::bdna() {

  nsteps = 0;
  bpstep = new int[maxseq];

//  for (int i = 0; i < nsteps; i++) v[i].setsize(6,1);  

  for (int i = 0; i < number_of_bases*number_of_bases; i++) {
    F[i].setsize(6,6);
    I[i].setsize(6,1);
  }

  readFparameters("FGH.dat", F);
  readIparameters("tp0.dat", I);

} 

bdna::bdna(int ns, int q) {
  nsteps = ns;
  v = new matrix[nsteps];
  bpstep = new int[nsteps];
  for (int i = 0; i < nsteps; i++) v[i].setsize(6,1);
}

bdna::bdna(int ns) {

  nsteps = ns;
  bpstep = new int[nsteps];

  beta = 1.0;
  
  for (int i = 0; i < number_of_bases*number_of_bases; i++) {
    F[i].setsize(6,6);
    I[i].setsize(6,1);
  }

  readFparameters("FGH.dat", F);
  readIparameters("tp0.dat", I);

  v = new matrix[nsteps];
  for (int i = 0; i < nsteps; i++) v[i].setsize(6,1);

}

bdna::~bdna() {
  delete [] bpstep;
  delete [] v;
  delete [] seq;
}


void bdna::initialize_bdna(dna D) {

  char step[2];

  ion = 0.10;
  
  nsteps = D.nbp-1;

  seq = new char[D.nbp+1];
  strcpy(seq, D.sequence);

  for (int i = 0; i < nsteps; i++) {
    step[0]=D.sequence[i];
    step[1]=D.sequence[i+1];
    bpstep[i] = translate(step[0])*number_of_bases+translate(step[1]);    
  }

  v = new matrix[nsteps];

  for (int i = 0; i < nsteps; i++) v[i] = I[bpstep[i]];

}

void bdna::add_phosphates() {
  phosphate = new matrix[2*nsteps];
  for (int i = 0; i < 2*nsteps; i++) {
    phosphate[i].setsize(4,1);
	phosphate[i].setv(4,1, -0.24);
  }
}

void bdna::delete_phosphates() {
  delete [] phosphate;
}

void bdna::set_beta(FLOAT be) {
  beta = be;
}

void bdna::resize(int ns) {
  if (v != NULL) delete [] v;
  nsteps = ns;
  v = new matrix[nsteps];
}

void bdna::addbase(char c) {
  if (nsteps < maxseq) {
    bpstep[nsteps] = (bpstep[nsteps-1]%number_of_bases)*number_of_bases
      + translate(c);
    nsteps++;
    return;
  }
  printf("Could not add any more bases!\n\n\n");
}


void bdna::set_straight() {
  for (int i = 0; i < nsteps; i++) {
    v[i].setv(1,1, 0.0);
    v[i].setv(2,1, 0.0);
    v[i].setv(3,1, 36.0);
    v[i].setv(4,1, 0.0);
    v[i].setv(5,1, 0.0);
    v[i].setv(6,1, 3.4);
  }
}


void bdna::stress_free_state(int write) {
  matrix W = identity(4);
  FLOAT tw = 0.0;
  for (int i = 0; i < nsteps; i++) {
    v[i] = I[bpstep[i]];
    W = W * calculateW(v[i]);
    tw += v[i](3,1);
  }

  if (write) {
    printf("\nStress-free A_0:N matrix:\n");
    writematrix(stdout, W);
    printf("\nEnd-to-end [theta, rho]\n");
    writematrix(stdout, calculatetp(W));
    printf("\n");
    printf("Tw (number of turns) ~= %f\n", tw/360.0);
  }

}


FLOAT bdna::calculateR() {

  FLOAT R = 0.0;

  FLOAT x, y, z;

  x = 0.0; y = 0.0; z = 0.0;

  matrix temp = identity(4);

  for (int i = 0; i < nsteps; i++) {
    temp = temp*calculateW(v[i]);
    x += temp(1,4);
    y += temp(2,4);
    z += temp(3,4);
  }

  x /= ((FLOAT)nsteps+1.0);
  y /= ((FLOAT)nsteps+1.0);
  z /= ((FLOAT)nsteps+1.0);

  R = x*x+y*y+z*z;

  temp = identity(4);
  for (int i = 0; i < nsteps; i++) {
    temp = temp * calculateW(v[i]);
    R += (temp(1,4)-x)*(temp(1,4)-x)+(temp(2,4)-y)*(temp(2,4)-y)+(temp(3,4)-z)*(temp(3,4)-z);
  }

  R /= (FLOAT)nsteps+1.0;

  return sqrt(R);

}


FLOAT bdna::calculate_elenergy(int a, int b) {

  FLOAT E = 0.0;
  
  matrix tp(6,1);
  matrix tpt(1,6);

  for (int i = a; i <= b; i++) {

    //v[i].writematrix(stdout);

    tp = v[i] - I[bpstep[i]];

    //      tp.writematrix(stdout);

    for (int j = 0; j < 6; j++) tpt.setv(1,j+1, tp(j+1,1));

    E += 0.5 * (tpt*F[bpstep[i]]*tp)(1,1);

  }

  return E;

}

FLOAT bdna::calculate_elenergy(int a) {
  
  matrix tp(6,1);
  matrix tpt(1,6);

    tp = v[a] - I[bpstep[a]];

    if (tp.m == 0) {
      int m;
      printf("%d \n", a);
      writematrix(stdout, v[a]);
      writematrix(stdout, I[bpstep[a]]);
      scanf("%d", &m);
    }
    
    for (int j = 0; j < 6; j++) tpt.setv(1,j+1, tp(j+1,1));

  return 0.5 * (tpt*F[bpstep[a]]*tp)(1,1);

}


FLOAT bdna::calculate_elenergy(matrix tp0, int m) {
  
  matrix tpt(1,6);

  matrix tp = tp0 - I[bpstep[m]];

  for (int j = 0; j < 6; j++) tpt.setv(1,j+1, tp(j+1,1));

  return 0.5 * (tpt*F[bpstep[m]]*tp)(1,1);
 
}


FLOAT bdna::calculate_elenergy() {

  FLOAT E = 0.0;
  
  matrix tp;
  matrix temp;

  for (int i = 0; i < nsteps; i++) {

    tp = v[i] - I[bpstep[i]];

    temp = F[bpstep[i]]*tp;

    for (int j = 1; j <= 6; j++) E += temp(j,1)*tp(j,1);

  }

  E *= 0.5;

  return E;

}

void bdna::print_end() {
  matrix temp = identity(4);
  matrix tp;
  for (int i = 0; i < nsteps; i++) temp = temp*calculateW(v[i]);
  tp = calculatetp(temp);
  printf("\nEnd-to-end theta/rho\n");
  writematrix(stdout, tp);
}

void bdna::print_energy() {
  printf("Total Energy = %lf kT\n", calculate_elenergy()+calculate_penergy());
  printf("Elastic Energy = %lf kT,  Contact (r-12) Energy = %lf kT\n",
	 calculate_elenergy(), calculate_penergy());
#ifdef ELECTROSTATICS
  matrix W;
  printf("Electrostatic Energy = %lf kT\n", calculate_esenergy(W));
#endif
}

FLOAT bdna::calculate_energy(matrix &W) {
  FLOAT eng = calculate_elenergy();
#ifdef ELECTROSTATICS
  eng += calculate_esenergy(W);
#endif
#ifdef SELF_CONTACT
  eng += calculate_penergy();
#endif
#ifndef ELECTROSTATICS
  W = identity(4);
  for (int i = 0 ; i < nsteps; i++)
    W = W * calculateW(v[i]);
#endif
  return eng;
}


// self-contact repulsion
FLOAT bdna::calculate_penergy() {

  FLOAT E = 0.0;

  matrix temp = identity(4);
  int nst = (int)(nsteps/3)+1;
  FLOAT *xp = new FLOAT[nst];
  FLOAT *yp = new FLOAT[nst];
  FLOAT *zp = new FLOAT[nst];

  FLOAT x,y,z,r;
  
  int m = 1;
  
  for (int i = 0; i < nsteps; i++) {
    if (i % 3 == m) { 
      int ni = i / 3;
      xp[ni] = temp(1,4);
      yp[ni] = temp(2,4);
      zp[ni] = temp(3,4);
      for (int j = 1; j < ni-4; j++) {
	x = xp[ni]-xp[j];
	y = yp[ni]-yp[j];
	z = zp[ni]-zp[j];
	r = sqrt(x*x+y*y+z*z);
	if ((r < 30.0) && ((ni-j+4) < nst)) {
	  E += 12855.0*pow((r/10.0), -12.0);
	}
      }
    }
    temp = temp*calculateW(v[i]);
  }
  
  delete [] xp;
  delete [] yp;
  delete [] zp;
  
  return E;
   
}


// electrostatic (dampened by ion screening) energy
FLOAT bdna::calculate_esenergy(matrix &W) {

  FLOAT E = 0.0;

  /*
  const FLOAT eps = 77.7;
  const FLOAT e = 0.24*1.6021773e-19;
  const FLOAT eps0 = 8.854188e-32;
  const FLOAT kT = 1.380658e-13*300.0;
  const FLOAT coeff = e*e/(4*3.14159*eps0*eps*kT);
  */

  FLOAT kap = 0.329*sqrt(ion);

  matrix W1 = identity(4);
  FLOAT r,x,y,z;

  FLOAT q1, q2;

  for (int i = 0; i < nsteps-1; i++) {
 
    FLOAT X1 = W1(1,4);
    FLOAT Y1 = W1(2,4);
    FLOAT Z1 = W1(3,4);

    phosphate[i].setv(1,1, X1);
    phosphate[i].setv(2,1, Y1);
    phosphate[i].setv(3,1, Z1);

    q1 = phosphate[i](4,1);

    W1 = W1 * calculateW(v[i]);

    for (int j = (i == (nsteps - 2) ? 1 : 0); j < i-1; j++) {
      x = X1-phosphate[j](1,1);
      y = Y1-phosphate[j](2,1);
      z = Z1-phosphate[j](3,1); 
	  q2 = phosphate[j](4,1);
      r = sqrt(x*x+y*y+z*z);
      E += q1*q2*exp(-kap*r)/r;
    }

  }

  W = W1;

  E *= coeff;

  return E;

}


FLOAT bdna::calculate_twist_closed() {

  FLOAT twist = 0.0;
  matrix temp = identity(4);
                                                                                                              
  FLOAT *si_x, *si_y, *si_z;
  FLOAT *o_x, *o_y, *o_z;
  FLOAT ti_x, ti_y, ti_z;
  FLOAT bi_x, bi_y, bi_z;
  FLOAT ti2_x, ti2_y, ti2_z;
  FLOAT ts_x, ts_y, ts_z;
  FLOAT p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  FLOAT bp1_x, bp1_y, bp1_z, bp2_x, bp2_y, bp2_z;
                                                                                                              
  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;
                                                                                                              
  si_x = new FLOAT[nsteps+2];
  si_y = new FLOAT[nsteps+2];
  si_z = new FLOAT[nsteps+2];
                                                                                                              
  o_x = new FLOAT[nsteps+2];
  o_y = new FLOAT[nsteps+2];
  o_z = new FLOAT[nsteps+2];
                                                                                                              
                                                                                                              
  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    temp = temp * calculateW(v[i]);
  }
  o_x[nsteps] = temp(1,4); o_y[nsteps] = temp(2,4); o_z[nsteps] = temp(3,4);
  si_x[nsteps] = temp(1,1); si_y[nsteps] = temp(2,1); si_z[nsteps] = temp(3,1); 
  o_x[nsteps+1] = 0.0; o_y[nsteps+1] = 0.0; o_z[nsteps+1] = 0.0;
  si_x[nsteps+1] = 1.0; si_y[nsteps+1] = 0.0; si_z[nsteps+1] = 0.0;                                                                                   
/* 
 printf("origins:\n");
  for (int i = 0; i <= nsteps; i++) { 
    printf("%lf %lf %lf\n", o_x[i], o_y[i], o_z[i]);
  }
  printf("\nshort axes:\n");
  for (int i = 0; i <= nsteps; i++) {    
    printf("%lf %lf %lf\n", si_x[i], si_y[i], si_z[i]);  
  }
*/

  for (int i = 0; i <= nsteps; i++) {
    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];
    normalize(ti_x, ti_y, ti_z);
    if (i != nsteps) {
      ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];
    } else {
       ti2_x = o_x[1] - o_x[i+1]; ti2_y = o_y[1] - o_y[i+1]; ti2_z = o_z[1] - o_z[i+1];
    }      
                                                                                                        
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                              
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
//    normalize(bp2_x, bp2_y, bp2_z);
                                                                                                              
    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;
    FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
 
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;
    if (dotp2 > 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;                                                                                                          
    FLOAT th1 = acos(dotp);
    if (ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z < 0.0) th1 = -th1;
                                                                                                              
    FLOAT th2 = acos(dotp2);
    if (ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z < 0.0) th2 = -th2;
                                                                                                              
    FLOAT diff = th2-th1;
                                                                                                              
    while (diff < -M_PI) diff += 2.0*M_PI;
    while (diff > M_PI) diff -= 2.0*M_PI;
                                                                                                              
//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);
                                                                                                              
    twist += diff;

    //printf("Step %d, twist (forward) %lf\n", i+1, diff*180.0/M_PI);

  }

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;
                                                                                                              
  temp = identity(4);
                                                                                                              
  for (int i = nsteps+1; i > 0; i--) {
                                                                                                              
    ti_x = o_x[i-1]-o_x[i]; ti_y = o_y[i-1]-o_y[i]; ti_z = o_z[i-1]-o_z[i];
                                                                                                              
    normalize(ti_x, ti_y, ti_z);
                                                                                                              
    if (i > 1) {
      ti2_x = o_x[i-2]-o_x[i-1]; ti2_y = o_y[i-2]-o_y[i-1]; ti2_z = o_z[i-2]-o_z[i-1];
    } else {
      ti2_x = o_x[nsteps]-o_x[i-1]; ti2_y = o_y[nsteps]-o_y[i-1]; ti2_z = o_z[nsteps]-o_z[i-1];
    }
                                                                                                              
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i-1], si_y[i-1], si_z[i-1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                              
                                                                                                              
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
 //   normalize(bp1_x, bp1_y, bp1_z);
                                                                                                              
                                                                                                              
    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
//    normalize(bp2_x, bp2_y, bp2_z);
                                                                                                             
                          
    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;                                     FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;
    if (dotp2 > 1.0) dotp2 = 1.0;                                                   if (dotp2 < -1.0) dotp2 = -1.0;                                                                                                                                 
    FLOAT th1 = acos(dotp);
    if (ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z < 0.0) th1 = -th1;
                                                                                                                                               
    FLOAT th2 = acos(dotp2);
    if (ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z < 0.0) th2 = -th2;
                                                                                                              
    FLOAT diff = th2 - th1;
                                                                                                              
    while (diff < -M_PI) diff += 2*M_PI;
    while (diff > M_PI) diff -= 2*M_PI;
                                                                                                              
//    printf("%lf %lf %lf %lf %lf\n", diff, acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z), ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z, acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z), ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z);
                                                                                                              
    twist += diff;

//    printf("Step %d, twist (reverse) %lf\n", i, diff*180.0/M_PI);                                                                                     
  }
                                                                                                              
  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;
                                                                                                              
//  twist = twist/360.0;
  return twist/(4.0*M_PI);

}


FLOAT bdna::calculate_twist_open() {

  FLOAT twist = 0.0;
  matrix temp = identity(4);
                                                                                                              
  FLOAT *si_x, *si_y, *si_z;
  FLOAT *o_x, *o_y, *o_z;
  FLOAT ti_x, ti_y, ti_z;
  FLOAT bi_x, bi_y, bi_z;
  FLOAT ti2_x, ti2_y, ti2_z;
  FLOAT ts_x, ts_y, ts_z;
  FLOAT p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  FLOAT bp1_x, bp1_y, bp1_z, bp2_x, bp2_y, bp2_z;

  matrix tptemp(6,1);
  tptemp.setv(1,1, 0.0);
  tptemp.setv(2,1, 0.0);
  tptemp.setv(3,1, 34.28);
  tptemp.setv(4,1, 0.0);
  tptemp.setv(5,1, 0.0);
  tptemp.setv(6,1, 3.4);
                       
  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;
                                                                                                              
  si_x = new FLOAT[nsteps+1];
  si_y = new FLOAT[nsteps+1];
  si_z = new FLOAT[nsteps+1];
                                                                                                              
  o_x = new FLOAT[nsteps+1];
  o_y = new FLOAT[nsteps+1];
  o_z = new FLOAT[nsteps+1];
                                                                                                              
                                                                                                              
  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    temp = temp * calculateW(v[i]);
  }
  o_x[nsteps] = 0.0; o_y[nsteps] = 0.0; o_z[nsteps] = 0.0;
  si_x[nsteps] = 1.0; si_y[nsteps] = 0.0; si_z[nsteps] = 0.0;
 
                                                                                                              
  for (int i = 0; i < nsteps; i++) {
    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];
    normalize(ti_x, ti_y, ti_z);
    if (i != nsteps-1) {
      ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];
    } else {
	  ti2_x = o_x[1] - o_x[0]; ti2_y = o_y[1] - o_y[0]; ti2_z = o_z[1] - o_z[0];
    }
                                                                                                              
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                              
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
//    normalize(bp2_x, bp2_y, bp2_z);
                                                                                                              
                                                                                                              
    FLOAT th1 = acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z);
    if (ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z < 0.0) th1 = -th1;
                                                                                                              
    FLOAT th2 = acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z);
    if (ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z < 0.0) th2 = -th2;
                                                                                                              
    FLOAT diff = th2-th1;
                                                                                                              
    while (diff < -M_PI) diff += 2.0*M_PI;
    while (diff > M_PI) diff -= 2.0*M_PI;
                                                                                                              
//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);
                                                                                                              
    twist += diff;

  //  printf("Step %d, twist (forward) %lf\n", i+1, diff*180.0/M_PI);

  }

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;

//  printf("%lf\n", twist/(2.0*M_PI));
//  return twist/(2.0*M_PI);

  FLOAT twist2;

  temp = identity(4);
  matrix vb(6,1);
  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    vb = v[nsteps-1-i];
    vb.setv(1,1,-vb(1,1));
    vb.setv(4,1,-vb(4,1));
    temp = temp * calculateW(vb);
  }
  o_x[nsteps] = 0.0; o_y[nsteps] = 0.0; o_z[nsteps] = 0.0;
  si_x[nsteps] = 1.0; si_y[nsteps] = 0.0; si_z[nsteps] = 0.0;
  
                                                                                                              
  for (int i = 0; i < nsteps; i++) {    
    
    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];

    normalize(ti_x, ti_y, ti_z);

    if (i != nsteps-1) {
      ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];
    } else {
      ti2_x = o_x[1] - o_x[0]; ti2_y = o_y[1] - o_y[0]; ti2_z = o_z[1] - o_z[0];
    }
                     
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                              
                                                                                                              
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
 //   normalize(bp1_x, bp1_y, bp1_z);
                                                                                                              
                                                                                                              
    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
//    normalize(bp2_x, bp2_y, bp2_z);
                                                                                                             
                                                                                                              
    FLOAT th1 = acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z);
    if (ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z < 0.0) th1 = -th1;
                                                                                                              
                                                                                                              
    FLOAT th2 = acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z);
    if (ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z < 0.0) th2 = -th2;
                                                                                                              
    FLOAT diff = th2 - th1;
                                                                                                              
    while (diff < -M_PI) diff += 2.0*M_PI;
    while (diff > M_PI) diff -= 2.0*M_PI;
                                                                                                              
//    printf("%lf %lf %lf %lf %lf\n", diff, acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z), ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z, acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z), ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z);
                                                                                                              
    twist += diff;

  //  printf("Step %d, twist (reverse) %lf\n", i, diff*180.0/M_PI);                                                                                                              
    twist2 += diff;                                                                                                          
  }
                                                                                                              
  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;
                
//  printf("%lf\n", twist/(4.0*M_PI));                                                                                              
//  twist = twist/360.0;
  return twist/(4.0*M_PI);
}



FLOAT bdna::calculate_twist_print() {
  FLOAT twist = 0.0;
  matrix temp = identity(4);

  FLOAT *si_x, *si_y, *si_z;
  FLOAT *o_x, *o_y, *o_z;
  FLOAT ti_x, ti_y, ti_z;
  FLOAT bi_x, bi_y, bi_z;
  FLOAT ti2_x, ti2_y, ti2_z;
  FLOAT ts_x, ts_y, ts_z;
  FLOAT p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  FLOAT bp1_x, bp1_y, bp1_z, bp2_x, bp2_y, bp2_z;

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;

  FLOAT ts21_x, ts21_y, ts21_z;
  FLOAT ts12_x, ts12_y, ts12_z;
  FLOAT ps21_x, ps21_y, ps21_z;
  FLOAT ps12_x, ps12_y, ps12_z;

  si_x = new FLOAT[nsteps+2];
  si_y = new FLOAT[nsteps+2];
  si_z = new FLOAT[nsteps+2];

  o_x = new FLOAT[nsteps+2];
  o_y = new FLOAT[nsteps+2];
  o_z = new FLOAT[nsteps+2];

  FLOAT sum1 = 0.0; 
  FLOAT sum2 = 0.0;

  FLOAT dir = 0.0;
  FLOAT undir = 0.0;

  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    temp = temp * calculateW(v[i]);
  }
  si_x[nsteps] = temp(1,1); si_y[nsteps] = temp(2,1); si_z[nsteps] = temp(3,1);
  o_x[nsteps] = temp(1,4); o_y[nsteps] = temp(2,4); o_z[nsteps] = temp(3,4);

  matrix temp2(6,1); temp2.setv(3,1, 34.28); temp2.setv(6,1,3.4);
         temp2.setv(1,1, 0.0); temp2.setv(2,1, 0.0); temp2.setv(4,1, 0.0); temp2.setv(5,1, 0.0);
  temp = temp * calculateW(temp2);
  si_x[nsteps+1] = temp(1,1); si_y[nsteps+1] = temp(2,1); si_z[nsteps+1] = temp(3,1);
  o_x[nsteps+1] = temp(1,4); o_y[nsteps+1] = temp(2,4); o_z[nsteps+1] = temp(3,4);


  for (int i = 0; i < nsteps; i++) {
    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];
    normalize(ti_x, ti_y, ti_z);
    ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];

    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);

    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    normalize(bp1_x, bp1_y, bp1_z);

    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
    normalize(bp2_x, bp2_y, bp2_z);

                                                                                
    crossr(ts12_x, ts12_y, ts12_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(ps12_x, ps12_y, ps12_z, ti_x, ti_y, ti_z, ts12_x, ts12_y, ts12_z);
    normalize(ps12_x, ps12_y, ps12_z);
                                                                                
    crossr(ts21_x, ts21_y, ts21_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);    crossr(ps21_x, ps21_y, ps21_z, ti_x, ti_y, ti_z, ts21_x, ts21_y, ts21_z);
    normalize(ps21_x, ps21_y, ps21_z);
                                                                                


    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;

    FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
    if (dotp2 >= 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;

    FLOAT dotp3 = ps12_x*ps21_x + ps12_y*ps21_y + ps12_z*ps21_z;
    if (dotp3 >= 1.0) dotp3 = 1.0;
    if (dotp3 < -1.0) dotp3 = -1.0;

    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, ps21_x, ps21_y, ps21_z, ps12_x, ps12_y, ps12_z);
                                                                                
    FLOAT th = acos(dotp3);
    if ((lx*ti_x+ly*ti_y+lz*ti_z) < 0.0) th = -th;
                                                 

    FLOAT th1 = acos(dotp);
    if ((ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z) < 0.0) th1 = -th1;

    FLOAT th2 = acos(dotp2);
    if ((ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z) < 0.0) th2 = -th2;

   // if (th2 < th1) th2 += 2.0*M_PI;

    FLOAT dth = th2-th1;
	if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;

    printf("Step %d, dir = %lf deg, undir = %lf deg\n", i+1, dth*180.0/M_PI, th*180.0/M_PI); 


//    printf("th1 = %lf th2 = %lf,  %lf %lf\n", th1, th2, dotp, dotp2);

//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);

    twist += diff;

  }

  printf("\n\nTwist one way = %lf\n\n", twist/(2.0*M_PI));

  return twist/(2.0*M_PI);

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;

  dir = 0.0;
  undir = 0.0;
                                                                                
  temp = identity(4);

  FLOAT twist2 = 0.0;


  for (int i = nsteps; i > 0; i--) {
                                                                                                              
    ti_x = o_x[i-1]-o_x[i]; ti_y = o_y[i-1]-o_y[i]; ti_z = o_z[i-1]-o_z[i];
                                                                                                              
    normalize(ti_x, ti_y, ti_z);
                                                                                                              
    if (i > 1) {
      ti2_x = o_x[i-2]-o_x[i-1]; ti2_y = o_y[i-2]-o_y[i-1]; ti2_z = o_z[i-2]-o_z[i-1];
    } else {
      ti2_x = o_x[nsteps-1]-o_x[i-1]; ti2_y = o_y[nsteps-1]-o_y[i-1]; ti2_z = o_z[nsteps-1]-o_z[i-1];
    }
                                                                                                              
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i-1], si_y[i-1], si_z[i-1], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i-1], si_y[i-1], si_z[i-1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                              
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    normalize(bp1_x, bp1_y, bp1_z);

    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
    normalize(bp2_x, bp2_y, bp2_z);

                                                                                
    crossr(ts12_x, ts12_y, ts12_z, si_x[i-1], si_y[i-1], si_z[i-1], ti_x, ti_y, ti_z);
    crossr(ps12_x, ps12_y, ps12_z, ti_x, ti_y, ti_z, ts12_x, ts12_y, ts12_z);
    normalize(ps12_x, ps12_y, ps12_z);
                                                                                
    crossr(ts21_x, ts21_y, ts21_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);    crossr(ps21_x, ps21_y, ps21_z, ti_x, ti_y, ti_z, ts21_x, ts21_y, ts21_z);
    normalize(ps21_x, ps21_y, ps21_z);
                                                                                


    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;

    FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
    if (dotp2 > 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;

    FLOAT dotp3 = ps12_x*ps21_x + ps12_y*ps21_y + ps12_z*ps21_z;
    if (dotp3 > 1.0) dotp3 = 1.0;
    if (dotp3 < -1.0) dotp3 = -1.0;

    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, ps21_x, ps21_y, ps21_z, ps12_x, ps12_y, ps12_z);
                                                                                
    FLOAT th = acos(dotp3);
    if ((lx*ti_x+ly*ti_y+lz*ti_z) < 0.0) th = -th;
                                                                      
																	
    FLOAT th1 = acos(dotp);
    if ((ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z) < 0.0) th1 = -th1;

    FLOAT th2 = acos(dotp2);
    if ((ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z) < 0.0) th2 = -th2;

   // if (th2 < th1) th2 += 2.0*M_PI;
    FLOAT dth = th2-th1;
	if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;

    undir += th;
    dir += dth;


//    printf("th1 = %lf th2 = %lf,  %lf %lf\n", th1, th2, dotp, dotp2);

//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);

    printf("Step %d, dir = %lf, undir = %lf\n", i, dth*180.0/M_PI, th*180.0/M_PI);

    twist += diff;
    twist2 += diff;
  }

  printf("Twist second way = %lf\n", twist2/(2.0*M_PI));

  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;

  return twist/(4.0*M_PI);

}



FLOAT bdna::calculate_twist() {
  FLOAT twist = 0.0;
  FLOAT ctwist = 0.0;
  matrix temp = identity(4);


  FLOAT *si_x, *si_y, *si_z;
  FLOAT *o_x, *o_y, *o_z;
  FLOAT ti_x, ti_y, ti_z;
  FLOAT bi_x, bi_y, bi_z;
  FLOAT ti2_x, ti2_y, ti2_z;
  FLOAT ts_x, ts_y, ts_z;
  FLOAT p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  FLOAT bp1_x, bp1_y, bp1_z, bp2_x, bp2_y, bp2_z;
  FLOAT ts21_x, ts21_y, ts21_z;
  FLOAT ts12_x, ts12_y, ts12_z;
  FLOAT ps21_x, ps21_y, ps21_z;
  FLOAT ps12_x, ps12_y, ps12_z;

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;

  matrix tptend;

  si_x = new FLOAT[nsteps+1];
  si_y = new FLOAT[nsteps+1];
  si_z = new FLOAT[nsteps+1];

  o_x = new FLOAT[nsteps+1];
  o_y = new FLOAT[nsteps+1];
  o_z = new FLOAT[nsteps+1];
  FLOAT dir = 0.0;
  FLOAT undir = 0.0;
  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    temp = temp * calculateW(v[i]);
  }
  o_x[nsteps] = 0.0; o_y[nsteps] = 0.0; o_z[nsteps] = 0.0;
  si_x[nsteps] = 1.0; si_y[nsteps] = 0.0; si_z[nsteps] = 0.0;
  temp = invert(temp*(invert(calculateW(v[nsteps-1]))));
  tptend = calculatetp(temp);


  //printf("%lf degree closing\n", tptend(3,1));

  //printf("forward twists: ");

  for (int i = 0; i < nsteps-3; i++) {

    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];
    normalize(ti_x, ti_y, ti_z);
    if (i != nsteps-1) {
      ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];
    } else {
      ti2_x = o_x[1] - o_x[i+1]; ti2_y = o_y[1] - o_y[i+1]; ti2_z = o_z[1] - o_z[i+1];
    }

    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);

    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    normalize(bp1_x, bp1_y, bp1_z);

    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
    normalize(bp2_x, bp2_y, bp2_z);

                                                                                
    crossr(ts12_x, ts12_y, ts12_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(ps12_x, ps12_y, ps12_z, ti_x, ti_y, ti_z, ts12_x, ts12_y, ts12_z);
    normalize(ps12_x, ps12_y, ps12_z);
                                                                                
    crossr(ts21_x, ts21_y, ts21_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(ps21_x, ps21_y, ps21_z, ti_x, ti_y, ti_z, ts21_x, ts21_y, ts21_z);
    normalize(ps21_x, ps21_y, ps21_z);
                                                                                


    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;

    FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
    if (dotp2 >= 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;

    FLOAT dotp3 = ps12_x*ps21_x + ps12_y*ps21_y + ps12_z*ps21_z;
    if (dotp3 >= 1.0) dotp3 = 1.0;
    if (dotp3 < -1.0) dotp3 = -1.0;

    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, ps21_x, ps21_y, ps21_z, ps12_x, ps12_y, ps12_z);
                                                                                

    FLOAT th = acos(dotp3);
    if ((lx*ti_x+ly*ti_y+lz*ti_z) < 0.0) th = -th;
                                                 

    FLOAT th1 = acos(dotp);
    if ((ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z) < 0.0) th1 = -th1;

    FLOAT th2 = acos(dotp2);
    if ((ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z) < 0.0) th2 = -th2;

   // if (th2 < th1) th2 += 2.0*M_PI;

    FLOAT dth = th2-th1;
	if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;

//    printf("Step %d, dir = %lf, undir = %lf, total = %lf\n", i+1, dth*180.0/M_PI, th*180.0/M_PI, diff*180.0/M_PI); 


//    printf("th1 = %lf th2 = %lf,  %lf %lf\n", th1, th2, dotp, dotp2);

//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);

    ctwist += diff;
	
  }
  ctwist += (v[nsteps-3](3,1)+v[nsteps-2](3,1)+tptend(3,1))*(M_PI/180.0);


  for (int i = 0; i < nsteps; i++) {

    ti_x = o_x[i+1]-o_x[i]; ti_y = o_y[i+1]-o_y[i]; ti_z = o_z[i+1]-o_z[i];
    normalize(ti_x, ti_y, ti_z);
    if (i != nsteps-1) {
      ti2_x = o_x[i+2] - o_x[i+1]; ti2_y = o_y[i+2] - o_y[i+1]; ti2_z = o_z[i+2] - o_z[i+1];
    } else {
      ti2_x = o_x[1] - o_x[i+1]; ti2_y = o_y[1] - o_y[i+1]; ti2_z = o_z[1] - o_z[i+1];
    }

    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i+1], si_y[i+1], si_z[i+1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);

    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
    normalize(bp1_x, bp1_y, bp1_z);

    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
    normalize(bp2_x, bp2_y, bp2_z);

                                                                                
    crossr(ts12_x, ts12_y, ts12_z, si_x[i+1], si_y[i+1], si_z[i+1], ti_x, ti_y, ti_z);
    crossr(ps12_x, ps12_y, ps12_z, ti_x, ti_y, ti_z, ts12_x, ts12_y, ts12_z);
    normalize(ps12_x, ps12_y, ps12_z);
                                                                                
    crossr(ts21_x, ts21_y, ts21_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);
    crossr(ps21_x, ps21_y, ps21_z, ti_x, ti_y, ti_z, ts21_x, ts21_y, ts21_z);
    normalize(ps21_x, ps21_y, ps21_z);
                                                                                


    FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;

    FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;
    if (dotp2 >= 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;

    FLOAT dotp3 = ps12_x*ps21_x + ps12_y*ps21_y + ps12_z*ps21_z;
    if (dotp3 >= 1.0) dotp3 = 1.0;
    if (dotp3 < -1.0) dotp3 = -1.0;

    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, ps21_x, ps21_y, ps21_z, ps12_x, ps12_y, ps12_z);
                                                                                

    FLOAT th = acos(dotp3);
    if ((lx*ti_x+ly*ti_y+lz*ti_z) < 0.0) th = -th;
                                                 

    FLOAT th1 = acos(dotp);
    if ((ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z) < 0.0) th1 = -th1;

    FLOAT th2 = acos(dotp2);
    if ((ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z) < 0.0) th2 = -th2;

   // if (th2 < th1) th2 += 2.0*M_PI;

    FLOAT dth = th2-th1;
	if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;

//    printf("Step %d, dir = %lf, undir = %lf, total = %lf\n", i+1, dth*180.0/M_PI, th*180.0/M_PI, diff*180.0/M_PI); 


//    printf("th1 = %lf th2 = %lf,  %lf %lf\n", th1, th2, dotp, dotp2);

//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);

    twist += diff;


  }

  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;

  if (fabs(twist-ctwist) > M_PI) {
//	printf("Bad twists %lf approx, %lf actual, will correct to %lf\n", ctwist/(2.0*M_PI), twist/(2.0*M_PI), (twist > ctwist ? (twist/(2.0*M_PI))-1.0 : (twist/(2.0*M_PI))+1.0));
	return (twist > ctwist ? (twist/(2.0*M_PI))-1.0 : (twist/(2.0*M_PI))+1.0);
  }

  return twist/(2.0*M_PI);
  

//  printf("\ntwist one way = %lf, directional = %lf, non-directional = %lf\n", twist/(2.0*M_PI), dir/(2.0*M_PI), undir/(2.0*M_PI));

  bi_x = 0.0;
  bi_y = 1.0;
  bi_z = 0.0;
                                                                                
  temp = identity(4);

  FLOAT twist2 = 0.0;
  dir = 0.0;
  undir = 0.0;
  for (int i = nsteps; i > 0; i--) {

    ti_x = o_x[i-1]-o_x[i]; ti_y = o_y[i-1]-o_y[i]; ti_z = o_z[i-1]-o_z[i];

    normalize(ti_x, ti_y, ti_z);

    if (i > 1) {
      ti2_x = o_x[i-2]-o_x[i-1]; ti2_y = o_y[i-2]-o_y[i-1]; ti2_z = o_z[i-2]-o_z[i-1];
    } else {
      ti2_x = o_x[nsteps-1]-o_x[i-1]; ti2_y = o_y[nsteps-1]-o_y[i-1]; ti2_z = o_z[nsteps-1]-o_z[i-1];
    }
 
    normalize(ti2_x, ti2_y, ti2_z);
    if (!vequal(ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z)) crossr(bi_x, bi_y, bi_z, ti_x, ti_y, ti_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi_x, bi_y, bi_z);
    crossr(ts_x, ts_y, ts_z, si_x[i-1], si_y[i-1], si_z[i-1], ti_x, ti_y, ti_z);
    crossr(p1_x, p1_y, p1_z, ti_x, ti_y, ti_z, ts_x, ts_y, ts_z);
    normalize(p1_x, p1_y, p1_z);
    crossr(ts_x, ts_y, ts_z, si_x[i-1], si_y[i-1], si_z[i-1], ti2_x, ti2_y, ti2_z);
    crossr(p2_x, p2_y, p2_z, ti2_x, ti2_y, ti2_z, ts_x, ts_y, ts_z);
    normalize(p2_x, p2_y, p2_z);
                                                                                                                  
    crossr(bp1_x, bp1_y, bp1_z, bi_x, bi_y, bi_z, p1_x, p1_y, p1_z);
  //  normalize(bp1_x, bp1_y, bp1_z);
                                                                                                                  
    crossr(bp2_x, bp2_y, bp2_z, bi_x, bi_y, bi_z, p2_x, p2_y, p2_z);
   // normalize(bp2_x, bp2_y, bp2_z);
    crossr(ts12_x, ts12_y, ts12_z, si_x[i-1], si_y[i-1], si_z[i-1], ti_x, ti_y, ti_z);
    crossr(ps12_x, ps12_y, ps12_z, ti_x, ti_y, ti_z, ts12_x, ts12_y, ts12_z);
    normalize(ps12_x, ps12_y, ps12_z);
                                                                                
    crossr(ts21_x, ts21_y, ts21_z, si_x[i], si_y[i], si_z[i], ti_x, ti_y, ti_z);    crossr(ps21_x, ps21_y, ps21_z, ti_x, ti_y, ti_z, ts21_x, ts21_y, ts21_z);
    normalize(ps21_x, ps21_y, ps21_z);
                                                                                
                      

   FLOAT dotp = bi_x*p1_x+bi_y*p1_y+bi_z*p1_z;    
   if (dotp > 1.0) dotp = 1.0;    
   if (dotp < -1.0) dotp = -1.0;

   FLOAT dotp2 = bi_x*p2_x+bi_y*p2_y+bi_z*p2_z;    
   if (dotp2 > 1.0) dotp2 = 1.0;                                                   if (dotp2 < -1.0) dotp2 = -1.0;                                                                                                                  

    FLOAT dotp3 = ps12_x*ps21_x + ps12_y*ps21_y + ps12_z*ps21_z;
    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, ps21_x, ps21_y, ps21_z, ps12_x, ps12_y, ps12_z);                                                                                
    FLOAT th = acos(dotp3);
    if ((lx*ti_x+ly*ti_y+lz*ti_z) < 0.0) th = -th;
                                                                                
    FLOAT th1 = acos(dotp);
    if ((ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z) < 0.0) th1 = -th1;
  
//    FLOAT th1 = atan2(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z, bi_x*p1_x+bi_y*p1_y+bi_z*p1_z);                                                                                                                
    FLOAT th2 = acos(dotp2);
    if ((ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z) < 0.0) th2 = -th2;

//    FLOAT th2 = atan2(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z, bi_x*p2_x+bi_y*p2_y+bi_z*p2_z);

    //if (th2 < th1) th2 += 2.0*M_PI;

    FLOAT dth = th2-th1;
	if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;


//    if (diff < -M_PI) diff += 2.0*M_PI;
//    if (diff > M_PI) diff -= 2.0*M_PI;

 //   printf(" %lf (%lf + %lf), ", diff*180.0/M_PI, th*180.0/M_PI, (th2-th1)*180.0/M_PI);

//    printf("%lf th1=%lf sin=%lf th2=%lf sin=%lf\n", diff, th1, asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z), th2, asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z));

    dir += th2-th1;
    undir += th;
    twist += diff;
    twist2 += diff;
//    if ((diff < 0.0) || (diff > M_PI)) printf("%lf\n", 180.0*diff/M_PI);

  }

// printf("\ntwist second way = %lf, directional = %lf, non-directional = %lf\n", twist2/(2.0*M_PI), dir/(2.0*M_PI), undir/(2.0*M_PI));

//  printf("calculate_twist_open says: ");
//  calculate_twist_open();

 /* FLOAT EH = 0.0;
  for (int i = 0; i < nsteps; i++) {
    EH += v[i](3,1);
  }

  printf("El Hassan says: %lf\n", EH/(360.0));
  */
//  printf("Twist = %lf\n", twist/(4.0*M_PI));

  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;

//  twist = twist/360.0;
  return twist/(4.0*M_PI);

}

/*
FLOAT bdna::calculate_twist2() {
  FLOAT twist = 0.0;
  matrix temp = identity(4);


  FLOAT *si_x, *si_y, *si_z;
  FLOAT *o_x, *o_y, *o_z;
  FLOAT ti1_x, ti1_y, ti1_z, ti3_x, ti3_y, ti3_z;
  FLOAT bi1_x, bi1_y, bi1_z, bi2_x, bi2_y, bi2_z;
  FLOAT ti2_x, ti2_y, ti2_z;
  FLOAT ts1_x, ts1_y, ts1_z, ts2_x, ts2_y, ts2_z;
  FLOAT p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  FLOAT bp1_x, bp1_y, bp1_z, bp2_x, bp2_y, bp2_z;
  FLOAT ps1_x, ps1_y, ps1_z;
  FLOAT ps2_x, ps2_y, ps2_z;

  bi1_x = 0.0;
  bi1_y = 1.0;
  bi1_z = 0.0;

  bi2_x = 0.0;
  bi2_y = 1.0;
  bi2_z = 0.0;

  si_x = new FLOAT[nsteps+1];
  si_y = new FLOAT[nsteps+1];
  si_z = new FLOAT[nsteps+1];

  o_x = new FLOAT[nsteps+1];
  o_y = new FLOAT[nsteps+1];
  o_z = new FLOAT[nsteps+1];
  FLOAT dir = 0.0;
  FLOAT undir = 0.0;
  for (int i = 0; i < nsteps; i++) {
    si_x[i] = temp(1,1); si_y[i] = temp(2,1); si_z[i] = temp(3,1);
    o_x[i] = temp(1,4); o_y[i] = temp(2,4); o_z[i] = temp(3,4);
    temp = temp * calculateW(v[i]);
  }
  o_x[nsteps] = 0.0; o_y[nsteps] = 0.0; o_z[nsteps] = 0.0;
  si_x[nsteps] = 1.0; si_y[nsteps] = 0.0; si_z[nsteps] = 0.0;

  //printf("forward twists: ");

  for (int i = 0; i < nsteps; i++) {

    if (i == 0) {
	   ti1_x = o_x[i]-o_x[nsteps-1];  ti1_y = o_y[i]-o_y[nsteps-1];  ti1_z = o_z[i]-o_z[nsteps-1];
    } else {
      ti1_x = o_x[i]-o_x[i-1]; ti1_y = o_y[i]-o_y[i-1]; ti1_z = o_z[i]-o_z[i-1];
    }
    normalize(ti1_x, ti1_y, ti1_z);


    ti2_x = o_x[i+1]-o_x[i]; ti2_y = o_y[i+1]-o_y[i]; ti2_z = o_z[i+1]-o_z[i];
    normalize(ti2_x, ti2_y, ti2_z);


    if (i != nsteps-1) {
      ti3_x = o_x[i+2] - o_x[i+1]; ti3_y = o_y[i+2] - o_y[i+1]; ti3_z = o_z[i+2] - o_z[i+1];
    } else {
      ti3_x = o_x[1] - o_x[i+1]; ti3_y = o_y[1] - o_y[i+1]; ti3_z = o_z[1] - o_z[i+1];
    }
    normalize(ti3_x, ti3_y, ti3_z);

    if (!vequal(ti1_x, ti1_y, ti1_z, ti2_x, ti2_y, ti2_z)) crossr(bi1_x, bi1_y, bi1_z, ti1_x, ti1_y, ti1_z, ti2_x, ti2_y, ti2_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi1_x, bi1_y, bi1_z);

    if (!vequal(ti3_x, ti3_y, ti3_z, ti2_x, ti2_y, ti2_z)) crossr(bi2_x, bi2_y, bi2_z, ti2_x, ti2_y, ti2_z, ti3_x, ti3_y, ti3_z);
    else {
      printf("bad binormal\n");
    }
    normalize(bi2_x, bi2_y, bi2_z);


    ts1_x = ti1_x+ti2_x;  ts1_y = ti1_y+ti2_y;  ts1_z = ti1_z+ti2_z;
    ts2_x = ti3_x+ti2_x;  ts2_y = ti3_y+ti2_y;  ts2_z = ti3_z+ti2_z;
    normalize(ts1_x, ts1_y, ts1_z);
    normalize(ts2_x, ts2_y, ts2_z);


    crossr(bp2_x, bp2_y, bp2_z, si_x[i+1], si_y[i+1], si_z[i+1], ts2_x, ts2_y, ts2_z);
    crossr(ps2_x, ps2_y, ps2_z, ts2_x, ts2_y, ts2_z, bp2_x, bp2_y, bp2_z);
    normalize(ps2_x, ps2_y, ps2_z);

    crossr(bp1_x, bp1_y, bp1_z, si_x[i], si_y[i], si_z[i], ts1_x, ts1_y, ts1_z);
    crossr(ps1_x, ps1_y, ps1_z, ts1_x, ts1_y, ts1_z, bp1_x, bp1_y, bp1_z);
    normalize(ps1_x, ps1_y, ps1_z);

                                                                               

    FLOAT dotp = bi1_x*ps1_x+bi1_y*ps1_y+bi1_z*ps1_z;
    if (dotp > 1.0) dotp = 1.0;
    if (dotp < -1.0) dotp = -1.0;

    FLOAT dotp2 = bi2_x*ps2_x+bi2_y*ps2_y+bi2_z*ps2_z;
    if (dotp2 >= 1.0) dotp2 = 1.0;
    if (dotp2 < -1.0) dotp2 = -1.0;

    FLOAT dotp3 = bi1_x*bi2_x+bi1_y*bi2_y+bi1_z*bi2_z;
    if (dotp3 >= 1.0) dotp3 = 1.0;
    if (dotp3 < -1.0) dotp3 = -1.0;


    FLOAT lx, ly, lz;
    crossr(lx, ly, lz, bi1_x, bi1_y, bi1_z, bi2_x, bi2_y, bi2_z);
                                                                               
    FLOAT th = acos(dotp3);
    if ((lx*ti2_x+ly*ti2_y+lz*ti2_z) < 0.0) {
	  if (dotp3 < 0.0) th = 2.0*M_PI-acos(dotp3);
	  else th = -th;
    }
                                                 

    FLOAT th1 = acos(dotp);
    crossr(lx, ly, lz, bi1_x, bi1_y, bi1_z, ps1_x, ps1_y, ps1_z);

    if ((ts1_x*lx+ts1_y*ly+ts1_z*lz) < 0.0) {
      if (dotp < 0.0) th1 = 2.0*M_PI-acos(dotp);
      else th1 = -th1;
    }

    FLOAT th2 = acos(dotp2);
    crossr(lx, ly, lz, bi2_x, bi2_y, bi2_z, ps2_x, ps2_y, ps2_z);

    if ((ts2_x*lx+ts2_y*ly+ts2_z*lz) < 0.0) {
	  if (dotp2 < 0.0) th2 = 2.0*M_PI-acos(dotp2);
      else th2 = -th2;
    }

   // if (th2 < th1) th2 += 2.0*M_PI;

    FLOAT dth = th2-th1;
    if (dth < -M_PI) dth += 2.0*M_PI;
	if (dth > M_PI) dth -= 2.0*M_PI;
    FLOAT diff = th + dth;

    if (diff > 1.5*M_PI) diff -= 2.0*M_PI;

 //   printf("Step %d, %lf %lf\n", i+1, diff*180.0/M_PI, v[i](3,1)); 


//    printf("th1 = %lf th2 = %lf,  %lf %lf\n", th1, th2, dotp, dotp2);

//    printf("%lf %lf %lf %lf %lf\n", diff, 180.0*acos(bi_x*p1_x+bi_y*p1_y+bi_z*p1_z)/M_PI, 180.0*asin(ti_x*bp1_x+ti_y*bp1_y+ti_z*bp1_z)/M_PI, 180.0*acos(bi_x*p2_x+bi_y*p2_y+bi_z*p2_z)/M_PI, 180.0*asin(ti2_x*bp2_x+ti2_y*bp2_y+ti2_z*bp2_z)/M_PI);

    twist += diff;


  }

  delete [] o_x;
  delete [] o_y;
  delete [] o_z;
  delete [] si_x;
  delete [] si_y;
  delete [] si_z;


//  printf("%lf\n", twist/(2.0*M_PI));

  return twist/(2.0*M_PI);

}
*/

int bdna::overlap() {
  FLOAT *x, *y, *z;
  x = new FLOAT[nsteps-15];
  y = new FLOAT[nsteps-15];
  z = new FLOAT[nsteps-15];
  matrix temp = identity(4);
  for (int i = 0; i < nsteps; i++) {
	temp = temp * calculateW(v[i]);
	if ((i >= 10) && (i < nsteps-10)) {
      x[i-10] = temp(1,4);
	  y[i-10] = temp(2,4);
	  z[i-10] = temp(3,4);
	  for (int j = i-30; j >= 0; j--)
		if (sqrt((x[i-10]-x[j])*(x[i-10]-x[j])+(y[i-10]-y[j])*(y[i-10]-y[j])+(z[i-10]-z[j])*(z[i-10]-z[j])) < DNA_DNA_DIST) {
		  delete [] x;
		  delete [] y;
		  delete [] z;
		  return 1;
		}
   }
  }

  delete [] x;
  delete [] y;
  delete [] z;
  return 0;
}


int bdna::calculate_twist0() {
  FLOAT twist = 0;
  for (int i = 0; i < nsteps; i++) {
    twist += I[bpstep[i]](3,1);
  }
  twist = twist/360.0;
  return (int)rint(twist);
}


FLOAT bdna::calculate_link3() {
  FLOAT *xx = new FLOAT[nsteps+1];
  FLOAT *xy = new FLOAT[nsteps+1];
  FLOAT *xz = new FLOAT[nsteps+1];

  FLOAT *sx = new FLOAT[nsteps+1];
  FLOAT *sy = new FLOAT[nsteps+1];
  FLOAT *sz = new FLOAT[nsteps+1];

  FLOAT *yx = new FLOAT[3*nsteps+1];
  FLOAT *yy = new FLOAT[3*nsteps+1];
  FLOAT *yz = new FLOAT[3*nsteps+1];

  FLOAT *dxx = new FLOAT[nsteps];
  FLOAT *dxy = new FLOAT[nsteps];
  FLOAT *dxz = new FLOAT[nsteps];
  FLOAT *dyx = new FLOAT[3*nsteps];
  FLOAT *dyy = new FLOAT[3*nsteps];
  FLOAT *dyz = new FLOAT[3*nsteps];

  FLOAT stx, sty, stz;
  FLOAT px, py, pz;
  FLOAT bx, by, bz;

  xx[0] = 0.0; xy[0] = 0.0; xz[0] = 0.0;
  sx[0] = 1.0; sy[0] = 0.0; sz[0] = 0.0;

  matrix temp = identity(4);

  for (int i = 0; i < nsteps; i++) {

// fill out curve one and the short axes
    dxx[i] = -temp(1,4);
    dxy[i] = -temp(2,4);
    dxz[i] = -temp(3,4);

    temp = temp*calculateW(v[i]);
    xx[i+1] = temp(1,4); xy[i+1] = temp(2,4); xz[i+1] = temp(3,4);
    sx[i+1] = temp(1,1); sy[i+1] = temp(2,1); sz[i+1] = temp(3,1);

    if (i+1 == nsteps) {
      xx[i+1] = 0.0; xy[i+1] = 0.0; xz[i+1] = 0.0;
	  sx[i+1] = 1.0; sy[i+1] = 0.0; sz[i+1] = 0.0;
    }

    dxx[i] += xx[i+1];
    dxy[i] += xy[i+1];
    dxz[i] += xz[i+1];

// fill out curve two

   if (i > 0) {
     yx[3*(i-1)] = xx[i-1]; 
     yy[3*(i-1)] = xy[i-1];
     yz[3*(i-1)] = xz[i-1];
	 crossr(stx, sty, stz, sx[i-1], sy[i-1], sz[i-1], dxx[i-1], dxy[i-1], dxz[i-1]);
	 crossr(px, py, pz, dxx[i-1], dxy[i-1], dxz[i-1], stx, sty, stz);
	 yx[3*(i-1)] += px/10.0;
	 yy[3*(i-1)] += py/10.0;
	 yz[3*(i-1)] += pz/10.0;

     yx[3*(i-1)+1] = xx[i]; 
     yy[3*(i-1)+1] = xy[i];
     yz[3*(i-1)+1] = xz[i];
	 crossr(stx, sty, stz, sx[i], sy[i], sz[i], dxx[i-1], dxy[i-1], dxz[i-1]);
	 crossr(px, py, pz, dxx[i-1], dxy[i-1], dxz[i-1], stx, sty, stz);
	 normalize(px, py, pz);
     yx[3*(i-1)+1] += px/10.0;
     yy[3*(i-1)+1] += py/10.0;
     yz[3*(i-1)+1] += pz/10.0;

     yx[3*(i-1)+2] = xx[i]; 
     yy[3*(i-1)+2] = xy[i];
     yz[3*(i-1)+2] = xz[i];
	 crossr(bx, by, bz, dxx[i-1], dxy[i-1], dxz[i-1], dxx[i], dxy[i], dxz[i]);
	 normalize(bx, by, bz);
     yx[3*(i-1)+2] += bx/10.0; 
     yy[3*(i-1)+2] += by/10.0;
     yz[3*(i-1)+2] += bz/10.0;

   }

  }

  yx[3*nsteps-3] = xx[nsteps-1];
  yy[3*nsteps-3] = xy[nsteps-1];
  yz[3*nsteps-3] = xz[nsteps-1];

  crossr(stx, sty, stz, sx[nsteps-1], sy[nsteps-1], sz[nsteps-1], dxx[nsteps-1], dxy[nsteps-1], dxz[nsteps-1]);
  crossr(px, py, pz, dxx[nsteps-1], dxy[nsteps-1], dxz[nsteps-1], stx, sty, stz);

  yx[3*nsteps-3] += px/10.0;
  yy[3*nsteps-3] += py/10.0;
  yz[3*nsteps-3] += pz/10.0;

  yx[3*nsteps-2] = 0.0;
  yy[3*nsteps-2] = 0.0;
  yz[3*nsteps-2] = 0.0;

  crossr(stx, sty, stz, sx[0], sy[0], sz[0], dxx[nsteps-1], dxy[nsteps-1], dxz[nsteps-1]);
  crossr(px, py, pz, dxx[nsteps-1], dxy[nsteps-1], dxz[nsteps-1], stx, sty, stz);

  yx[3*nsteps-2] += px/10.0;
  yy[3*nsteps-2] += py/10.0;
  yz[3*nsteps-2] += pz/10.0;

  yx[3*nsteps-1] = 0.0; 
  yy[3*nsteps-1] = 0.0;
  yz[3*nsteps-1] = 0.0;

  crossr(bx, by, bz, dxx[nsteps-1], dxy[nsteps-1], dxz[nsteps-1], dxx[0], dxy[0], dxz[0]);
  normalize(bx, by, bz);

  yx[3*nsteps-1] += bx/10.0; 
  yy[3*nsteps-1] += by/10.0;
  yz[3*nsteps-1] += bz/10.0;

  yx[3*nsteps] = yx[0];
  yy[3*nsteps] = yy[0];
  yz[3*nsteps] = yz[0];

  for (int i = 0; i < nsteps; i++) {
    dyx[3*i] = yx[3*i+1]-yx[3*i];
    dyy[3*i] = yy[3*i+1]-yy[3*i];
    dyz[3*i] = yz[3*i+1]-yz[3*i];
    dyx[3*i+1] = yx[3*i+2]-yx[3*i+1];
    dyy[3*i+1] = yy[3*i+2]-yy[3*i+1];
    dyz[3*i+1] = yz[3*i+2]-yz[3*i+1];
    dyx[3*i+2] = yx[3*i+3]-yx[3*i+2];
    dyy[3*i+2] = yy[3*i+3]-yy[3*i+2];
    dyz[3*i+2] = yz[3*i+3]-yz[3*i+2];
  }


  FLOAT LkN = 0.0;

  FLOAT t1x, t1y, t1z, t2x, t2y, t2z;
  FLOAT mu, wr, LkL;

  for (int i = 0; i < nsteps; i++) {
    t1x = dxx[i]; t1y = dxy[i]; t1z = dxz[i];
	LkL = 0.0;
    for (int j = 0; j < 3*nsteps; j++) {
	  t2x = dyx[j]; t2y = dyy[j]; t2z = dyz[j];
      mu = dih(t1x, t1y, t1z, 
               (xx[i]+xx[i+1])/2.0-(yx[j]+yx[j+1])/2.0, (xy[i]+xy[i+1])/2.0-(yy[j]+yy[j+1])/2.0, (xz[i]+xz[i+1])/2.0-(yz[j]+yz[j+1])/2.0, t2x, t2y, t2z);

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j], xy[i]-yy[j], xz[i]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j], xy[i+1]-yy[j], xz[i+1]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr); 
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j+1], xy[i]-yy[j+1], xz[i]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j+1], xy[i+1]-yy[j+1], xz[i+1]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;
	}
	LkN += LkL;
  }

  delete [] dxx; delete [] dxy; delete [] dxz;
  delete [] xx; delete [] xy; delete [] xz;
  delete [] dyx; delete [] dyy; delete [] dyz;
  delete [] yx; delete [] yy; delete [] yz;

  return LkN / (4.0*M_PI);

}


FLOAT bdna::calculate_link_new() {
  FLOAT *xx = new FLOAT[nsteps+1];
  FLOAT *xy = new FLOAT[nsteps+1];
  FLOAT *xz = new FLOAT[nsteps+1];
  FLOAT *yx = new FLOAT[nsteps+1];
  FLOAT *yy = new FLOAT[nsteps+1];
  FLOAT *yz = new FLOAT[nsteps+1];

  FLOAT *sx = new FLOAT[nsteps+1];
  FLOAT *sy = new FLOAT[nsteps+1];
  FLOAT *sz = new FLOAT[nsteps+1];

  FLOAT *dxx = new FLOAT[nsteps];
  FLOAT *dxy = new FLOAT[nsteps];
  FLOAT *dxz = new FLOAT[nsteps];
  FLOAT *dyx = new FLOAT[nsteps];
  FLOAT *dyy = new FLOAT[nsteps];
  FLOAT *dyz = new FLOAT[nsteps];

  xx[0] = 0.0; xy[0] = 0.0; xz[0] = 0.0;
  sx[0] = 1.0; sy[0] = 0.0; sz[0] = 0.0;
  

  matrix temp = identity(4);

  for (int i = 0; i < nsteps; i++) {
    dxx[i] = -temp(1,4);
    dxy[i] = -temp(2,4);
    dxz[i] = -temp(3,4);


    temp = temp*calculateW(v[i]);
    xx[i+1] = temp(1,4); xy[i+1] = temp(2,4); xz[i+1] = temp(3,4);

    sx[i+1] = temp(1,1);
    sy[i+1] = temp(2,1);
    sz[i+1] = temp(3,1);
   

    if (i+1 == nsteps) {
      xx[i+1] = 0.0; xy[i+1] = 0.0; xz[i+1] = 0.0;
	  sx[i+1] = 1.0; sy[i+1] = 0.0; sz[i+1] = 0.0;
    }

    dxx[i] += xx[i+1];
    dxy[i] += xy[i+1];
    dxz[i] += xz[i+1];
    dyx[i] += yx[i+1];
    dyy[i] += yy[i+1];
    dyz[i] += yz[i+1];

  }

  for (int i = 0; i < nsteps; i++) {
    FLOAT px, py, pz, ptx, pty, ptz, tsx, tsy, tsz;
	FLOAT bx, by, bz;
    tsx = dxx[i]+dxx[i+1];
    tsy = dxy[i]+dxy[i+1];
    tsz = dxz[i]+dxz[i+1];
    normalize(tsx, tsy, tsz);
    crossr(bx, by, bz, dxx[i], dxy[i], dxz[i], dxx[i+1], dxy[i+1], dxz[i+1]);
    normalize(bx, by, bz);
    
	crossr(ptx, pty, ptz, sx[i+1], sy[i+1], sz[i+1], tsx, tsy, tsz);
    crossr(px, py, pz, tsx, tsy, tsz, ptx, pty, ptz);
	normalize(px, py, pz);

    if (i+1 != nsteps) {
      yx[i] = xx[i+1] + px;
      yy[i] = xy[i+1] + py;
      yz[i] = xz[i+1] + pz;

	} else {
	  yx[i] = xx[0] + px;
	  yy[i] = xy[0] + py;
	  yz[i] = xz[0] + pz;
   }
   if (i > 0) {
     dyx[i-1] = yx[i]-yx[i-1];
     dyy[i-1] = yy[i]-yy[i-1];
     dyz[i-1] = yz[i]-yz[i-1];
   }
 }  

  yx[nsteps] = yx[0];
  yy[nsteps] = yy[0];
  yz[nsteps] = yz[0];

  dyx[nsteps-1] = yx[nsteps]-yx[nsteps-1];
  dyy[nsteps-1] = yy[nsteps]-yy[nsteps-1];
  dyz[nsteps-1] = yz[nsteps]-yz[nsteps-1];


  FLOAT LkN = 0.0;

  FLOAT t1x, t1y, t1z, t2x, t2y, t2z;
  FLOAT mu, wr, LkL;

  for (int i = 0; i < nsteps; i++) {
    t1x = dxx[i]; t1y = dxy[i]; t1z = dxz[i];
	LkL = 0.0;
    for (int j = 0; j < nsteps; j++) {
	  t2x = dyx[j]; t2y = dyy[j]; t2z = dyz[j];
      mu = dih(t1x, t1y, t1z, 
               (xx[i]+xx[i+1])/2.0-(yx[j]+yx[j+1])/2.0, (xy[i]+xy[i+1])/2.0-(yy[j]+yy[j+1])/2.0, (xz[i]+xz[i+1])/2.0-(yz[j]+yz[j+1])/2.0, t2x, t2y, t2z);

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j], xy[i]-yy[j], xz[i]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j], xy[i+1]-yy[j], xz[i+1]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr); 
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j+1], xy[i]-yy[j+1], xz[i]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j+1], xy[i+1]-yy[j+1], xz[i+1]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;
	}
	LkN += LkL;
  }

  delete [] dxx; delete [] dxy; delete [] dxz;
  delete [] xx; delete [] xy; delete [] xz;
  delete [] dyx; delete [] dyy; delete [] dyz;
  delete [] yx; delete [] yy; delete [] yz;


  return LkN / (4.0*M_PI);
}



FLOAT bdna::calculate_link() {
  FLOAT *xx = new FLOAT[nsteps+1];
  FLOAT *xy = new FLOAT[nsteps+1];
  FLOAT *xz = new FLOAT[nsteps+1];
  FLOAT *yx = new FLOAT[nsteps+1];
  FLOAT *yy = new FLOAT[nsteps+1];
  FLOAT *yz = new FLOAT[nsteps+1];

  FLOAT *dxx = new FLOAT[nsteps];
  FLOAT *dxy = new FLOAT[nsteps];
  FLOAT *dxz = new FLOAT[nsteps];
  FLOAT *dyx = new FLOAT[nsteps];
  FLOAT *dyy = new FLOAT[nsteps];
  FLOAT *dyz = new FLOAT[nsteps];

  xx[0] = 0.0; xy[0] = 0.0; xz[0] = 0.0;
  yx[0] = 1.0; yy[0] = 0.0; yz[0] = 0.0;

  matrix temp = identity(4);

  for (int i = 0; i < nsteps; i++) {
    dxx[i] = -temp(1,4);
    dxy[i] = -temp(2,4);
    dxz[i] = -temp(3,4);
    dyx[i] = -temp(1,4)-temp(1,1);
    dyy[i] = -temp(2,4)-temp(2,1);
    dyz[i] = -temp(3,4)-temp(3,1);


    temp = temp*calculateW(v[i]);
    xx[i+1] = temp(1,4); xy[i+1] = temp(2,4); xz[i+1] = temp(3,4);

    yx[i+1] = temp(1,4)+temp(1,1);
    yy[i+1] = temp(2,4)+temp(2,1);
    yz[i+1] = temp(3,4)+temp(3,1);

    if (i+1 == nsteps) {
      xx[i+1] = 0.0; xy[i+1] = 0.0; xz[i+1] = 0.0;
	  yx[i+1] = 1.0; yy[i+1] = 0.0; yz[i+1] = 0.0;
    }

    dxx[i] += xx[i+1];
    dxy[i] += xy[i+1];
    dxz[i] += xz[i+1];
    dyx[i] += yx[i+1];
    dyy[i] += yy[i+1];
    dyz[i] += yz[i+1];

  }

  FLOAT LkN = 0.0;

  FLOAT t1x, t1y, t1z, t2x, t2y, t2z;
  FLOAT mu, wr, LkL;

  for (int i = 0; i < nsteps; i++) {
    t1x = dxx[i]; t1y = dxy[i]; t1z = dxz[i];
	LkL = 0.0;
    for (int j = 0; j < nsteps; j++) {
	  t2x = dyx[j]; t2y = dyy[j]; t2z = dyz[j];
      mu = dih(t1x, t1y, t1z, 
               (xx[i]+xx[i+1])/2.0-(yx[j]+yx[j+1])/2.0, (xy[i]+xy[i+1])/2.0-(yy[j]+yy[j+1])/2.0, (xz[i]+xz[i+1])/2.0-(yz[j]+yz[j+1])/2.0, t2x, t2y, t2z);

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j], xy[i]-yy[j], xz[i]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j], xy[i+1]-yy[j], xz[i+1]-yz[j], t2x, t2y, t2z);
      wr = correct(mu, wr); 
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i]-yx[j+1], xy[i]-yy[j+1], xz[i]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL -= wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-yx[j+1], xy[i+1]-yy[j+1], xz[i+1]-yz[j+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      LkL += wr;
	}
	LkN += LkL;
  }

  delete [] dxx; delete [] dxy; delete [] dxz;
  delete [] xx; delete [] xy; delete [] xz;
  delete [] dyx; delete [] dyy; delete [] dyz;
  delete [] yx; delete [] yy; delete [] yz;


  return LkN / (4.0*M_PI);
}

FLOAT bdna::calculate_writhe() {
  FLOAT taul = 0.0,  
    wrest = 0.0,
    wrl;
  FLOAT tau;

  matrix temp = identity(4);

  FLOAT *dx = new FLOAT[nsteps];
  FLOAT *dy = new FLOAT[nsteps];
  FLOAT *dz = new FLOAT[nsteps];
  FLOAT *xx = new FLOAT[nsteps+2];
  FLOAT *xy = new FLOAT[nsteps+2];
  FLOAT *xz = new FLOAT[nsteps+2];

  xx[0] = 0.0; xy[0] = 0.0; xz[0] = 0.0;

    // calculate dXYZ
  for (int i = 0; i < nsteps; i++) {
    dx[i] = -temp(1,4);
    dy[i] = -temp(2,4);
    dz[i] = -temp(3,4);
    temp = temp*calculateW(v[i]);
    xx[i+1] = temp(1,4); xy[i+1] = temp(2,4); xz[i+1] = temp(3,4);
    if (i+1 == nsteps) {
      xx[i+1] = 0.0;
      xy[i+1] = 0.0;
      xz[i+1] = 0.0;
    }
    dx[i] += xx[i+1];
    dy[i] += xy[i+1];
    dz[i] += xz[i+1];
  }


  int j, k;

  // calculate tau integral
  for (int i = 0; i < nsteps; i++) {
    j = (i+1 < nsteps ? i+1 : 0);
    k = (i+2 < nsteps ? i+2 : i+2-nsteps);
    tau = dih(dx[i], dy[i], dz[i], dx[j], dy[j], dz[j], dx[k], dy[k], dz[k]);
    taul -= tau;
  }

  // writhe integral estimate

  FLOAT t1x, t1y, t1z, t2x, t2y, t2z;
  FLOAT mu, wr;

  for (int i = 0; i < nsteps-2; i++) {
    t1x = dx[i]; t1y = dy[i]; t1z = dz[i];
    wrl = 0.0;
    for (int q = i+2; q < nsteps; q++) {
      t2x = dx[q]; t2y = dy[q]; t2z = dz[q];
      mu = dih(t1x, t1y, t1z, 
               (xx[i]+xx[i+1])/2.0-(xx[q]+xx[q+1])/2.0, (xy[i]+xy[i+1])/2.0-(xy[q]+xy[q+1])/2.0, (xz[i]+xz[i+1])/2.0-(xz[q]+xz[q+1])/2.0,
	       t2x, t2y, t2z);

      wr = dih(t1x, t1y, t1z, xx[i]-xx[q], xy[i]-xy[q], xz[i]-xz[q], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl += wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-xx[q], xy[i+1]-xy[q], xz[i+1]-xz[q], t2x, t2y, t2z);
      wr = correct(mu, wr); 
      wrl -= wr;

      wr = dih(t1x, t1y, t1z, xx[i]-xx[q+1], xy[i]-xy[q+1], xz[i]-xz[q+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl -= wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-xx[q+1], xy[i+1]-xy[q+1], xz[i+1]-xz[q+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl += wr;
    }
    wrest += wrl;
  }

  wrest = wrest / (2.0*M_PI);
//  FLOAT taui = taul / (2.0*M_PI);  

//  printf("Writhe Estimate = %lf, Tau_Integral = %lf,  SLk = %lf\n", wrest, taui, SLk);

  delete [] dx;
  delete [] dy;
  delete [] dz;
  delete [] xx;
  delete [] xy;
  delete [] xz;

  return wrest;

}



FLOAT bdna::calculate_writhe_closed() {
  FLOAT taul = 0.0,  
    wrest = 0.0,
    wrl;
  FLOAT tau;

  matrix temp = identity(4);

  FLOAT *dx = new FLOAT[nsteps+1];
  FLOAT *dy = new FLOAT[nsteps+1];
  FLOAT *dz = new FLOAT[nsteps+1];
  FLOAT *xx = new FLOAT[nsteps+2];
  FLOAT *xy = new FLOAT[nsteps+2];
  FLOAT *xz = new FLOAT[nsteps+2];

  xx[0] = 0.0; xy[0] = 0.0; xz[0] = 0.0;

    // calculate dXYZ
  for (int i = 0; i < nsteps; i++) {
    dx[i] = -temp(1,4);
    dy[i] = -temp(2,4);
    dz[i] = -temp(3,4);
    temp = temp*calculateW(v[i]);
    xx[i+1] = temp(1,4); xy[i+1] = temp(2,4); xz[i+1] = temp(3,4);
    dx[i] += xx[i+1];
    dy[i] += xy[i+1];
    dz[i] += xz[i+1];
  }
  xx[nsteps+1] = 0.0; xy[nsteps+1] = 0.0; xz[nsteps+1] = 0.0;
  dx[nsteps] = xx[nsteps+1]-xx[nsteps]; dy[nsteps] = xy[nsteps+1]-xy[nsteps]; dz[nsteps] = xz[nsteps+1]-xz[nsteps]; 


  int j, k;

  // calculate tau integral
  for (int i = 0; i < nsteps; i++) {
    j = (i+1 < nsteps ? i+1 : 0);
    k = (i+2 < nsteps ? i+2 : i+2-nsteps);
    tau = dih(dx[i], dy[i], dz[i], dx[j], dy[j], dz[j], dx[k], dy[k], dz[k]);
    taul -= tau;
  }

  // writhe integral estimate

  FLOAT t1x, t1y, t1z, t2x, t2y, t2z;
  FLOAT mu, wr;

  for (int i = 0; i < nsteps; i++) {
    t1x = dx[i]; t1y = dy[i]; t1z = dz[i];
    wrl = 0.0;
    for (int q = i+2; q < nsteps+2; q++) {
      t2x = dx[q]; t2y = dy[q]; t2z = dz[q];
      mu = dih(t1x, t1y, t1z, 
               (xx[i]+xx[i+1])/2.0-(xx[q]+xx[q+1])/2.0, (xy[i]+xy[i+1])/2.0-(xy[q]+xy[q+1])/2.0, (xz[i]+xz[i+1])/2.0-(xz[q]+xz[q+1])/2.0,
	       t2x, t2y, t2z);

      wr = dih(t1x, t1y, t1z, xx[i]-xx[q], xy[i]-xy[q], xz[i]-xz[q], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl += wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-xx[q], xy[i+1]-xy[q], xz[i+1]-xz[q], t2x, t2y, t2z);
      wr = correct(mu, wr); 
      wrl -= wr;

      wr = dih(t1x, t1y, t1z, xx[i]-xx[q+1], xy[i]-xy[q+1], xz[i]-xz[q+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl -= wr;

      wr = dih(t1x, t1y, t1z, xx[i+1]-xx[q+1], xy[i+1]-xy[q+1], xz[i+1]-xz[q+1], t2x, t2y, t2z);
      wr = correct(mu, wr);
      wrl += wr;
    }
    wrest += wrl;
  }

  wrest = wrest / (2.0*M_PI);
//  FLOAT taui = taul / (2.0*M_PI);  
//  FLOAT SLk = wrest+taui;

//  printf("Writhe Estimate = %lf, Tau_Integral = %lf,  SLk = %lf\n", wrest, taui, SLk);

  delete [] dx;
  delete [] dy;
  delete [] dz;
  delete [] xx;
  delete [] xy;
  delete [] xz;

  return wrest;

}



bdna::bdna(const bdna &b) {
  nsteps = b.nsteps;
  beta = b.beta;
  ion = b.ion;
  v = new matrix[b.nsteps];
  bpstep = new int[b.nsteps];
  //F = b.F;
  //I = b.I;  
  for (int i = 0; i < nsteps; i++) {
    v[i] = b.v[i];
    bpstep[i] = b.bpstep[i];
  }
}

bdna &bdna::operator=(const bdna &b) {
  if (this != &b) {
    delete [] v;
    delete [] bpstep;
    beta = b.beta;
    ion = b.ion;
    nsteps = b.nsteps;
    v = new matrix[b.nsteps];
    bpstep = new int[b.nsteps];   

    for (int i = 0; i < nsteps; i++) {
      v[i] = b.v[i];
      bpstep[i] = b.bpstep[i];
    }
  }
  return *this;
}

void bdna::print3dna(char *filename) {
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("PDB output file %s could not be opened for writing\n", filename);
    return;
  }
  fprintf(f, "  %d base-pairs\n", nsteps+1);
  fprintf(f, "   0  step parameters\n");
  fprintf(f, "      shift   slide    rise    tilt    roll   twist\n");
  fprintf(f, "%c-%c    0.00    0.00    0.00    0.00    0.00    0.00\n", (seq[0] == 'Z' ? 'A' : seq[0]), complement(seq[0]));
  for (int i = 0; i < nsteps; i++) {
	matrix vi = v[i];
    fprintf(f, "%c-%c   %6.4lf %6.4lf %6.4lf %6.4lf %6.4lf %6.4lf\n", (seq[i+1] == 'Z' ? 'A' : seq[i+1]), complement(seq[i+1]), vi(4,1), vi(5,1), vi(6,1), vi(1,1), vi(2,1), vi(3,1));
  }
  fclose(f);
}

void bdna::printpdb(char *filename) {

  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("PDB output file %s could not be opened for writing\n", filename);
    return;
  }

  phosphate = new matrix[2*nsteps];
  for (int i = 0; i < 2*nsteps; i++)
    phosphate[i].setsize(3,1);

  matrix temp = identity(4);
  matrix M;

  FLOAT X1, Y1, Z1, X2, Y2, Z2;

  
  for (int i = 0; i < nsteps; i++) {
 
    M = temp * calculateM(v[i]);

    X1 = M(1,4)-3.0*M(1,1)+8.9*M(1,2)-0.4*M(1,3);
    Y1 = M(2,4)-3.0*M(2,1)+8.9*M(2,2)-0.4*M(2,3);
    Z1 = M(3,4)-3.0*M(3,1)+8.9*M(3,2)-0.4*M(3,3);

    phosphate[2*i].setv(1,1, X1);
    phosphate[2*i].setv(2,1, Y1);
    phosphate[2*i].setv(3,1, Z1);

    X2 = M(1,4)-3.0*M(1,1)-8.9*M(1,2)+0.4*M(1,3);
    Y2 = M(2,4)-3.0*M(2,1)-8.9*M(2,2)+0.4*M(2,3);
    Z2 = M(3,4)-3.0*M(3,1)-8.9*M(3,2)+0.4*M(3,3);

    phosphate[2*i+1].setv(1,1, X2);
    phosphate[2*i+1].setv(2,1, Y2);
    phosphate[2*i+1].setv(3,1, Z2);

    fprintf(f, "ATOM  %5d  N   UNK     1    %8.3lf%8.3lf%8.3lf\n", 2*i+1, phosphate[2*i](1,1), phosphate[2*i](2,1), phosphate[2*i](3,1));
    fprintf(f, "ATOM  %5d  C   UNK     1    %8.3lf%8.3lf%8.3lf\n", 2*i+2, phosphate[2*i+1](1,1), phosphate[2*i+1](2,1), phosphate[2*i+1](3,1));

    temp = temp * calculateW(v[i]);

  }

  fclose(f);

}


void bdna::printconfig(char *filename) {

  FILE *f = fopen(filename, "w");

  if (f == NULL) {
    printf("Output file %s could not be opened for writing\n", filename);
    return;
  }

  fprintf(f, "%d\n", nsteps);

  for (int i = 0; i < nsteps; i++) {
    for (int j = 1; j <= 6; j++) {
      fprintf(f, "%6.3f", (float)v[i](j,1));
      fprintf(f, "  ");
    }
    fprintf(f, "\n");
  }

  fclose(f);

}


void bdna::readconfig(char *filename) {

  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    printf("Input coordinate file %s could not be opened!\n\n", filename);
    exit(0);
  }

  char *p;
  char s[2048];
  int ns;
  do {
    p = fgets(s, 2048, f);
    if (p == NULL) {
      printf("Bad file format in structure file!\n");
      return;
    }
  } while (sscanf(s, "%d", &ns) != 1);
  if (ns != nsteps) {
    printf("Not enough allocated space to read configuration!\n");
    return;
  }
  for (int i = 0; i < nsteps; i++) {
    FLOAT a,b,c,d,e,g;
    do {
      p = fgets(s, 2048, f);
      if (p == NULL) {
	printf("Bad file format in structure file!\n");
        return;
      }
    }
    while (sscanf(s, FLOAT_ESTR6, &a, &b, &c, &d, &e, &g) != 6);
    v[i].setv(1,1,a);
    v[i].setv(2,1,b);
    v[i].setv(3,1,c);
    v[i].setv(4,1,d);
    v[i].setv(5,1,e);
    v[i].setv(6,1,g);
  }
  printf("Read %d step parameters\n", ns);
  fclose(f);
  return;
}
