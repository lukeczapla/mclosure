
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "../translate_tools.h"
#include "../file.h"
#include "closure-radial.h"
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <sys/resource.h>

FLOAT Lk0 = 0.0;


FLOAT randg(unsigned short *ckey) {

  static int iset = 0;
  static FLOAT gset;

  FLOAT v1, v2, fac, rsq;

  if (iset == 0) {
    do {
      v1 = 2.0*erand48(ckey) - 1.0;
      v2 = 2.0*erand48(ckey) - 1.0;
      rsq = v1*v1+v2*v2;
    } while ((rsq >= 1.0) || (rsq == 0));
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }

}


closure::closure() {
  unsigned short A;
  A = (int)getpid()+(int)time(0);
  ckey[0] = A;
  ckey[1] = 0;
  ckey[2] = 0;
  printf("erand48 key = %d %d %d, not using pair file\n", ckey[0], ckey[1], ckey[2]);
  skey = ckey[0];
  beta = 1.0;
}

closure::closure(unsigned short key[3]) {
  ckey[0] = key[0];
  ckey[1] = key[1];
  ckey[2] = key[2];
  skey = key[0];
  beta = 1.0;
}

closure::~closure() {
  delete [] z;
  delete [] lambda;
  delete [] factor;
}


void closure::generate_config1() {
  for (int i = 0; i < nsteps/2; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(ckey));
    }
  }  
}


void closure::generate_config2() {
    for (int i = nsteps/2; i < nsteps; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(ckey));
    }
  }  
}


void closure::printeigenvalues() {
  for (int i = 0 ; i < number_of_bases*number_of_bases; i++) {
    printf("%c%c\n", btranslate(i / number_of_bases), btranslate(i % number_of_bases));
    writematrix(stdout, lambda[i]);
  }
}


void closure::initialize_params() {

  lambda = new matrix[number_of_bases*number_of_bases];

  FLOAT *d = new FLOAT[6];

  z = new matrix[number_of_bases*number_of_bases];

  factor = new matrix[number_of_bases*number_of_bases];

  int i;
  int j;

  for (i = 0; i < number_of_bases*number_of_bases; i++) {

    z[i].setsize(6,6); 

    z[i] = jeigen(F[i], d);


    lambda[i].setsize(6, 1);
    for (j = 1; j <= 6; j++) lambda[i].setv(j, 1, d[j-1]);

    factor[i].setsize(6,1);
    for (j = 1; j <= 6; j++) { 
      if (lambda[i](j,1) > 0.0)
        factor[i].setv(j, 1, sqrt(1/(beta*lambda[i](j,1))));
      else
        factor[i].setv(j, 1, sqrt(-1/(beta*lambda[i](j,1))));
    }

  }

  delete [] d;

}

void closure::read_constraints() {
  FILE *fp = fopen("fixed-proteins.dat", "r");
  if (fp == NULL) {
    printf("No fixed-proteins.dat file, skipping\n");
    nproteins = 0;
    return;
  }
  char buf[80];

  fgets(buf, 80, fp);
  sscanf(buf, "%d", &nproteins);
  printf("%d proteins, reading\n", nproteins);
  if (nproteins == 0) return;
  proteins = new bdna[nproteins];
  pposition = new int[nproteins];
  proteinPDBs = new proteinPDB[nproteins];
  char pname[20];
  int position;
  int psize;
  char temp[25];
  for (int i = 0; i < nproteins; i++) {
    fgets(buf, 80, fp);
    sscanf(buf, "%s %d %d", pname, &psize, &position);
    strcpy(pnames[i], pname);
    if (psize > 0) printf("Protein %d, %s (%d steps) will be placed at position %d\n", i+1, pname, psize, position);
	else printf("Protein %d, %s will be unbound, make sure position is less than zero\n", i+1, pname);
	sprintf(temp, "%s.pdb", pname);
    proteinPDBs[i].readPDB(temp);
    pposition[i] = position;
	if (psize > 0) {
  	  proteins[i].nsteps = psize;
	  proteins[i].v = new matrix[psize];
	  for (int k = 0; k < psize; k++) proteins[i].v[k].setsize(6,1);
	  sprintf(temp, "%s.dat", pname);
      proteins[i].readconfig(temp);
	  printf("Loaded %s - %s DNA step parameters\n", temp, pname);
    }
  }
}

void closure::place_constraints() {
  for (int i = 0; i < nproteins; i++) {
    if (pposition[i] >= 0) {
    for (int j = pposition[i]; j < pposition[i]+proteins[i].nsteps; j++) {
      v[j] = proteins[i].v[j-pposition[i]];
    } 
    }
  }
}


void closure::setbinsize(FLOAT sz) {
  bin_width = sz;
}

int closure::overlap1() {
  FLOAT *x, *y, *z;
  x = new FLOAT[nsteps/2];
  y = new FLOAT[nsteps/2];
  z = new FLOAT[nsteps/2];
  matrix temp = identity(4);
  for (int i = 0; i < nsteps/2; i++) {
	temp = temp * calculateW(v[i]);
    x[i] = temp(1,4);
	y[i] = temp(2,4);
	z[i] = temp(3,4);
	for (int j = i-20; j >= 0; j--)
		if (sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])) < 20.0) {
		  delete [] x;
		  delete [] y;
		  delete [] z;
		  return 1;
		}
  }
 

  if (nproteins == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  matrix *Ws = new matrix[nproteins];

  matrix W = identity(4);
  int start = 0;
  for (int i = 0; i < nproteins; i++) {
	if (pposition[i] >= 0) {
      for (int j = start; j < pposition[i]; j++) W = W * calculateW(v[j]);
  	  start = pposition[i];
	  Ws[i] = W;
	} else Ws[i] = identity(4);
  }

  for (int i = 0; i < nproteins; i++) {
    if ((pposition[i] < nsteps/2) || (pposition[i] < 0)) {
	for (int j = i+1; j < nproteins; j++) {
	  if (((pposition[j] < nsteps/2) || (pposition[j] < 0)) && (proteinoverlap(proteinPDBs[i], proteinPDBs[j], Ws[i], Ws[j]))) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;
	  }
	  if (pDNAoverlap(proteinPDBs[i], Ws[i], pposition[i], proteins[i].nsteps, x, y, z, 0, nsteps/2)) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;	
	  }
	} }
  }

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] Ws;

  return 0;

}

int closure::overlap2() {
  FLOAT *x, *y, *z;
  x = new FLOAT[nsteps/2+1];
  y = new FLOAT[nsteps/2+1];
  z = new FLOAT[nsteps/2+1];
  matrix temp = identity(4);
  for (int i = nsteps/2; i < nsteps; i++) {
	int ti = i-nsteps/2;
	temp = temp * calculateW(v[i]);
    x[ti] = temp(1,4);
	y[ti] = temp(2,4);
	z[ti] = temp(3,4);
	for (int j = ti-20; j >= 0; j--)
		if (sqrt((x[ti]-x[j])*(x[ti]-x[j])+(y[ti]-y[j])*(y[ti]-y[j])+(z[ti]-z[j])*(z[ti]-z[j])) < 20.0) {
		  delete [] x;
		  delete [] y;
		  delete [] z;
		  return 1;
		}
     }

  if (nproteins == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  matrix *Ws = new matrix[nproteins];

  matrix W = identity(4);
  int start = nsteps/2;
  for (int i = 0; i < nproteins; i++) {
	if (pposition[i] >= nsteps/2) {
      for (int j = start; j < pposition[i]; j++) W = W * calculateW(v[j]);
  	  start = pposition[i];
	  Ws[i] = W;
	} else Ws[i] = identity(4);
  }

  for (int i = 0; i < nproteins; i++) {
    if (pposition[i] >= nsteps/2) {
	for (int j = i+1; j < nproteins; j++) {
	  if ((pposition[j] >= nsteps/2) && (proteinoverlap(proteinPDBs[i], proteinPDBs[j], Ws[i], Ws[j]))) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;
	  }
	  if (pDNAoverlap(proteinPDBs[i], Ws[i], pposition[i], proteins[i].nsteps, x, y, z, nsteps/2, nsteps)) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;	
	  }
	} }
  }

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] Ws;

  return 0;

}


int closure::overlap() {
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
		if (sqrt((x[i-10]-x[j])*(x[i-10]-x[j])+(y[i-10]-y[j])*(y[i-10]-y[j])+(z[i-10]-z[j])*(z[i-10]-z[j])) < 20.0) {
		  delete [] x;
		  delete [] y;
		  delete [] z;
		  return 1;
		}
   }
  }

  if (nproteins == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  matrix *Ws = new matrix[nproteins];

  matrix W = identity(4);
  int start = 0;
  for (int i = 0; i < nproteins; i++) {
	if (pposition[i] >= 0) {
      for (int j = start; j < pposition[i]; j++) W = W * calculateW(v[j]);
  	  start = pposition[i];
	  Ws[i] = W;
	} else Ws[i] = identity(4);
  }

  for (int i = 0; i < nproteins; i++) {
    for (int j = i+1; j < nproteins; j++) {
	  if (proteinoverlap(proteinPDBs[i], proteinPDBs[j], Ws[i], Ws[j])) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;
	  }
	  if (pDNAoverlap(proteinPDBs[i], Ws[i], pposition[i], proteins[i].nsteps, this, x, y, z)) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] Ws;
		return 1;	
	  }
	}
  }

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] Ws;

  return 0;

}

void closure::combination(int N, char *rfile) {

  WA = new matrix[N];

  read_constraints();

  FLOAT RL = 3.4*nsteps;

  histogram *rad = new histogram((int)(ceil(RL/bin_width))+10, bin_width);

  matrix tpt;

  matrix temp, tempW;

  printf("Calculating: generating first half-chains\n");

  for (int i = 0; i < N; i++) {
    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {
       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;
    }

    place_constraints();

	temp = identity(4);
    for (int j = 0; j < nsteps/2; j++) {
       tempW = calculateW(v[j]);
       temp = temp * tempW;
    }

	WA[i] = temp;
  }

  matrix temp2;

  printf("Part Two: generating second half-chains and combining\n");

  for (int i = 0; i < N; i++) {
	generate_config2();

    for (int j = nsteps/2; j < nsteps; j++) {
       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;
    }

    place_constraints();

	temp = identity(4);
    for (int j = nsteps/2; j < nsteps; j++) {
       tempW = calculateW(v[j]);
       temp = temp * tempW;
    }

	FLOAT rv;

	for (int j = 0; j < N; j++) {
	  temp2 = WA[j]*temp;
	  rv = sqrt(temp2(1,4)*temp2(1,4)+temp2(2,4)*temp2(2,4)+temp2(3,4)*temp2(3,4));
	  rad->add_data(rv);
	}

  }

  rad->printhistogram_pm(rfile);

  delete rad;
  delete [] WA;

}
