
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "../translate_tools.h"
#include "../file.h"
#include "closure-parallel.h"
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <sys/resource.h>

FLOAT Lk0 = 0.0;


matrix closure::memW(int index) {
  matrix tpt(6,1);
  tpt.setv(1,1, Wmem[6*index]);
  tpt.setv(2,1, Wmem[6*index+1]);
  tpt.setv(3,1, Wmem[6*index+2]);
  tpt.setv(4,1, Wmem[6*index+3]);
  tpt.setv(5,1, Wmem[6*index+4]);
  tpt.setv(6,1, Wmem[6*index+5]);
  return tpt;
}


void closure::setmemW(matrix tpt, int index) {
  Wmem[6*index] = tpt(1,1);
  Wmem[6*index+1] = tpt(2,1);
  Wmem[6*index+2] = tpt(3,1);
  Wmem[6*index+3] = tpt(4,1);
  Wmem[6*index+4] = tpt(5,1);
  Wmem[6*index+5] = tpt(6,1);
}

void closure::pushkey(unsigned short *tempkey, unsigned short *key) {
  tempkey[0] = key[0];
  tempkey[1] = key[1];
  tempkey[2] = key[2];
}

void closure::popkey(unsigned short *tempkey, unsigned short *key) {
  key[0] = tempkey[0];
  key[1] = tempkey[1];
  key[2] = tempkey[2];
}


element_item *closure::getelement(element_top_himem e, int n) {
  element_item *x = new element_item;
  x->list0 = e.data[n].list0;
  x->list1 = e.data[n].list1;
  x->list2 = e.data[n].list2;
  x->index = e.data[n].index;
  return x;
}

unsigned short *closure::getelement(struct element_top e, int n) {
  unsigned short *x = new unsigned short[3];
  x[0] = e.list0[n];
  x[1] = e.list1[n];
  x[2] = e.list2[n];
  return x;
}


int max = 0;

void closure::addelement(struct element_top &e, unsigned short *index) {
  
  e.n++;
  
  if (e.n > max) max = e.n;

  if (e.n == 1) {
    e.list0 = new unsigned short[1];
    e.list1 = new unsigned short[1];
    e.list2 = new unsigned short[1];
    e.list0[0] = index[0];
    e.list1[0] = index[1];
    e.list2[0] = index[2];
    return;
  }

  unsigned short *temp0 = e.list0;
  unsigned short *temp1 = e.list1;
  unsigned short *temp2 = e.list2;

  e.list0 = new unsigned short[e.n];
  e.list1 = new unsigned short[e.n];
  e.list2 = new unsigned short[e.n];

  memcpy(e.list0, temp0, sizeof(unsigned short)*(e.n-1));
  memcpy(e.list1, temp1, sizeof(unsigned short)*(e.n-1));
  memcpy(e.list2, temp2, sizeof(unsigned short)*(e.n-1));


  delete [] temp0;
  delete [] temp1;
  delete [] temp2;

  e.list0[e.n-1] = index[0];
  e.list1[e.n-1] = index[1];
  e.list2[e.n-1] = index[2];

  return;
}


void closure::addelement(struct element_top_himem &e, unsigned short *keyi, int index) {

  e.n++;
  
  if (e.n > max) max = e.n;

  if (e.n == 1) {
    e.data = new element_item[1];
   
    e.data[0].list0 = keyi[0];
    e.data[0].list1 = keyi[1];
    e.data[0].list2 = keyi[2];
    e.data[0].index = index;
    return;
  }

  element_item *temp = e.data;


//  memcpy(temp, e.data, sizeof(element_item)*(e.n

//  unsigned short *temp0 = e.list0;
//  unsigned short *temp1 = e.list1;
//  unsigned short *temp2 = e.list2;

  e.data = new element_item[e.n];

  memcpy(e.data, temp, sizeof(element_item)*(e.n-1));

  delete [] temp;

  e.data[e.n-1].list0 = keyi[0];
  e.data[e.n-1].list1 = keyi[1];
  e.data[e.n-1].list2 = keyi[2];
  e.data[e.n-1].index = index;


  return;
}



FLOAT randg(unsigned short *ckey, drand48_data *buffer) {

  static int iset = 0;
  static FLOAT gset;

  FLOAT r1, r2, v1, v2, fac, rsq;

  if (iset == 0) {
    do {
	  erand48_r(ckey, buffer, &r1);
	  erand48_r(ckey, buffer, &r2);
      v1 = 2.0*r1 - 1.0;
      v2 = 2.0*r2 - 1.0;
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


matrix ringclosureboundary() {
  matrix r(6,1);
  for (int i = 1; i <= 6; i++)
    r.setv(i,1, 0.0);
  return r;
}


closure::closure() {
  unsigned short A;
  A = (int)getpid()+(int)time(0);
  ckey[0] = A;
  ckey[1] = 0;
  ckey[2] = 0;
  printf("erand48 key = %d %d %d, not using pair file\n", ckey[0], ckey[1], ckey[2]);
  restrictends = 1;
}

closure::closure(unsigned short key[3]) {
  ckey[0] = key[0];
  ckey[1] = key[1];
  ckey[2] = key[2];
  restrictends = 1;
}

closure::~closure() {
  delete [] boundary;
  delete [] z;
  delete [] lambda;
  delete [] factor;
}


void closure::m_generate_config1(bdna *b, unsigned short *key, drand48_data *buffer) {
  for (int i = 0; i < nsteps/2; i++) {
    for (int j = 1; j <= 6; j++) {
      b->v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(key, buffer));
    }
  }  
}


void closure::m_generate_config2(bdna *b, unsigned short *key, drand48_data *buffer) {
    for (int i = nsteps/2; i < nsteps; i++) {
    for (int j = 1; j <= 6; j++) {
      b->v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(key, buffer));
    }
  }  
}



void closure::generate_config1(unsigned short *key) {
  drand48_data *buffer = new struct drand48_data;
  for (int i = 0; i < nsteps/2; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(key, buffer));
    }
  }  
}


void closure::generate_config2(unsigned short *key) {
  drand48_data *buffer = new struct drand48_data;
    for (int i = nsteps/2; i < nsteps; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg(key, buffer));
    }
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


void closure::setmem(int himemory) {
  himem = himemory;
}


void closure::setboundarynum(int n) {
  boundary = new matrix[n];
  nbounds = 0;
}

void closure::setepsilon(FLOAT e) {
  epsilon = e;
}

void closure::setrestrictions(int r) {
  restrictends = r;
}

void closure::setpics(int pics, int sc_pics, FLOAT sc_width) {
  make_pics = pics;
  make_pics_sc = sc_pics;
  sc_factor = sc_width;
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

void closure::m_place_constraints(bdna *b) {
  for (int i = 0; i < nproteins; i++) {
    if (pposition[i] >= 0) {
    for (int j = pposition[i]; j < pposition[i]+proteins[i].nsteps; j++) {
      b->v[j] = proteins[i].v[j-pposition[i]];
    } 
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

void closure::addboundary(matrix bc) {
  boundary[nbounds++] = bc;
  printf("boundary %d: \n", nbounds);
  writematrix(stdout, bc);
  printf("\n");
}


void closure::setparams(FLOAT a, FLOAT b, FLOAT c) {
  radi = a;
  gam = b;
  twi = c;
  printf("r=%f gam=%f twi=%f\n", radi, gam, twi);
  return;
}


void closure::printproteins(char *s) {
  char fn[25];
  int append = 0;
  for (int i = 0; i < nproteins; i++) {
    sprintf(fn, "%s.pdb", pnames[i]);
	matrix Wt = identity(4);
    if (pposition[i] >= 0) for (int j = 0; j < pposition[i]; j++) Wt = Wt*calculateW(v[j]);
	rewrite_pdb(fn, s, identity(4), Wt, append);
	append = 1;
  }
}

void closure::setoverlap(int o) {
  checkoverlap = o;
}

void closure::setrandomsamples(int o) {
  randomcombinations = o;
}

void closure::setbinsize(FLOAT sz) {
  bin_width = sz;
}

int closure::overlap1() {
  FLOAT *x, *y, *z;
  printf("Warning, non-threadsafe function called!\n");
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
  printf("Warning, non-threadsafe function called!\n");
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


int closure::overlap(bdna *b) {
  FLOAT *x, *y, *z;
  x = new FLOAT[nsteps-15];
  y = new FLOAT[nsteps-15];
  z = new FLOAT[nsteps-15];
  matrix temp = identity(4);
  for (int i = 0; i < nsteps; i++) {
	temp = temp * calculateW(b->v[i]);
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
      for (int j = start; j < pposition[i]; j++) W = W * calculateW(b->v[j]);
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
	  if (pDNAoverlap(proteinPDBs[i], Ws[i], pposition[i], proteins[i].nsteps, b, x, y, z)) {
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


int closure::overlap() {
  FLOAT *x, *y, *z;
  x = new FLOAT[nsteps-15];
  y = new FLOAT[nsteps-15];
  z = new FLOAT[nsteps-15];
  printf("Warning, non-threadsafe function called!\n");
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

void closure::setnodes(int n) {
  nnodes=n;
  printf("nnodes = %d\n", n);
}

void closure::reducedcombination(int bkey) {

  matrix *W0 = new matrix[nbounds];

  bdna *b = new bdna(nsteps, 1);

  int *Ntw = new int[nbounds];

  unsigned short key[3];
  key[0] = ckey[0]+1;
  key[1] = ckey[1]+5;
  key[2] = bkey;
  unsigned short tempkey[3];

  struct drand48_data *buffer;
  buffer = (drand48_data *)calloc(1, sizeof(drand48_data));

  printf("thread %d, key = %d %d %d\n", bkey, key[0], key[1], key[2]);

  for (int q = 0; q < nbounds; q++) {
    W0[q] = invert(calculateW(boundary[q]));
    Ntw[q] = 0;
  }

  FLOAT Lk0 = 0.0;

  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/(360.0);
  }

  Lk0 = round(Lk0);
//  printf("Lk_0 = %lf\n", Lk0);

  FILE *df;

  if (make_pics || make_pics_sc) {
	if (make_pics_sc) printf("Supercoiling dLk >= %lf\n", sc_factor);
    char descfile[80];
    sprintf(descfile, "structures-MC%dbp-ID%d-node%d.info", nsteps, ckey[0], bkey);
    df = fopen(descfile, "w");
  }

  FLOAT *avgwrithe = new FLOAT[nbounds];
  FLOAT *avgtwist = new FLOAT[nbounds];

  histogram *hwrithe[nbounds];
  histogram *htwist[nbounds];
  for (int i = 0; i < nbounds; i++) {
	hwrithe[i] = new histogram(lround((FLOAT)nsteps/(30.0*bin_width)), bin_width);
	htwist[i] = new histogram(lround((FLOAT)nsteps/(5.0*bin_width)), bin_width);
	avgwrithe[i] = 0.0;
	avgtwist[i] = 0.0;
  }
  histogram *hlink[nbounds];
  for (int i = 0; i < nbounds; i++) hlink[i] = new histogram(nsteps/5, 1.0);

  long long ncombinations = 0;

  int nZ = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;

  int i;

  int dis = (int)ceil(radi/epsilon);


  FLOAT rc;

  FLOAT factor = radi*radi;
  FLOAT trfactor = gam;
  FLOAT twfactor = twi;
  long long Pr = 0;
  int Ptr = 0;
  int Ptw = 0;

  matrix WW;
  matrix W2;
  matrix newW;
  matrix tpt;
  matrix tptp;
  matrix W;

  matrix WM;

  element_item *elementID_himem;   


  printf("\ngenerating and combining, configurations %d to %d...\n", bkey*(N/nnodes), (bkey+1)*(N/nnodes));

  int twsum = 0;
  int sum = 0;
  int pT = 0;
  char s[40];

  for (i = bkey*(N/nnodes); i < (bkey+1)*(N/nnodes); i++) {

//    if (i % 1000 == 0) printf("step %d\n", i);
	do {
    m_generate_config2(b, key, buffer);

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*b->v[ff])+I[bpstep[ff]];
      b->v[ff] = tptp;
    }

	m_place_constraints(b);

    WW = identity(4);
    for (int ff = nsteps/2; ff < nsteps; ff++) {
      WW = WW * calculateW(b->v[ff]);
    }
	sum++;
	} while ((checkoverlap) && (overlap2()));

	if (i < randomcombinations) {
	  int elementID = (int)((FLOAT)randomcombinations*erand48(key));
	  pushkey(tempkey, key);
      key[0] = corr[elementID].list0;
      key[1] = corr[elementID].list1;
      key[2] = corr[elementID].list2;
      m_generate_config1(b, key, buffer);
      for (int j = 0; j < nsteps/2; j++) {
        tptp = (z[bpstep[j]]*b->v[j])+I[bpstep[j]];
        b->v[j] = tptp;
	  }
	  m_place_constraints(b); 
	  popkey(tempkey, key);
	  if (!overlap(b)) pT++;
	}

    for (int q = 0; q < nbounds; q++) {

      WM = WW * W0[q];

      FLOAT Z = WM(3,4);
      FLOAT Y = WM(2,4);
      FLOAT X = WM(1,4);
      int Zindex = (int)(Z/epsilon+nZ/2);
      int Yindex = (int)(Y/epsilon+nY/2);
      int Xindex = (int)(X/epsilon+nX/2);

      for (int m = (Zindex-dis < 0 ? 0: Zindex-dis); m <= (Zindex+dis > nZ-1 ? nZ-1 : Zindex+dis); m++)
      for (int n = (Yindex-dis < 0 ? 0: Yindex-dis); n <= (Yindex+dis > nY-1 ? nY-1 : Yindex+dis); n++)
      for (int o = (Xindex-dis < 0 ? 0: Xindex-dis); o <= (Xindex+dis > nX-1 ? nX-1 : Xindex+dis); o++) {
        if (epsilon*epsilon*((Zindex-m)*(Zindex-m)+(Yindex-n)*(Yindex-n)+(Xindex-o)*(Xindex-o)) < (radi*radi+3.0*epsilon*epsilon))
        for (int l = 0; l < (himem == 1 ? EM[m+nZ*n+nZ*nY*o].n : E[m+nZ*n+nZ*nY*o].n); l++) {
		if (himem == 0) {
	  	  unsigned short *elementID = getelement(E[m+nZ*n+nZ*nY*o], l);
          pushkey(tempkey, key);
          key[0] = elementID[0];
          key[1] = elementID[1];
          key[2] = elementID[2];
          delete [] elementID;
          m_generate_config1(b, key, buffer);
          for (int j = 0; j < nsteps/2; j++) {
            tptp = (z[bpstep[j]]*b->v[j])+I[bpstep[j]];
            b->v[j] = tptp;
	      }
	  	  m_place_constraints(b); 
          W = identity(4);
          for (int j = 0; j < nsteps/2; j++) W = W*calculateW(b->v[j]);
          popkey(tempkey, key);

	  	  newW = W*WM;
        } else {
	  	  elementID_himem = getelement(EM[m+nZ*n+nZ*nY*o], l);
	      newW = calculateW(memW(elementID_himem->index))*WM;
        }

	    tpt = calculatetp(newW);
        rc = tpt(4,1)*tpt(4,1)+tpt(5,1)*tpt(5,1)+tpt(6,1)*tpt(6,1);
        if ((rc < factor) && ((restrictends == 0) || (newW(3,4) < 0.0))) {
	      Pr++;
	      if (newW(3,3) > trfactor) {
	        Ptr++;
	        if ((fabsf(tpt(3,1)) < twfactor)) {
			twsum++;
            if (himem == 1) {
              pushkey(tempkey, key);
              key[0] = elementID_himem->list0;
              key[1] = elementID_himem->list1;
              key[2] = elementID_himem->list2;
              m_generate_config1(b, key, buffer);
              for (int j = 0; j < nsteps/2; j++) {
                tptp = (z[bpstep[j]]*b->v[j])+I[bpstep[j]];
                b->v[j] = tptp;
	          }
	          m_place_constraints(b); 
//              W = identity(4);
//              for (int j = 0; j < nsteps/2; j++) W = W*calculateW(v[j]);
              popkey(tempkey, key);
			}

            if (!overlap(b)) {
	  	    Ntw[q]++;
		    Nrgt[q]++;
		    FLOAT tw, wr;
		    tw = b->calculate_twist();
		    wr = b->calculate_writhe();
			avgtwist[q] += tw;
			avgwrithe[q] += wr;
  	        hwrithe[q]->add_data(wr);
			htwist[q]->add_data(tw);
	        hlink[q]->add_data(lround(tw+wr));
		    if (make_pics || make_pics_sc) fprintf(df, "%d(B%d) writhe = %lf, twist = %lf, Rg = %lf\n", Ntw[q], q, wr, tw, b->calculateR());

		    if (make_pics || (make_pics_sc && (fabs(Lk0 - tw - wr) >= sc_factor))) {
	            sprintf(s, "structure-MC%dbp-B%d-node%d-%d.dat", nsteps, q, bkey, Ntw[q]);
	            b->print3dna(s);
		   	    sprintf(s, "proteins-MC%dbp-B%d-node%d-%d.pdb", nsteps, q, bkey, Ntw[q]);
		       // b->printproteins(s);
	        }
		  }
		}
	    }
	    }
        ncombinations++;
		if (himem == 1) delete elementID_himem;
        }
      }
    }
  }

  FLOAT wr, wtr, wtw;

  printf("%d discarded second configurations\n", sum-N/nnodes);
  FLOAT randomfraction = 1.0;
  if (randomcombinations > 0) printf("%lf random accepted fraction\n", randomfraction = (FLOAT)pT/(FLOAT)randomcombinations);

  printf("\nN_tw (by boundary) = ");

  for (int q = 0; q < nbounds; q++) {
    printf("%d ", Ntw[q]);
    printf("(<tw> = %lf, <wr> = %lf)   ", avgtwist[q]/(FLOAT)Ntw[q], avgwrithe[q]/(FLOAT)Ntw[q]);
    Ptw += Ntw[q];
  }
  printf("\n%lf circular accepted fraction\n", (FLOAT)Ptw/(FLOAT)twsum);

  printf("\nFractional occupancy (by boundary): ");
  for (int q = 0; q < nbounds; q++) {
    printf("%lf ", (FLOAT)Ntw[q]/(FLOAT)Ptw);
  }
  printf("\n\n");
  wr = (FLOAT)Pr/((FLOAT)N*(FLOAT)N*4.0*M_PI*sqrt(factor)*sqrt(factor)*sqrt(factor)/(3.0*(FLOAT)nnodes));
  if (restrictends) wr *= 2.0;
  printf("N_r   = %lld %le\n", Pr, wr);
  printf("N_gam = %d %le\n", Ptr, wtr = (FLOAT)Ptr/((1.0-trfactor)*(FLOAT)Pr));
  printf("N_tw  = %d %le\n", Ptw, wtw = (FLOAT)Ptw*360.0/(4.0*twfactor*M_PI*(FLOAT)Ptr));

  FLOAT J;
  printf("J[M] = %le\n", J=4.0*M_PI*10000.0*wr*wtr*wtw/6.022);
  if (randomcombinations > 0) {
	printf("J_overlap_corrected[M] = %le\n", J/randomfraction);
  }

  char sJ[50];
  sprintf(sJ, "Jfactor-MC%dbp_ID%d-node%d", nsteps, ckey[0], bkey); 
  FILE *fJ = openfwrite(sJ);

  fprintf(fJ, "J = %le M, total closed configurations = %d\n", J, Ptw);
  fprintf(fJ, "sample size = %d squared, %lld total configurations\n\n", N, (long long)N*(long long)N);

  fprintf(fJ, "Epsilon values: R < %lf Ang, cos(y) > %lf, tau < %lf deg\n", radi, gam, twi);
  fprintf(fJ, "Proteins: %s\n", (nproteins == 0 ? "none" : ""));
  for (int q = 0; q < nproteins; q++) {
    fprintf(fJ, "%s (%d steps), position=%d\n", pnames[q], proteins[q].nsteps, pposition[q]);
  }
  for (int q = 0; q < nbounds; q++) {
    fprintf(fJ, "Boundary %d:\n", q);
    writematrix(fJ, boundary[q]);
  }
  fprintf(fJ, "\nJ by boundary = ");
  for (int q = 0; q < nbounds; q++) {
    fprintf(fJ, "%le ", J*(FLOAT)Ntw[q]/(FLOAT)Ptw);
  }
  fprintf(fJ, "\n");
  fclose(fJ);

  printf("\n\n%lld combinations done. max bin equal %d\n", ncombinations, max);

  for (int i = 0; i < nbounds; i++) {
    if (Ntw[i] > 0) {
      char s[80];
      sprintf(s, "writhe_distribution-MC%dbp_ID%d-node%d_B%d", nsteps, ckey[0], bkey, i);
      hwrithe[i]->printhistogram_mm(s);
	  delete hwrithe[i];
	  sprintf(s, "twist_distribution-MC%dbp_ID%d-node%d_B%d", nsteps, ckey[0], bkey, i);
	  htwist[i]->printhistogram_pm(s);
	  delete htwist[i];
      sprintf(s, "link_distribution-MC%dbp_ID%d-node%d_B%d", nsteps, ckey[0], bkey, i);
      hlink[i]->printhistogram_pi(s);
	  delete hlink[i];
    }
  }

  if (make_pics || make_pics_sc) fclose(df);

  delete [] avgtwist;
  delete [] avgwrithe;
  delete [] Ntw;
  delete [] W0;

}

void closure::writeJ(FILE *out) {
  long long Ntotal = (long long)N*(long long)N;
  FLOAT J;
  for (int i = 0; i < nbounds; i++) {
    fprintf(out, "Boundary %d:\n", i+1);
    writematrix(out, boundary[i]);
    fprintf(out, "\nN(closed) = %d, N(total) = %lld, J = %le M\n", Nrgt[i], Ntotal, J = (3.0*360.0*10000.0/(6.022*radi*radi*radi*(1.0-gam)*2.0*M_PI*twi))*Nrgt[i]/Ntotal);
    fprintf(out, "J[%d] = %le M\ndG[%d] = %lf kT\n\n", i+1, J, i+1, -log(J/(8.0*M_PI*M_PI)));
  }
}


void closure::initE(int N1) {

  N = N1;

  read_constraints();

  int nZ = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;

  if (himem == 1) {
    printf("storing half chain end vectors for speed\n");
    Wmem = new float[6*N];
    EM = new element_top_himem[nZ*nY*nX];
    for (int i = 0; i < nZ*nY*nX; i++) {
      EM[i].n = 0;
    }
  } else {

    E = new element_top[nZ*nY*nX];
    for (int i = 0; i < nZ*nY*nX; i++) {
      E[i].n = 0;
    }
  }
  Nrgt = new int[nbounds];
  for (int i = 0; i < nbounds; i++) Nrgt[i] = 0;
  printf("%d bytes in grid (eps=%lf)\n", 4*nX*nY*nZ, epsilon);
  printf("step 1: generating %d first half chains\n", N);
  pthread_mutex_init(&mutex, NULL);
}


void closure::solveclosureproblem(int node, int offset) {

  matrix W;
  
  unsigned short key[3];

  drand48_data *buffer;

  buffer = (drand48_data *)calloc(1, sizeof(drand48_data));

  printf("Thread %d, starting with %d, generating %d elements\n", node, offset, N/nnodes);

  key[0] = ckey[0];
  key[1] = ckey[1];
  key[2] = node;

  printf("Thread %d, key = %d %d %d\n", node, key[0], key[1], key[2]);

  matrix temp;
  matrix tempW, Wall;
  matrix tpt;
  matrix endth;

//  printf("%d random combinations to take\n", randomcombinations);
  if (randomcombinations > 0) corr = new element_item[randomcombinations];

  Lk0 = 0.0;
  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/360.0;
  }
  Lk0 = round(Lk0);
//  printf("Lk0 = %lf\n", Lk0);

//  printf("\nstep 1: generating %d half-configurations and binning keys\n", N);


  FLOAT avgenergy = 0.0;

  int sum = 0;
  unsigned short Tkey[3];
  FLOAT energy = 0.0;


  int nZ = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;

  bdna *b = new bdna(nsteps, 1);

  for (int i = 0; i < N/nnodes; i++) {

	do {
    Tkey[0] = key[0];
    Tkey[1] = key[1];
    Tkey[2] = key[2];

    temp = identity(4);
    m_generate_config1(b, key, buffer);

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*b->v[j])+I[bpstep[j]];
       b->v[j] = tpt;

    }

    m_place_constraints(b);

	sum++;
	} while ((checkoverlap) && (overlap1()));

	if (i < randomcombinations) {
	  corr[i].list0 = Tkey[0];
	  corr[i].list1 = Tkey[1];
	  corr[i].list2 = Tkey[2];
	}

    for (int j = 0; j < nsteps/2; j++) {
       energy += calculate_elenergy(b->v[j], j);
       tempW = calculateW(b->v[j]);
       temp = temp * tempW;
    }

    FLOAT z = -(temp(1,3)*temp(1,4)+temp(2,3)*temp(2,4)+temp(3,3)*temp(3,4));
    FLOAT y = -(temp(1,2)*temp(1,4)+temp(2,2)*temp(2,4)+temp(3,2)*temp(3,4));
    FLOAT x = -(temp(1,1)*temp(1,4)+temp(2,1)*temp(2,4)+temp(3,1)*temp(3,4));

    if (himem) {
	  pthread_mutex_lock( &mutex );
//		printf("Adding data %d\n", node);
      addelement(EM[int(z/epsilon+nZ/2)+nZ*int(y/epsilon+nY/2)+nZ*nY*int(x/epsilon+nX/2)], Tkey, i+offset);
	  setmemW(calculatetp(temp), i+offset);
	  pthread_mutex_unlock( &mutex );

    } else {
      addelement(E[int(z/epsilon+nZ/2)+nZ*int(y/epsilon+nY/2)+nZ*nY*int(x/epsilon+nX/2)], Tkey);
    }


    avgenergy += energy;
	energy = 0.0;

  }

//  avgenergy /= (FLOAT)N;

  printf("%d discarded configurations\n", sum-N/nnodes);

  printf("average energy = %lf kT, approximate total average ~= %lf kT\n", nnodes*avgenergy/(FLOAT)N, 2.0*nnodes*avgenergy/(FLOAT)N);

  printf("done.\n");

/*  if (randomcombinations > 0) delete [] corr;

  if (himem == 1) {
    delete [] EM;
	delete [] Wmem;
  }
  else delete [] E; */
}



