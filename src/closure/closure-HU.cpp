
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "closure-HU.h"
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

FLOAT Lk0;

void sortlocs(short *l, short nhu) {
  for (int i = 0; i < nhu; i++)
  for (int j = 0; j < i; j++)
    if (l[i] < l[j]) {
      short q = l[i];
      l[i] = l[j];
      l[j] = q;
    }
}

int partition(struct pair *p, int left, int right, int index) {
  int pivot = p[index].a;
  struct pair temp = p[index];
  p[index] = p[right];
  p[right] = temp;
  int sindex = left;
  for (int i = left; i < right; i++) {
    if (p[i].a <= pivot) {
      temp = p[sindex];
      p[sindex] = p[i];
      p[i] = temp;
      sindex++;
    }
  }
  temp = p[right];
  p[right] = p[sindex];
  p[sindex] = temp;
  return sindex;
}

void quick_sort_pairlist(struct pair *p, int left, int right) {
  if (right > left) {
    int pivotindex = (left+right)/2;
    int pivotnewindex = partition(p, left, right, pivotindex);
    quick_sort_pairlist(p, left, pivotnewindex-1);
    quick_sort_pairlist(p, pivotnewindex+1, right);
  }
}

int nopairfile = 0;
int pairfile;

void closure::init_pairfile() {
  char s[40];
  sprintf(s, "pair_%ld", ckey);
  pairfile = open(s, O_WRONLY | O_CREAT | O_TRUNC, 0755); 
  write(pairfile, &ckey, sizeof(long int)); 
}

int closure::open_pairfile() {
  char s[40];
  long int temp;
  struct stat *tbuf;
  sprintf(s, "pair_%ld", ckey);
  stat(s, tbuf);
  printf("%d\n", (tbuf->st_size-4)/12);

  pairfile = open(s, O_RDONLY);
  read(pairfile, &temp, sizeof(long int));
  if (temp != ckey) {
    printf("Fatal error, keys in pair file do not match!\n");
    exit(0);
  }
  return (int)((tbuf->st_size-4)/12);
}


void closure::write_pair(int a, int b, int bound) {
  int p[3];
  p[0] = a;
  p[1] = b;
  p[2] = bound;
  int i = write(pairfile, p, 3*sizeof(int));
}
                                                                                


int closure::read_pair(int &a, int &b, int &bound) {
  int p[3];
  int i = read(pairfile, p, 3*sizeof(int));
  if (i < 3*sizeof(int)) return 0;
  a = p[0];
  b = p[1];
  bound = p[2];
  return 1;
}


void closure::close_pairfile() {
  close(pairfile);
}


struct element_top {
  int n;
  int *list;
};


char *abd(char *s, int d) {
  char *tmp = new char[30];
  sprintf(tmp, "%s_%d", s, d);
  return tmp;
}

int getelement(struct element_top e, int n) {
  return e.list[n];
}


int max = 0;

void addelement(struct element_top &e, int index) {

  e.n++;
  
  if (e.n > max) max = e.n;

  if (e.n == 1) {
    e.list = new int[1];
    e.list[0] = index;
    return;
  }

  int *temp = e.list;

  e.list = new int[e.n];

  memcpy(e.list, temp, sizeof(int)*(e.n-1));

  delete [] temp;

  e.list[e.n-1] = index;

  return;
}


FLOAT randg() {

  static int iset = 0;
  static FLOAT gset;

  FLOAT v1, v2, fac, rsq;

  if (iset == 0) {
    do {
      v1 = 2.0*drand48() - 1.0;
      v2 = 2.0*drand48() - 1.0;
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
  long int A = (long int)getpid()+(long int)time(0);
  ckey = A;
  printf("srand48 key = %ld, not using pair file\n", A);
  srand48(A);
  nopairfile = 1;
}

closure::closure(long int key) {
  ckey = key;
  srand48(key);
  printf("srand48 key = %ld\n", key);
}

closure::~closure() {
  delete [] boundary;
  delete [] z;
  delete [] lambda;
  delete [] factor;
}


void closure::init_HU() {
  char s[20];
  HU = new bdna[4];
  for (int k = 0; k < 4; k++) {
    HU[k].nsteps = 14;
    HU[k].v = new matrix[HU[k].nsteps];
    for (int i = 0; i < 14; i++) HU[k].v[i].setsize(6,1);
    sprintf(s, "HU-%d.dat", k+1);
    HU[k].readconfig(s);
  }
}

void closure::set_HU(FLOAT x) {
  PHU = x;
}

void closure::place_HU(bdna *b, int a) {
  int c = (int)(4*drand48());
  for (int i = 0; i < 14; i++) 
    b->v[i+a] = HU[c].v[i];
}

short *closure::place_HU(bdna *b, short a, short a2, short pref, short &nhu) {
  short *locs;
  locs = place_HU_forward(b, a, a2, pref, nhu);
  return locs;
}

short *closure::place_HU_forward(bdna *b, short a, short a2, short pref, short &nhu) {

  short s[1000];
  short *locations;
  nhu = 0;

  short first = 0;

  for (short x = a; x < a2; x++) {
    if ((drand48() < PHU) && ((pref == 0) || (x-pref < 20) || (x-pref > 22)) && ((nhu == 0) || (x-s[nhu-1] < 20) || (x-s[nhu-1] > 22))) {
      if (first == 0) first = x;
      s[nhu++] = x;
      int c = (int)(4*drand48());
      for (int j = 0; j < 8; j++) {
        b->v[j+x] = HU[c].v[j+3];
      }
      x += 8;
    }
  }

  locations = new short[nhu];

  for (int i = 0; i < nhu; i++) locations[i] = s[i];

  return locations;

}


short *closure::place_HU_reverse(bdna *b, short a, short a2, short pref, short &nhu) {

  short s[1000];
  short *locations;
  nhu = 0;
  
  short first = 0;

  short q = (short)(drand48()*(FLOAT)(a2-a))+a;

  for (short x = q; x >= a; x--) {
    if ((drand48() < PHU) && ((pref == 0) || (x-pref < 20) || (x-pref > 22)) && ((nhu == 0) || (s[nhu-1]-x < 20) || (s[nhu-1]-x > 22))) {
      if (first == 0) first = x;
      s[nhu++] = x;
      int c = (int)(4*drand48());
      for (int j = 0; j < 14; j++) {
        b->v[j+x] = HU[c].v[j];
      }
      x -= 13;
    }
  }

  for (short x = a2-1; x > q; x--) {
    if ((drand48() < PHU) && ((pref == 0) || (x-pref < 20) || (x-pref > 22)) && ((nhu == 0) || (s[nhu-1]-x < 20) || (s[nhu-1]-x > 22)) && (((x > first+14) && ((x > first+22) || (x < first+20))) || (first == 0))) {
      s[nhu++] = x;
      int c = (int)(4*drand48());
      for (int j = 0; j < 14; j++) {
        b->v[j+x] = HU[c].v[j];
      }
      x -= 13;
    }

  }
  
  locations = new short[nhu];
  
  for (int i = 0; i < nhu; i++) locations[i] = s[i];
  
  return locations;

}


void closure::generate_config1() {
  for (int i = 0; i < nsteps/2; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg());
    }
  }  
}


void closure::generate_config2() {
    for (int i = nsteps/2; i < nsteps; i++) {
    for (int j = 1; j <= 6; j++) {
      v[i].setv(j, 1, factor[bpstep[i]](j,1)*randg());
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


void closure::setboundarynum(int n) {
  boundary = new matrix[n];
  nbounds = 0;
}

void closure::setepsilon(FLOAT e) {
  epsilon = e;
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



void closure::reducedcombination() {

  matrix *W0 = new matrix[nbounds];

  int *Ntw = new int[nbounds];

  if (strlen(nhuH) > 0) { 
    nh = new histogram[nbounds];
    for (int i = 0; i < nbounds; i++) nh[i].construct_histogram(20, 1.0);
  }
  if (strlen(positionH) > 0) {
    position = new histogram[nbounds];
    for (int i = 0; i < nbounds; i++) position[i].construct_histogram(nsteps-15, 1.0);
  }
  if (strlen(spacingH) > 0) {
    spacing = new histogram[nbounds];
    for (int i = 0; i < nbounds; i++) spacing[i].construct_histogram(nsteps-15, 1.0);
  }

  for (int q = 0; q < nbounds; q++) {
    W0[q] = invert(calculateW(boundary[q]));
    Ntw[q] = 0;
  }

  matrix W1;

  ncombinations = 0;

  int nZ = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;

  int i;

  int dis = (int)ceil(radi/epsilon);

  element_top *E = new element_top[nZ*nY*nX];
  for (i = 0; i < nZ*nY*nX; i++) {
    E[i].n = 0;
  }

  printf("step 2: binning configurations\n");

  for (i = 0; i < N; i++) {
    W1 = calculateW(WA[i].W);
    FLOAT z = -(W1(1,3)*W1(1,4)+W1(2,3)*W1(2,4)+W1(3,3)*W1(3,4));
    FLOAT y = -(W1(1,2)*W1(1,4)+W1(2,2)*W1(2,4)+W1(3,2)*W1(3,4));
    FLOAT x = -(W1(1,1)*W1(1,4)+W1(2,1)*W1(2,4)+W1(3,1)*W1(3,4));

    addelement(E[int(z/epsilon+nZ/2)+nZ*int(y/epsilon+nY/2)+nZ*nY*int(x/epsilon+nX/2)], i);
  }


  FLOAT rc;

  FLOAT factor = radi*radi;
  FLOAT trfactor = gam;
  FLOAT twfactor = twi;
  long long Pr = 0;
  int Ptr = 0;
  int Ptw = 0;

  matrix W2;
  matrix newW;
  matrix tpt;
  matrix tptp;
  matrix W;

  matrix WM;

  if (nopairfile != 1) init_pairfile();

  matrix *WW = new matrix[23];
  bdna *b = new bdna(nsteps);
  short nhus[23];
  short *locs[23];

  printf("\nstep 3: generating second half segments and combining\n");

  FLOAT avgHU = 0.0;

  for (i = 0; i < N; i++) {

    generate_config2();

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
      v[ff] = tptp;
      b->v[ff] = tptp;
    }

    for (int Q = 0; Q <= 22; Q++) {
      if ((Q < 16) && (Q != 0)) place_HU(b, nsteps/2-Q);
      locs[Q] = place_HU(b, ((Q != 0) && (Q < 16) ? nsteps/2-Q+16 : nsteps/2), nsteps-16, (Q != 0 ? nsteps/2-Q : 0), nhus[Q]);
      if (Q == 0) avgHU += (FLOAT)nhus[Q];
      WW[Q] = identity(4);
      for (int ff = nsteps/2; ff < nsteps; ff++) {
	WW[Q] = WW[Q]*calculateW(b->v[ff]);
	b->v[ff] = v[ff];
      }
    }

    for (int q = 0; q < nbounds; q++) {

      for (int Q = 0; Q <= 22; Q++) {
      WM = WW[Q] * W0[q];

      FLOAT Z = WM(3,4);
      FLOAT Y = WM(2,4);
      FLOAT X = WM(1,4);
      int Zindex = (int)(Z/epsilon+nZ/2);
      int Yindex = (int)(Y/epsilon+nY/2);
      int Xindex = (int)(X/epsilon+nX/2);
   
      for (int m = (Zindex-dis < 0 ? 0: Zindex-dis); m <= (Zindex+dis > nZ-1 ? nZ-1 : Zindex+dis); m++)
      for (int n = (Yindex-dis < 0 ? 0: Yindex-dis); n <= (Yindex+dis > nY-1 ? nY-1 : Yindex+dis); n++)
      for (int o = (Xindex-dis < 0 ? 0: Xindex-dis); o <= (Xindex+dis > nX-1 ? nX-1 : Xindex+dis); o++) {
        if (epsilon*epsilon*((Zindex-m)*(Zindex-m)+(Yindex-n)*(Yindex-n)+(Xindex-o)*(Xindex-o)) < (radi*radi+2.0*epsilon*epsilon))
          for (int l = 0; l < E[m+nZ*n+nZ*nY*o].n; l++) {

            int elementID = getelement(E[m+nZ*n+nZ*nY*o], l);
            W = WA[elementID].W;
            int pos = 0;
            if ((WA[elementID].nhu != 0)) 
	      pos = nsteps/2-WA[elementID].locs[WA[elementID].nhu-1];

	    if (pos > 22) pos = 0;

	    if (pos == Q) {

	    newW = calculateW(W)*WM;
	    tpt = calculatetp(newW);
            rc = tpt(4,1)*tpt(4,1)+tpt(5,1)*tpt(5,1)+tpt(6,1)*tpt(6,1);
            if (rc < factor) {
	      Pr++;
	      if (newW(3,3) > trfactor) {
	        Ptr++;
	        if (fabsf(tpt(3,1)) < twfactor) {
	  	  Ntw[q]++;
//		    int p[3]; p[0] = elementID; p[1] = i; p[2] = q;
//		    write(pairfile, p, 3*sizeof(int));
                  if (nopairfile != 1) write_pair(elementID, i, q);
	//	  printf("Q=%d    (%d, %d)\n", Q, elementID, i);
	          short *lls = new short[WA[elementID].nhu+nhus[Q]];
		
	          for (int i = 0 ; i < WA[elementID].nhu; i++) lls[i] = WA[elementID].locs[i];
                  for (int i = 0 ; i < nhus[Q]; i++) lls[i+WA[elementID].nhu] = locs[Q][i];
	          for (int i = 0; i < nhus[Q]+WA[elementID].nhu; i++) {
		    if (strlen(positionH) > 0) position[q].add_data((FLOAT)(lls[i]+1));
		    if ((strlen(spacingH) > 0) && (i != 0)) spacing[q].add_data((FLOAT)(lls[i]-lls[i-1]));
		  }
		  if (strlen(nhuH) > 0) nh[q].add_data((FLOAT)(nhus[Q]+WA[elementID].nhu));
		 
		}
	      }
	    }

            ncombinations++;

	    }
          }
        }
      }
    }


    for (int Q = 0; Q <= 22; Q++)
      if (nhus[Q] > 0) delete [] locs[Q];
  }

  if (nopairfile != 1) close_pairfile();
  printf("Average HU (2nd half) = %f\n", avgHU/((FLOAT)N));

  FLOAT wr, wtr, wtw;

  printf("\nN_tw (by boundary) = ");

  for (int q = 0; q < nbounds; q++) {
    printf("%d ", Ntw[q]);
    Ptw += Ntw[q];
  }

  printf("\nFractional occupancy (by boundary): ");
  for (int q = 0; q < nbounds; q++) {
    printf("%f ", (FLOAT)Ntw[q]/(FLOAT)Ptw);
  }
  printf("\n\n");

  printf("N_r   = %lld %e\n", Pr, wr = (FLOAT)Pr/((FLOAT)N*(FLOAT)N*4.0*M_PI*sqrt(factor)*sqrt(factor)*sqrt(factor)/3.0));
  printf("N_gam = %d %e\n", Ptr, wtr = (FLOAT)Ptr/((1.0-trfactor)*(FLOAT)Pr));
  printf("N_tw  = %d %e\n", Ptw, wtw = (FLOAT)Ptw*360.0/(4.0*twfactor*M_PI*(FLOAT)Ptr));

  printf("J[M] = %e\n", 4.0*M_PI*10000*wr*wtr*wtw/6.022);

  printf("\n\n%lld combinations done. max bin equal %d\n", ncombinations, max);

  if (strlen(nhuH) > 0) {
    for (int q = 0; q < nbounds; q++) nh[q].printhistogram(abd(nhuH, q));
  }
  if (strlen(positionH) > 0) {
    for (int q = 0; q < nbounds; q++) position[q].printhistogram(abd(positionH,q));
  }
  if (strlen(spacingH) > 0) {
    for (int q = 0; q < nbounds; q++) spacing[q].printhistogram(abd(spacingH,q));
  }
  delete [] nh;
  delete [] position;
  delete [] spacing; 
  delete [] E;

  delete [] W0;

}


void closure::persistence(int N) {

  printf("N = %d, calculating persistence length, averages...\n", N);

  init_HU();

  FLOAT Nhu = 0.0;
  FLOAT energy = 0.0;
  FLOAT P = 0.0;

  matrix tpt;
  matrix temp, temp2;

  int *alocations = new int[nsteps];

  for (int i = 0; i < nsteps; i++) alocations[i] = 0;

/*
  for (int i = 0; i < N; i++) {
    temp = identity(4);
    temp2 = identity(4);
    for (int j = 0; j < nsteps; j++)
      v[j] = I[bpstep[j]];

    short *locations;
    short *locations2;
    short nhu, nhu2;

    locations = place_HU(this, 1, nsteps/2, 0, nhu);
    locations2 = place_HU(this, ((nhu != 0) && (locations[nhu-1] > nsteps/2-16) ? locations[nhu-1]+16 : nsteps/2), nsteps-16, (nhu != 0 ? locations[nhu-1] : 0), nhu2);

    for (int j = 0; j < nsteps; j++) {      
      temp = temp * calculateW(v[j]);      
      P += temp(3,4)-temp2(3,4);      
      temp2 = temp2 * calculateW(v[j]);    
    }

  }

  printf("Persistence length (static) = %f\n", P/(FLOAT)N);
*/
  P = 0.0;
 
  for (int i = 0; i < N; i++) {

    generate_config1();
    generate_config2();

    for (int j = 0; j < nsteps; j++) {
      tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
      v[j] = tpt;
      energy += calculate_elenergy(v[j], j);
    }

    short *locations;  
    short *locations2;
    short nhu, nhu2;

    locations = place_HU(this, 1, nsteps/2, 0, nhu);

    for (int k = 0; k < nhu; k++) {
      alocations[locations[k]]++;
    }

    sortlocs(locations, nhu);

    locations2 = place_HU(this, ((nhu != 0) && (locations[nhu-1] > nsteps/2-16) ? locations[nhu-1]+16 : nsteps/2), nsteps-16, (nhu != 0 ? locations[nhu-1] : 0), nhu2);

    sortlocs(locations2, nhu2);

    for (int k = 0; k < nhu2; k++) {
      alocations[locations2[k]]++;
    }

    Nhu += nhu+nhu2;

    temp = identity(4);
    temp2 = identity(4);

    if (i == 0) {
      for (int ll = 0; ll < nhu; ll++) printf("%d ", locations[ll]);
      for (int ll = 0; ll < nhu2; ll++) printf("%d ", locations2[ll]);
      printf("\n");
    }

    for (int j = 0; j < nsteps; j++) {
      temp = temp * calculateW(v[j]);
      P += temp(3,4)-temp2(3,4);
      temp2 = temp2 * calculateW(v[j]);
    }

  }

  printf("Persistence length (Ang) = %f, Avg HU = %f, Avg elastic energy = %f kT\n", P/(FLOAT)N, Nhu/(FLOAT)N, energy/(FLOAT)N);

  for (int i = 0; i < nsteps; i++) {
    printf("%d %d\n", i, alocations[i]);
  }

//  for (int i = 0; i < nsteps; i++) {
//    printf("N_HU[%d] = %d\n", i, alocations[i]);
//  } 

}



void closure::regeneratesystem(int N1) {

  N = N1;
  init_HU();
                                                                                
  Lk0 = 0.0;
  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/360.0;
  }
  Lk0 = round(Lk0);
  printf("Lk0 = %lf\n", Lk0);
                                                                                

  srand48(ckey);

  int npairs = open_pairfile();

  struct pair *pairlist = new pair[npairs];
  struct pair *pairlist2 = new pair[npairs];
  int *first, *second;

  for (int i = 0; i < npairs; i++) {
    read_pair(pairlist[i].a, pairlist[i].b, pairlist[i].bound);
    pairlist2[i].a = pairlist[i].a;
    pairlist2[i].b = pairlist[i].b;
    pairlist2[i].bound = pairlist[i].bound;
  }

  struct pair temp;

  // sort pairlist2 by a value
/*  for (int i = 0; i < npairs; i++) {
    printf("%d %d\n", pairlist[i].a, pairlist[i].b);
    for (int j = i; j > 0; j--) {
      if (pairlist2[j].a < pairlist2[j-1].a) {
        temp.a = pairlist2[j-1].a;
        temp.b = pairlist2[j-1].b;
        temp.bound = pairlist2[j-1].bound;
        pairlist2[j-1].a = pairlist2[j].a;
        pairlist2[j-1].b = pairlist2[j].b;
        pairlist2[j-1].bound = pairlist2[j].bound;
        pairlist2[j].a = temp.a;
        pairlist2[j].b = temp.b;
        pairlist2[j].bound = temp.bound;
      }
    }
  }
*/

  printf("%d bounds\n", nbounds);

  histogram *hwrithe[nbounds];
  for (int i = 0; i < nbounds; i++) hwrithe[i] = new histogram(1000, 0.005);
  histogram *hlink[nbounds];
  for (int i = 0; i < nbounds; i++) hlink[i] = new histogram(320, 0.25);


  printf("sortin\n");  
  quick_sort_pairlist(pairlist2, 0, npairs-1);  
  printf("done\n");

  int last = -1;
  int nfirst = 0;
  for (int i = 0; i < npairs; i++) {
    if (pairlist2[i].a != last) {
      nfirst++;
      last = pairlist2[i].a;
    }
  }

  first = new int[nfirst];

  int firstindex = 0;
  last = -1;
  for (int i = 0; i < npairs; i++) {
    if (pairlist2[i].a != last) {
      first[firstindex++] = pairlist2[i].a;
      last = pairlist2[i].a;
    }
  }

  
  last = -1;
  int nsecond = 0;
  for (int i = 0; i < npairs; i++) {
    if (pairlist[i].b != last) {
      nsecond++;
      last = pairlist[i].b;
    }
  }

  second = new int[nsecond];

  int secondindex = 0;
  last = -1;

  for (int i = 0; i < npairs; i++) {
      if (pairlist[i].b != last) {
      second[secondindex++] = pairlist[i].b;
      last = pairlist[i].b;
    }
  }

  close_pairfile();

  printf("%d unique first halves, %d unique second halves out of %d combinations\n", nfirst, nsecond, npairs);

  matrix tempW, Wall;
  matrix tpt;
  matrix endth;

  matrix *vf = new matrix[npairs*(nsteps/2)];
  WA = new conf[nfirst];
   

  int currentindex = 0;
  for (int i = 0; i < N; i++) {
    generate_config1();
      for (int j = 0; j < nsteps/2; j++) {
        tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
        v[j] = tpt;
      }

    WA[currentindex].locs = place_HU(this, 1, nsteps/2, 0, WA[currentindex].nhu);

    if (first[currentindex] == i) {
      sortlocs(WA[currentindex].locs, WA[currentindex].nhu);

      for (int j = 0; j < nsteps/2; j++) {
        vf[(nsteps/2)*currentindex+j] = v[j];
      }
      currentindex++;
    }
  }
 
  printf("done part 1\n");

  currentindex = 0;

  bdna *b[23];
  for (int lp = 0; lp < 23; lp++) b[lp] = new bdna(nsteps, 1);
  short nhus[23];
  short *locs[23];

  matrix tptp;

  printf("part 2\n");

  for (int i = 0; i < N; i++) {
    generate_config2();
    for (int j = nsteps/2; j < nsteps; j++) {
      tptp = (z[bpstep[j]]*v[j])+I[bpstep[j]];            
      v[j] = tptp; 
    }     
// ????????
    for (int Q = 0; Q <= 22; Q++) {
      for (int j = nsteps/2; j < nsteps; j++) {
        b[Q]->v[j] = v[j]; 
      }
      if ((Q < 16) && (Q != 0)) place_HU(b[Q], nsteps/2-Q);
      locs[Q] = place_HU(b[Q], ((Q != 0) && (Q < 16) ? nsteps/2-Q+16 : nsteps/2), nsteps-16, (Q != 0 ? nsteps/2-Q : 0), nhus[Q]);
      sortlocs(locs[Q], nhus[Q]);
    }

    if (second[currentindex] == i) {
      currentindex++;
      for (int k = 0; k < npairs; k++) {
        int q = 0;
        if (pairlist[k].b == i) {
          while ((q < nfirst) && (first[q] != pairlist[k].a)) q++;
          if ((q == nfirst) || (first[q] != pairlist[k].a)) {
            printf("something went wrong!\n");
            exit(0);
          }
          int QQ = nsteps/2 - WA[q].locs[WA[q].nhu-1];
	  if (QQ > 22) {
	    QQ = 0;
//	    printf("%d\n", QQ);
      }
//	  printf("%d\n", QQ);
          for (int l = 0; l < nsteps/2; l++) {
            v[l] = vf[(nsteps/2)*q+l];
	  }
	  for (int l = nsteps/2; l < nsteps; l++) {
	    v[l] = b[QQ]->v[l];
	  }
          if (boundary[pairlist[k].bound](1,1) == 0.0) {
            FLOAT writhe = calculate_writhe();
            FLOAT twist = calculate_twist();
//            printf("%d: ", QQ);
//            for (int dd = 0 ; dd < WA[q].nhu; dd++) printf("%d ", WA[q].locs[dd]);
//            for (int dd = 0 ; dd < nhus[QQ]; dd++) printf("%d ", locs[QQ][dd]);
//	    printf("\n");
//	    printf("twist = %lf, writhe = %lf\n", twist, writhe);
            if ((twist+writhe)-Lk0 < -1.01) printconfig("out");
            hlink[pairlist[k].bound]->add_data(twist+writhe);
            hwrithe[pairlist[k].bound]->add_data(writhe);

 	    
//	    if (currentindex == 1) printconfig("out");
          } else {
            FLOAT writhe = calculate_writhe();            
            FLOAT twist = calculate_twist();
	    printf("Lk = %lf", twist+writhe);
            if (fabs(round(twist+writhe)-(twist+writhe)) > 0.001) {
	        printf("BAD TWIST\n");
		calculate_twist(); 
	    }       
            hlink[pairlist[k].bound]->add_data(twist+writhe);
            hwrithe[pairlist[k].bound]->add_data(writhe);
          
          }
        }
      }
    }
  }


  for (int i = 0; i < nbounds; i++) {
    char s[50];
    sprintf(s, "writhe_distribution_%d", i);
    hwrithe[i]->printhistogram_pm(s);
    sprintf(s, "link_distribution_%d", i);
    hlink[i]->printhistogram_pm(s);
  }

  delete [] first;
  delete [] second;
  delete [] vf;
  delete [] pairlist;
  delete [] pairlist2;

}




void closure::solveclosureproblem(int N1, int M1, int reduced) {

  matrix W;

  init_HU();

  WA = new conf[N1];

  //  FLOAT *P = new FLOAT[3*N1*nsteps];  // positions of 2nd -> last base (1st at 0 0 0)
  
  N = N1;

  matrix temp;
  matrix Wall;
  matrix tpt;
  matrix endth;

  stress_free_state(1);

 
  Lk0 = 0.0;
  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/360.0;
  }
  Lk0 = round(Lk0);
  printf("Lk0 = %lf\n", Lk0);

//  for (int i = 0; i < nbounds; i++) 
//    for (int j = 0; j <= 100; j++) p[i][j] = 0;

  printf("\nstep 1: generating %d half-configurations\n", N);

  FLOAT avgHU = 0.0;

  for (int i = 0; i < N; i++) {

    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

    }

    WA[i].locs = place_HU(this, 1, nsteps/2, 0, WA[i].nhu);
    sortlocs(WA[i].locs, WA[i].nhu);

    avgHU += WA[i].nhu;

    temp = identity(4);
    for (int j = 0; j < nsteps/2; j++) {
       temp = temp * calculateW(v[j]);
    }

    endth = calculatetp(temp);

    W = endth;

    WA[i].W = W;

  }

  avgHU /= (FLOAT)N;

  printf("average HU (first half) = %f\n", avgHU);

  if (reduced) {

    reducedcombination();

  } 

  printf("done.\n");

  delete [] WA;
//  delete [] HU;


}


