
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "../translate_tools.h"
#include "closure-HU-end.h"
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

FLOAT Lk0;
/*
void closure::sort_locs(long int *l, short nP) {
  for (int i = 0; i < nP; i++)
  for (int j = 0; j < i; j++)
    if (l[i]/(nstructs*2) < l[j]/(nstructs*2)) {
      short q = l[i];
      l[i] = l[j];
      l[j] = q;
    }
}
*/

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
  //int i = write(pairfile, p, 3*sizeof(int));
  write(pairfile, p, 3*sizeof(int));
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


void closure::init_proteins() {
  FILE *fp = fopen("proteins.dat", "r");
  char buf[80];

  int nt;
  int ns;

  fgets(buf, 80, fp);
  sscanf(buf, "%d %d", &nt, &ns);
  ntypes = nt;
  nstructs = ns;

  if (nt == 0) return;

  nproteins = new int[nt];
  psteps = new int[nt];
//  pnames = new char[nt][20];
//  for (int m = 0; m < nt; m++) pnames[m] = new char[20];
  pstartindex = new int[nt];
  Pp = new FLOAT[nt];
  proteins = new bdna[ns];
//  pfilenames = new char[ns][20];
  structtype = new int[ns];

  int pindex = 0;

  for (int i = 0; i < nt; i++) {
    printf("Structure  %d\n", i+1);
    pstartindex[i] = pindex;
    fgets(buf, 80, fp);
	int x1, y1;
	sscanf(buf, "%s %d %d %lf", pnames[i], &x1, &y1, &Pp[i]);
	nproteins[i] = x1;
	psteps[i] = y1;
    printf("%s\n", pnames[i]);
    if (nproteins[i] != 0) 
	for (int j = 0; j < nproteins[i]; j++) {
      if (pindex == ns) {
		printf("error, too many structures!\n");
	  }
      fgets(buf, 80, fp);
	  sscanf(buf, "%s", pfilenames[pindex]);
      char temp[40];
      strcpy(temp, pfilenames[pindex]);
	  strcat(temp, ".dat");
	  printf("Loading %s\n", temp);
	  proteins[pindex].nsteps = psteps[i];
	  proteins[pindex].v = new matrix[psteps[i]];
	  for (int k = 0; k < psteps[i]; k++) proteins[pindex].v[k].setsize(6,1);
	  proteins[pindex].readconfig(temp);
      structtype[pindex] = j;
	  printf("Loaded %s\n", pfilenames[pindex]);
      pindex++;
    }
  }
  fclose(fp);
}

void closure::set_protein(int n, FLOAT x) { Pp[n] = x; }

void closure::place_protein(bdna *b, short int a, short int type, char reverse) {
  int c = (int)(nproteins[type]*drand48())+pstartindex[type];
  if (reverse) {
    for (int i = 0; i < psteps[type]; i++) {
      b->v[i+a] = proteins[c].v[psteps[type]-i-1];
      b->v[i+a].setv(1,1, -b->v[i+a](1,1));
	  b->v[i+a].setv(4,1, -b->v[i+a](4,1));
    }
  } else
  for (int i = 0; i < psteps[type]; i++) 
    b->v[i+a] = proteins[c].v[i];
}

short *q;
int done = 0;
void closure::init_positions(short a, short a2) {
  if (done == 1) delete [] q;
  q = new short[a2-a];
  done = 1;
  for (int i = 0; i < a2-a; i++) q[i] = a+i;
//  for (int i = 0; i < 10*(a2-a); i++) swap_positions(a, a2);
}

void closure::swap_positions(short a, short a2) {
  short p1 = (short)((a2-a)*drand48());
  short p2 = (short)((a2-a)*drand48());
  short temp = q[p1];
  q[p1] = q[p2];
  q[p2] = temp;
}

void closure::add_sort(long int *loc, int nlocs, long int X) {
  int factor = X/(nstructs*2);
  int pos = 0;
  while ((loc[pos]/(nstructs*2) < factor) && (pos < nlocs)) {
    pos++;
  }
  if (pos == nlocs) {
    loc[pos] = X;
    return;
  }
  long int temp = X;
  long int temp2;
  for (int i = pos; i < nlocs+1; i++) {
    if (i < nlocs) temp2 = loc[i];
	loc[i] = temp;
	temp = temp2;
  }
}

long int *closure::place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &nlocs, short exclude) {
  long int *locs;
  long int s[1000];
  nlocs = 0;

  for (int i = 0; i < b2-b1; i++) {
    int po = q[i];
	int leftn = -100;
	int leftindex = -1;
    int maxsize = 1000;
    for (int j = 0; j < nlocs; j++) {
      if (s[j]/(nstructs*2) < po) {
	    leftn = s[j]; leftindex = j;
	  }
	}
	if ((leftn != -100) || (nlocs > 0)) {
	  if ((leftn == -100) && (nlocs > 0)) {
	    maxsize = s[0]/(nstructs*2)-po;
	  } 
	  else if (leftindex != nlocs-1) {
	    maxsize = s[leftindex+1]/(nstructs*2)-po;
	  } else if (overhang) maxsize = a2-po;
	} else if (overhang) maxsize = a2-po;
	char notplaced = 1;
//	if (po < exclude) notplaced = 0;
    if ((leftn != -100) && (leftn/(nstructs*2)+proteins[leftn%(nstructs*2)%nstructs].nsteps > po)) notplaced = 0;
	for (int j = 0; (notplaced) && (j < ntypes); j++) {
          FLOAT pdice = drand48();
	  if ((pdice < Pp[j]) && (maxsize >= psteps[j])) {
        int pnum = (int)(nproteins[j]*drand48())+pstartindex[j];
		int direction = (int)(2.0*drand48());
	    if (((po >= a) && (po < a2)) && (po >= exclude)) {
		if (direction == 0) {
		  for (int k = 0; k < psteps[j]; k++) b->v[po+k] = proteins[pnum].v[k];
		}
		else {
          for (int k = 0; k < psteps[j]; k++) {
		    b->v[po+k] = proteins[pnum].v[psteps[j]-k-1];
		    b->v[po+k].setv(1,1, -v[po+k](1,1));
		    b->v[po+k].setv(4,1, -v[po+k](4,1));
          }
		}
	    }
		notplaced = 0;
   	    long int newloc = po*(nstructs*2)+nstructs*direction+pnum;
	    add_sort(s, nlocs++, newloc);
	  }
	}
  }

//  swap_positions(b1, b2);

  int nreal = 0;

  locs = new long int[nlocs];
  for (int i = 0; i < nlocs; i++) {
    if ((s[i] >= 0) && ((s[i]/(2*nstructs) >= a) && (s[i]/(2*nstructs) < a2)) && (s[i]/(2*nstructs) >= exclude))
	 locs[nreal++] = s[i];
  }
  nlocs = nreal;
  return locs;
}


long int *closure::place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &nlocs) {
  long int *locs;
  long int s[1000];
  nlocs = 0;

  for (int i = 0; i < b2-b1; i++) {
    int po = q[i];
    int leftn = -100;
    int leftindex = -1;
    int maxsize = 1000;
    for (int j = 0; j < nlocs; j++) {
      if (s[j]/(nstructs*2) < po) {
		leftn = s[j]; leftindex = j;
      }
    }
    if ((leftn != -100) || (nlocs > 0)) {
      if ((leftn == -100) && (nlocs > 0)) {
		maxsize = s[0]/(nstructs*2)-po;
      } else if ((leftindex != nlocs-1) && (leftindex != -1)) {
		maxsize = s[leftindex+1]/(nstructs*2)-po;
      } //else if (overhang) maxsize = a2-po;
    }  //
    char notplaced = 1;
    if ((leftn != -100) && ((leftn >= 0 ? leftn/(nstructs*2) : leftn/(nstructs*2)-1)+proteins[leftn%(nstructs*2)%nstructs].nsteps > po)) notplaced = 0;
    for (int j = 0; (notplaced) && (j < ntypes); j++) {
      FLOAT pdice = drand48();
      if ((pdice < Pp[j]) && (maxsize >= psteps[j])) {
      int pnum = (int)(nproteins[j]*drand48())+pstartindex[j];
      int direction = (int)(2.0*drand48());
      if (((po >= a) && (po < a2)) && (!overhang || (psteps[j]+po <= a2))) {
      if (direction == 0) {
		//printf("Placing forward type %d (%d step, %d total), number %d at %d\n", j+1, psteps[j], nproteins[j], pnum+1, po);
    	  for (int k = 0; k < psteps[j]-0; k++) b->v[po+k] = proteins[pnum].v[k];
      }
      else for (int k = 0; k < psteps[j]-0; k++) {
	    b->v[po+k] = proteins[pnum].v[psteps[j]-k-1];
	    b->v[po+k].setv(1,1, -v[po+k](1,1));
	    b->v[po+k].setv(4,1, -v[po+k](4,1));
      }
	  }
      notplaced = 0;
      long int newloc = po*(nstructs*2)+nstructs*direction+pnum;
	  add_sort(s, nlocs++, newloc);
    }
  }
  }

//  swap_positions(b1, b2);
  int nreal = 0;
                                                                                
  locs = new long int[nlocs];
  for (int i = 0; i < nlocs; i++) {
    if ((s[i] >= 0) && (s[i]/(2*nstructs) >= a) && (s[i]/(2*nstructs) < a2))
     if (!overhang || (proteins[s[i]%(2*nstructs)%nstructs].nsteps+s[i]/(2*nstructs) <= a2)) locs[nreal++] = s[i];
  }
  nlocs = nreal;
  return locs;

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

  int Qtotal = 0;
  for (int nn = 0; nn < ntypes; nn++) {
    Qtotal += 2*(psteps[nn]-1);
  }

  short *post = new short[Qtotal];
  short *typet = new short[Qtotal];
  char *directt = new char[Qtotal];

  int qtt = 0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < psteps[i]-1; j++) {
      post[qtt] = j+1;
      typet[qtt] = i;
	  directt[qtt] = 0;
	  qtt++;
	  post[qtt] = j+1;
	  typet[qtt] = i;
	  directt[qtt] = 1;
	  qtt++;
	}
  }

  matrix *WW = new matrix[Qtotal+1];
  bdna *b = new bdna(nsteps);
  short int *nPs = new short[Qtotal+1];
  long *locs;

  histogram *posh = new histogram(nsteps, 1.0);

  printf("\nstep 3: generating second half segments (%d for each first half) and combining\n", Qtotal+1);

  FLOAT avgP = 0.0;

  init_positions(nsteps/2, nsteps);

  matrix W1H;

  histogram *rad = new histogram(200, 5.0);

  for (i = 0; i < N; i++) {

    generate_config2();

    W1H = identity(4);

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
      v[ff] = tptp;
	  b->v[ff] = v[ff];
    }
    locs = place_proteins(this, nsteps/2, nsteps-14, nsteps/2, nsteps-14, 0, WA[i].nP);

    for (int j = nsteps/2; j < nsteps; j++) W1H = W1H*calculateW(v[j]);

    for (int j = 0; j < N; j++) {
	  matrix WTT = calculateW(WA[j].W)*W1H;
      rad->add_data(sqrt(WTT(1,4)*WTT(1,4)+WTT(2,4)*WTT(2,4)+WTT(3,4)*WTT(3,4)));
    }

  }

  rad->printhistogram_pi("radial-HU");

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

  printf("N_r   = %lld %e\n", Pr, wr = (FLOAT)Pr*2.0/((FLOAT)N*(FLOAT)N*4.0*M_PI*sqrt(factor)*sqrt(factor)*sqrt(factor)/3.0));
  printf("N_gam = %d %e\n", Ptr, wtr = (FLOAT)Ptr/((1.0-trfactor)*(FLOAT)Pr));
  printf("N_tw  = %d %e\n", Ptw, wtw = (FLOAT)Ptw*360.0/(4.0*twfactor*M_PI*(FLOAT)Ptr));

  printf("J[M] = %e\n", 4.0*M_PI*10000*wr*wtr*wtw/6.022);

  printf("\n\n%lld combinations done. max bin equal %d\n", ncombinations, max);

  char posf[60];
  sprintf(posf, "protein_distances-%d-ID%d", nsteps, ckey);
  posh->printhistogram_pi(posf);
  delete posh;


  delete [] E;

  delete [] W0;

}


void closure::persistence(int N) {

  printf("N = %d, calculating persistence length, averages...\n", N);

  init_proteins();

  FLOAT NPS = 0.0;
  FLOAT energy = 0.0;
  FLOAT P = 0.0;

  matrix tpt;
  matrix temp, temp2;

  int *alocations = new int[nsteps];

  for (int i = 0; i < nsteps; i++) alocations[i] = 0;
/*
  init_positions(-20, nsteps+20);

  for (int i = 0; i < N; i++) {
    temp = identity(4);
    temp2 = identity(4);
    for (int j = 0; j < nsteps; j++)
      v[j] = I[bpstep[j]];

    long int *locations;
    short nP;

    locations = place_proteins(this, 0, nsteps, -20, nsteps+20, 1, nP);

    for (int j = 0; j < nsteps; j++) {
      temp = temp * calculateW(v[j]);      
      P += temp(3,4)-temp2(3,4);      
      temp2 = temp2 * calculateW(v[j]);    
    }

    delete [] locations;

  }

  printf("Persistence length (static) = %f\n", P/(FLOAT)N);

  P = 0.0;
*/
  short nPs = 0.0;

//  init_positions(nsteps/2, nsteps);
  for (int i = 0; i < N; i++) {

    long int *locations = place_proteins(this, 0, nsteps/2, -30, nsteps/2+30, 0, nPs);

    for (int k = 0; k < nPs; k++) {
      alocations[locations[k]/(nstructs*2)]++;
    }

  } 

  init_positions(nsteps/2-30, nsteps+30);
  for (int i = 0; i < N; i++) {
    long int *locations = place_proteins(this, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs);

    for (int k = 0; k < nPs; k++) {
	  alocations[locations[k]/(nstructs*2)]++;
    }

  }

  for (int i = 0; i < nsteps; i++) {
    printf("%d %d\n", i, alocations[i]);
  }
/*
  init_positions(-20, nsteps+20);

  for (int i = 0; i < N; i++) {

    generate_config1();
    generate_config2();

    for (int j = 0; j < nsteps; j++) {
      tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
      v[j] = tpt;
      energy += calculate_elenergy(v[j], j);
    }

    long int *locations;  
    short int nPs;

   locations = place_proteins(this, 0, nsteps, -20, nsteps+20, 1, nPs);


    NPS += nPs;

    temp = identity(4);
    temp2 = identity(4);

    if (i < 1000) {
      for (int ll = 0; ll < nPs; ll++) printf("%d ", locations[ll]);
      printf("\n");
    }

    for (int j = 0; j < nsteps; j++) {
      temp = temp * calculateW(v[j]);
      P += temp(3,4)-temp2(3,4);
      temp2 = temp2 * calculateW(v[j]);
    }
    delete [] locations;
  }

  printf("Persistence length (Ang) = %f, Avg proteins = %f, Avg elastic energy = %f kT\n", P/(FLOAT)N, NPS/(FLOAT)N, energy/(FLOAT)N);


//  for (int i = 0; i < nsteps; i++) {
//    printf("N_HU[%d] = %d\n", i, alocations[i]);
//  } 
*/
}



void closure::regeneratesystem(int N1) {

  N = N1;
  init_proteins();
                                                                                
  Lk0 = 0.0;
  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/360.0;
  }
  Lk0 = round(Lk0);
  printf("Lk0 = %lf\n", Lk0);
                                                                                


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

  printf("%d bounds\n", nbounds);

  histogram *hwrithe[nbounds];
  for (int i = 0; i < nbounds; i++) hwrithe[i] = new histogram(1000, 0.005);
  histogram *hlink[nbounds];
  for (int i = 0; i < nbounds; i++) hlink[i] = new histogram(320, 1.0);


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
  WA = new conf[nfirst+1];
   

  int currentindex = 0;
  srand48(ckey);
  init_positions(-30, nsteps/2+30);

  printf("regenerating first halves\n");

  for (int i = 0; i < N; i++) {
    generate_config1();
    for (int j = 0; j < nsteps/2; j++) {
      tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
      v[j] = tpt;
    }

    WA[currentindex].locs = place_proteins(this, 0, nsteps/2, -30, nsteps/2+30, 0, WA[currentindex].nP);

    if (first[currentindex] == i) {
      for (int j = 0; j < nsteps/2; j++) {
        vf[(nsteps/2)*currentindex+j] = v[j];
      }
      currentindex++;
    }
  }
 
  printf("done part 1\n");

  int Qtotal = 0;
  for (int nn = 0; nn < ntypes; nn++) {
    Qtotal += 2*(psteps[nn]-1);
  }

  short *post = new short[Qtotal];
  short *typet = new short[Qtotal];
  char *directt = new char[Qtotal];
 
  int qtt = 0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < psteps[i]-1; j++) {
      post[qtt] = j+1;
      typet[qtt] = i;
      directt[qtt] = 0;
      qtt++;
      post[qtt] = j+1;
      typet[qtt] = i;
      directt[qtt] = 1;
      qtt++;
    }
  }
                                                                       
  init_positions(nsteps/2-30, nsteps+30);

  currentindex = 0;

  bdna *b[Qtotal+1];
  printf("making %d DNAs (for middle split positions)\n", Qtotal+1);
  for (int lp = 0; lp <= Qtotal; lp++) {
    b[lp] = new bdna(nsteps, 1);
  }
  short *nPs = new short[Qtotal+1];
  long int *locs[2000];

  matrix tptp;

  printf("part 2\n");

  int conf = 0;
  int negconf = 0;
  int nweird = 0;
  char s1[40], s2[40];
  char s3[40];

  FILE *fdesc;
  FILE *fsc;

  if (make_pics) {
    sprintf(s3, "structures-ID%d.info", ckey);
    fdesc = fopen(s3, "w");
  }

  if (make_pics_sc) {
	sprintf(s3, "supercoiled-ID%d.info", ckey);
	fsc = fopen(s3, "w");
  }

  long int Nct = 0;

  for (int i = 0; i < N; i++) {
    generate_config2();
    for (int j = nsteps/2; j < nsteps; j++) {
      tptp = (z[bpstep[j]]*v[j])+I[bpstep[j]];            
      v[j] = tptp; 
    }     
// ????????
    for (int Q = 0; Q <= Qtotal; Q++) {
      for (int j = nsteps/2; j < nsteps; j++) {
        b[Q]->v[j] = v[j]; 
      }
      if (Q != 0) {
        place_protein(b[Q], nsteps/2-post[Q-1], typet[Q-1], directt[Q-1]);
        locs[Q] = place_proteins(b[Q], nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs[Q], nsteps/2-post[Q-1]+psteps[typet[Q-1]]);
      } else
        locs[Q] = place_proteins(b[Q], nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs[Q]);
        
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
          int pos = 0;
          if ((WA[q].nP != 0)) {
            int lpos = nsteps/2-WA[q].locs[WA[q].nP-1]/(2*nstructs);
            int stx;
            int lnl = proteins[stx=(WA[q].locs[WA[q].nP-1]%(2*nstructs)%(nstructs))].nsteps;
            if (lnl - lpos > 0) {
            int tpx = 0;
            for (int xx = 0; xx < ntypes; xx++) {
              if (pstartindex[xx] <= stx) tpx = xx;
            }
            int reverse = WA[q].locs[WA[q].nP-1]%(2*nstructs)/nstructs;
            pos += 1;
            for (int xx = 0; xx < tpx; xx++) {
              pos += 2*(psteps[xx]-1);
            }
            pos += 2*(lpos-1)+reverse;
          } else pos = 0;
//	  printf("%d\n", QQ);
	      matrix Ww = identity(4);
          for (int l = 0; l < nsteps/2; l++) {
            v[l] = vf[(nsteps/2)*q+l];
		    Ww = Ww * calculateW(v[l]);
	      }
	      for (int l = nsteps/2; l < nsteps; l++) {
	        v[l] = b[pos]->v[l];
			Ww = Ww * calculateW(v[l]);
	      }
//	      writematrix(stdout, Ww);
          if (boundary[pairlist[k].bound](1,1) == 0.0) {
			printf("Linking number = %lf, Linking number (new) = %lf, Tw+Wr=%lf\n", calculate_link(), calculate_link_new(), calculate_twist()+calculate_writhe());
            FLOAT writhe = calculate_writhe();
            FLOAT twist = calculate_twist();

			calculate_twist2();
			printf("Regular twist = %lf\n", calculate_twist());

		    if (make_pics) {
            sprintf(s1, "proteins-%d.pdb", conf);
//          printf("Q=%d: ", pos);
			int append = 0;
            fprintf(fdesc, "Structure %5d: %ld Tw=%1.4lf, Wr=%1.4lf, %d proteins: ", conf, lround(twist+writhe), twist, writhe, WA[q].nP+nPs[pos]);
		    for (int lm = 0; lm < WA[q].nP+nPs[pos]; lm++) {
			  matrix Wplace = identity(4);
              if (lm < WA[q].nP) {
				sprintf(s2, "%s.pdb", pfilenames[(WA[q].locs[lm]%(2*nstructs))%(nstructs)]);
				int prot = WA[q].locs[lm]%(2*nstructs)%nstructs;
                int place = WA[q].locs[lm]/(2*nstructs);
				int rever = WA[q].locs[lm]%(2*nstructs)/nstructs;
                fprintf(fdesc, "%d(%d%s) ", place, prot, (rever ? "-": "+"));
				if (rever) place += proteins[prot].nsteps;
				for (int xy = 0; xy < place; xy++) Wplace = Wplace * calculateW(v[xy]);
				if (rever) {
				  Wplace.setv(1,2, -Wplace(1,2));
			      Wplace.setv(2,2, -Wplace(2,2));
                  Wplace.setv(3,2, -Wplace(3,2));
                  Wplace.setv(1,3, -Wplace(1,3));
                  Wplace.setv(2,3, -Wplace(2,3));
                  Wplace.setv(3,3, -Wplace(3,3));
			    }
				rewrite_pdb(s2, s1, identity(4), Wplace, append);
			  }
			  else {
				int place = locs[pos][lm-WA[q].nP]/(2*nstructs); 
				int prot = locs[pos][lm-WA[q].nP]%(2*nstructs)%nstructs;
				int rever = locs[pos][lm-WA[q].nP]%(2*nstructs)/nstructs;
                fprintf(fdesc, "%d(%d%s) ", place, prot, (rever ? "-": "+"));
			    sprintf(s2, "%s.pdb", pfilenames[prot]);
				if (rever) place += proteins[prot].nsteps;
                for (int xy = 0; xy < place; xy++) Wplace = Wplace * calculateW(v[xy]);
                if (rever) {
                  Wplace.setv(1,2, -Wplace(1,2));
                  Wplace.setv(2,2, -Wplace(2,2));
                  Wplace.setv(3,2, -Wplace(3,2));
                  Wplace.setv(1,3, -Wplace(1,3));
                  Wplace.setv(2,3, -Wplace(2,3));
                  Wplace.setv(3,3, -Wplace(3,3));
                }
                rewrite_pdb(s2, s1, identity(4), Wplace, append);
			  }
			  append = 1;
            }
//			printf("\n");
//	        printf("twist = %lf, writhe = %lf\n", twist, writhe);
			  fprintf(fdesc, "\n");
  		      sprintf(s1, "structure-%d.dat", conf);
			  print3dna(s1);
		    }

			conf++;

			if ((fabs(twist+writhe-Lk0) > 0.1) && make_pics_sc) {
			  fprintf(fsc, "Supercoiled structure %5d (structure %d), Tw=%1.4lf, Wr=%1.4lf, %d proteins\n", nweird, conf-1, twist, writhe, WA[q].nP+nPs[pos]);
			  sprintf(s1, "supercoiled-%d.dat", nweird++);
			  print3dna(s1);
			}

			if (fabs(twist+writhe-round(twist+writhe)) > 0.1) { 
			  calculate_twist_print();
			} else if (!overlap()) {
			  Nct++;
              hlink[pairlist[k].bound]->add_data(lround(twist+writhe));
              hwrithe[pairlist[k].bound]->add_data(writhe);
			}
 	    
//	    if (currentindex == 1) printconfig("out");
        } else {
		    if (make_pics) {
            sprintf(s1, "proteins-%d-bound-%d.pdb", conf, pairlist[k].bound);
//          printf("Q=%d: ", pos);
			int append = 0;
		    for (int lm = 0; lm < WA[q].nP+nPs[pos]; lm++) {
			  matrix Wplace = identity(4);
              if (lm < WA[q].nP) {
				sprintf(s2, "%s.pdb", pfilenames[(WA[q].locs[lm]%(2*nstructs))%(nstructs)]);
				int prot = WA[q].locs[lm]%(2*nstructs)%nstructs;
                int place = WA[q].locs[lm]/(2*nstructs);
				int rever = WA[q].locs[lm]%(2*nstructs)/nstructs;
				if (rever) place += proteins[prot].nsteps;
				for (int xy = 0; xy < place; xy++) Wplace = Wplace * calculateW(v[xy]);
				if (rever) {
				  Wplace.setv(1,2, -Wplace(1,2));
			      Wplace.setv(2,2, -Wplace(2,2));
                  Wplace.setv(3,2, -Wplace(3,2));
                  Wplace.setv(1,3, -Wplace(1,3));
                  Wplace.setv(2,3, -Wplace(2,3));
                  Wplace.setv(3,3, -Wplace(3,3));
			    }
				rewrite_pdb(s2, s1, identity(4), Wplace, append);
			  }
			  else {
				int place = locs[pos][lm-WA[q].nP]/(2*nstructs); 
				int prot = locs[pos][lm-WA[q].nP]%(2*nstructs)%nstructs;
				int rever = locs[pos][lm-WA[q].nP]%(2*nstructs)/nstructs;
			    sprintf(s2, "%s.pdb", pfilenames[prot]);
				if (rever) place += proteins[prot].nsteps;
                for (int xy = 0; xy < place; xy++) Wplace = Wplace * calculateW(v[xy]);
                if (rever) {
                  Wplace.setv(1,2, -Wplace(1,2));
                  Wplace.setv(2,2, -Wplace(2,2));
                  Wplace.setv(3,2, -Wplace(3,2));
                  Wplace.setv(1,3, -Wplace(1,3));
                  Wplace.setv(2,3, -Wplace(2,3));
                  Wplace.setv(3,3, -Wplace(3,3));
                }
                rewrite_pdb(s2, s1, identity(4), Wplace, append);
			  }
			  append = 1;
            }
//			printf("\n");
//	        printf("twist = %lf, writhe = %lf\n", twist, writhe);
	
  		      sprintf(s1, "structure-%d-bound-%d.dat", conf, pairlist[k].bound);
			  print3dna(s1);
		    }
			conf++;
	    } 
      }
    }
    }
    }          
  }

  if (make_pics) fclose(fdesc);
  if (make_pics_sc) fclose(fsc);

  printf("Overlap factor = %lf\n", (FLOAT)Nct/npairs);

  for (int i = 0; i < nbounds; i++) {
    char s[50];
    sprintf(s, "writhe_distribution_%d", i);
    hwrithe[i]->printhistogram_pm(s);
    sprintf(s, "link_distribution_%d", i);
    hlink[i]->printhistogram_pi(s);
  }

  delete [] first;
  delete [] second;
  delete [] vf;
  delete [] pairlist;
  delete [] pairlist2;

}




void closure::solveclosureproblem(int N1, int M1, int reduced) {

  matrix W;

  printf("Initializing proteins\n");
  init_proteins();
  printf("Done\n");

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

  FLOAT avgP = 0.0;
  srand48(ckey);
  init_positions(0, nsteps/2-14); 
  for (int i = 0; i < N; i++) {

    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

    }

    WA[i].locs = place_proteins(this, 0, nsteps/2-14, 0, nsteps/2-14, 0, WA[i].nP);
//    sortlocs(WA[i].locs, WA[i].nP);

    avgP += WA[i].nP;

    temp = identity(4);
    for (int j = 0; j < nsteps/2; j++) {
       temp = temp * calculateW(v[j]);
    }

    endth = calculatetp(temp);

    W = endth;

    WA[i].W = W;

  }

  avgP /= (FLOAT)N;

  printf("average protein (first half) = %f\n", avgP);

  if (reduced) {

    reducedcombination();

  } 

  printf("done.\n");

  delete [] WA;
//  delete [] HU;


}


