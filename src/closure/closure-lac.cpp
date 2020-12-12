
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "../proteinPDB/proteinPDB.h"
#include "../translate_tools.h"
#include "../file.h"
#include "closure-lac.h"
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

void closure::pushkey() {
  tempkey[0] = ckey[0];
  tempkey[1] = ckey[1];
  tempkey[2] = ckey[2];
}

void closure::popkey() {
  ckey[0] = tempkey[0];
  ckey[1] = tempkey[1];
  ckey[2] = tempkey[2];
}

struct element_top {
  int n;
  int *list;
};


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


matrix ringclosureboundary() {
  matrix r(6,1);
  for (int i = 1; i <= 6; i++)
    r.setv(i,1, 0.0);
  return r;
}


closure::closure() {
  unsigned short A = getpid()+time(0);
  ckey[0] = A;
  ckey[1] = 0;
  ckey[2] = 0;
  printf("key = %d %d %d\n", ckey[0], ckey[1], ckey[2]);
}

closure::closure(unsigned short *key) {
  ckey[0] = key[0];
  ckey[1] = key[1];
  ckey[2] = key[2];  
  printf("key = %d %d %d\n", key[0], key[1], key[2]);
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

  proteinPDBs = new proteinPDB[ns];

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
		printf("error, too many structures (i.e. change the first line of the protein data file)!\n");
	  }
      fgets(buf, 80, fp);
	  sscanf(buf, "%s", pfilenames[pindex]);
      char temp[40];
      char temp2[40];
      strcpy(temp, pfilenames[pindex]);
	  strcat(temp, ".dat");
	  printf("Loading %s\n", temp);
      strcpy(temp2, pfilenames[pindex]);
	  strcat(temp2, ".pdb");
	  proteinPDBs[pindex].readPDB(temp2);

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
  int c = (int)(nproteins[type]*erand48(ckey))+pstartindex[type];
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



void closure::printeigenvalues() {
  for (int i = 0 ; i < number_of_bases*number_of_bases; i++) {
    printf("%c%c\n", btranslate(i / number_of_bases), btranslate(i % number_of_bases));
    writematrix(stdout, lambda[i]);
    writematrix(stdout, F[i]);
    writematrix(stdout, z[i]);
  }
}


short *q;

int done = 0;
void closure::init_positions(short a, short a2) {
  if (done == 1) delete [] q;
  q = new short[a2-a];
  done = 1;
  for (int i = 0; i < a2-a; i++) q[i] = a+i;
  for (int i = 0; i < 10*(a2-a); i++) swap_positions(a, a2);
}

void closure::swap_positions(short a, short a2) {
  short p1 = (short)((a2-a)*erand48(ckey));
  short p2 = (short)((a2-a)*erand48(ckey));
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

  if (ntypes == 0) {
    return locs;
  }
  

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
          FLOAT pdice = erand48(ckey);
	  if ((pdice < Pp[j]) && (maxsize >= psteps[j])) {
        int pnum = (int)(nproteins[j]*erand48(ckey))+pstartindex[j];
		int direction = (int)(2.0*erand48(ckey));
	    if (((po >= a) && (po < a2)) && (po >= exclude)) {
		if (direction == 0) {
		  for (int k = 0; k < psteps[j]; k++) b->v[po+k] = proteins[pnum].v[k];
		}
		else {
          for (int k = 0; k < psteps[j]; k++) {
		    b->v[po+k] = proteins[pnum].v[psteps[j]-k-1];
		    b->v[po+k].setv(1,1, -b->v[po+k](1,1));
		    b->v[po+k].setv(4,1, -b->v[po+k](4,1));
          }
		}
	    }
		notplaced = 0;
   	    long int newloc = po*(nstructs*2)+nstructs*direction+pnum;
	    add_sort(s, nlocs++, newloc);
	  }
	}
  }

  swap_positions(b1, b2);

  int nreal = 0;

  locs = new long int[nlocs];
  for (int i = 0; i < nlocs; i++) {
    if ((s[i] >= 0) && ((s[i]/(2*nstructs) >= a) && (s[i]/(2*nstructs) < a2)) && (s[i]/(2*nstructs) >= exclude))
	 locs[nreal++] = s[i];
  }
  nlocs = nreal;
  return locs;
}

void closure::place_proteins(bdna *b, long int *locs, int nP) {
  for (int i = 0; i < nP; i++) {
    long int loc = locs[i];
	int pos = loc/(2*nstructs);
	int rev = loc%(2*nstructs)/nstructs;
	int pn = loc%(2*nstructs)%nstructs;
	if (rev) {
	  for (int j = 0; j < proteins[pn].nsteps; j++) {
		b->v[pos+j] = proteins[pn].v[proteins[pn].nsteps-j-1];
	        b->v[pos+j].setv(1,1, -b->v[pos+j](1,1));
		b->v[pos+j].setv(4,1, -b->v[pos+j](4,1)); 
	  }
    } else {
	  for (int j = 0; j < proteins[pn].nsteps; j++) {
		b->v[pos+j] = proteins[pn].v[j];
	  }
	}
  }
}

long int *closure::place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &nlocs) {
  long int *locs;
  long int s[1000];
  nlocs = 0;

  if (ntypes == 0) {
    locs = new long int;
    return locs;
  }


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
      FLOAT pdice = erand48(ckey);
      if ((pdice < Pp[j]) && (maxsize >= psteps[j])) {
      int pnum = (int)(nproteins[j]*erand48(ckey))+pstartindex[j];
      int direction = (int)(2.0*erand48(ckey));
      if (((po >= a) && (po < a2)) && (!overhang || (psteps[j]+po <= a2))) {
      if (direction == 0) {
		//printf("Placing forward type %d (%d step, %d total), number %d at %d\n", j+1, psteps[j], nproteins[j], pnum+1, po);
    	  for (int k = 0; k < psteps[j]; k++) b->v[po+k] = proteins[pnum].v[k];
      }
      else for (int k = 0; k < psteps[j]; k++) {
	    b->v[po+k] = proteins[pnum].v[psteps[j]-k-1];
	    b->v[po+k].setv(1,1, -b->v[po+k](1,1));
	    b->v[po+k].setv(4,1, -b->v[po+k](4,1));
      }
	  }
      notplaced = 0;
      long int newloc = po*(nstructs*2)+nstructs*direction+pnum;
	  add_sort(s, nlocs++, newloc);
    }
  }
  }

  swap_positions(b1, b2);
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

void closure::convert_steps() {
    matrix tpt;
    for (int j = 0; j < nsteps; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

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

void closure::setrandomsamples(int o) {
  randomcombinations = o;
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

  delete [] x;
  delete [] y;
  delete [] z;
  return 0;
}

void closure::illustrate_config(char *dnafile, char *proteinfile, long int *locs, int nP) {
	print3dna(dnafile);
    int append = 0;
    for (int i = 0; i < nP; i++) {
	  matrix Wplace = identity(4);
	  char s2[40];
  	  sprintf(s2, "%s.pdb", pfilenames[locs[i]%(2*nstructs)%(nstructs)]);
	  int prot = locs[i]%(2*nstructs)%nstructs;
      int place = locs[i]/(2*nstructs);
	  int rever = locs[i]%(2*nstructs)/nstructs;
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
				rewrite_pdb(s2, proteinfile, identity(4), Wplace, append);
		append = 1;
    }
}


void closure::setoverlap(int o) {
  checkoverlap = o;
}

int closure::overlap1(long int *locs, int nP) {
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

  if (nP == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  int *pos = new int[nP];
  int *rev = new int[nP];
  int *num = new int[nP];
  matrix *Ws = new matrix[nP];
  for (int i = 0; i < nP; i++) {
	pos[i] = locs[i]/(2*nstructs);
	rev[i] = locs[i]%(2*nstructs)/nstructs;
	num[i] = locs[i]%(2*nstructs)%nstructs;
	Ws[i] = identity(4);
    for (int q = 0; q < (rev[i] ? pos[i]+proteins[num[i]].nsteps : pos[i]); q++) Ws[i] = Ws[i]*calculateW(v[q]);
    if (rev[i]) {
	  Ws[i].setv(1,2, -Ws[i](1,2));
	  Ws[i].setv(2,2, -Ws[i](2,2));
	  Ws[i].setv(3,2, -Ws[i](3,2));
	  Ws[i].setv(1,3, -Ws[i](1,3));
	  Ws[i].setv(2,3, -Ws[i](2,3));
	  Ws[i].setv(3,3, -Ws[i](3,3));
    }
  }

  for (int i = 0; i < nP; i++) {
    for (int j = i+1; j < nP; j++) {
    if (proteinoverlap(proteinPDBs[num[i]], proteinPDBs[num[j]], Ws[i], Ws[j])) { 
//	  printf("PROTEIN overlap!!\n");
//	  illustrate_config(locs, nP);
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1; 
	}
    }
    if (pDNAoverlap(proteinPDBs[num[i]], Ws[i], pos[i], proteins[num[i]].nsteps, x, y, z, 0, nsteps/2)) {
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1;
	}
  }


  delete [] x;
  delete [] y;
  delete [] z;

  delete [] pos;
  delete [] rev;
  delete [] num;
  delete [] Ws;
  return 0;
 

}

int closure::overlap2(long int *locs, int nP) {
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

  if (nP == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  int *pos = new int[nP];
  int *rev = new int[nP];
  int *num = new int[nP];
  matrix *Ws = new matrix[nP];
  for (int i = 0; i < nP; i++) {
	pos[i] = locs[i]/(2*nstructs);
	rev[i] = locs[i]%(2*nstructs)/nstructs;
	num[i] = locs[i]%(2*nstructs)%nstructs;
	Ws[i] = identity(4);
    for (int q = nsteps/2; q < (rev[i] ? pos[i]+proteins[num[i]].nsteps : pos[i]); q++) Ws[i] = Ws[i]*calculateW(v[q]);
    if (rev[i]) {
	  Ws[i].setv(1,2, -Ws[i](1,2));
	  Ws[i].setv(2,2, -Ws[i](2,2));
	  Ws[i].setv(3,2, -Ws[i](3,2));
	  Ws[i].setv(1,3, -Ws[i](1,3));
	  Ws[i].setv(2,3, -Ws[i](2,3));
	  Ws[i].setv(3,3, -Ws[i](3,3));
    }
  }

  for (int i = 0; i < nP; i++) {
    for (int j = i+1; j < nP; j++) {
    if (proteinoverlap(proteinPDBs[num[i]], proteinPDBs[num[j]], Ws[i], Ws[j])) { 
//	  printf("PROTEIN overlap!!\n");
//	  illustrate_config(locs, nP);
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1; 
	}
    }
    if (pDNAoverlap(proteinPDBs[num[i]], Ws[i], pos[i], proteins[num[i]].nsteps, x, y, z, 0, nsteps/2)) {
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1;
	}
  }


  delete [] x;
  delete [] y;
  delete [] z;

  delete [] pos;
  delete [] rev;
  delete [] num;
  delete [] Ws;
  return 0;
 
}

int closure::overlap(long int *locs, int nP) {

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

  if (nP == 0) {
	delete [] x;
	delete [] y;
	delete [] z;
	return 0;
  }

  int *pos = new int[nP];
  int *rev = new int[nP];
  int *num = new int[nP];
  matrix *Ws = new matrix[nP];
  for (int i = 0; i < nP; i++) {
	pos[i] = locs[i]/(2*nstructs);
	rev[i] = locs[i]%(2*nstructs)/nstructs;
	num[i] = locs[i]%(2*nstructs)%nstructs;
	Ws[i] = identity(4);
    for (int q = 0; q < (rev[i] ? pos[i]+proteins[num[i]].nsteps : pos[i]); q++) Ws[i] = Ws[i]*calculateW(v[q]);
    if (rev[i]) {
	  Ws[i].setv(1,2, -Ws[i](1,2));
	  Ws[i].setv(2,2, -Ws[i](2,2));
	  Ws[i].setv(3,2, -Ws[i](3,2));
	  Ws[i].setv(1,3, -Ws[i](1,3));
	  Ws[i].setv(2,3, -Ws[i](2,3));
	  Ws[i].setv(3,3, -Ws[i](3,3));
    }
  }

  for (int i = 0; i < nP; i++) {
    for (int j = i+1; j < nP; j++) {
    if (proteinoverlap(proteinPDBs[num[i]], proteinPDBs[num[j]], Ws[i], Ws[j])) { 
//	  printf("PROTEIN overlap!!\n");
//	  illustrate_config(locs, nP);
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1; 
	}
    }
    if (pDNAoverlap(proteinPDBs[num[i]], Ws[i], pos[i], proteins[num[i]].nsteps, this, x, y, z)) {
	  delete [] pos;
	  delete [] rev;
	  delete [] num;
	  delete [] Ws;
	  delete [] x;
	  delete [] y;
	  delete [] z;

	  return 1;
	}
  }


  delete [] x;
  delete [] y;
  delete [] z;

  delete [] pos;
  delete [] rev;
  delete [] num;
  delete [] Ws;
  return 0;
}


void closure::combination() {

  printf("Generating second halves and combining\n");

  int Qtotal = 0;
  for (int nn = 0; nn < ntypes; nn++) {
    Qtotal += 2*(psteps[nn]-1);
  }

  short *post = new short[Qtotal];
  short *typet = new short[Qtotal];
  char *directt = new char[Qtotal];

  long *locs;

  long int ncomb = 0;

  short int nPs;

  matrix tptp;
  bdna *b = new bdna(nsteps);
  

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

  matrix WW;

  histogram *rdist = new histogram(3.4*(FLOAT)nsteps/5+5.0, 5.0);
  extended_histogram eh(3.4*(FLOAT)nsteps/5.0+5.0, 5.0, 50);
  extended_histogram eh2(200, 0.01, 50);
  extended_histogram eh3(180, 1.0, 50);

  init_positions(nsteps/2-30, nsteps+30);

  for (int Q = 0; Q <= Qtotal; Q++) {
	for (int i = 0; i < N; i++) {
      generate_config2();

	  WW = identity(4);

      for (int ff = nsteps/2; ff < nsteps; ff++) {
        tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
        v[ff] = tptp;
  	    b->v[ff] = v[ff];
      }

      if (Q != 0) {
		place_protein(b, nsteps/2-post[Q-1], typet[Q-1], directt[Q-1]);
		locs = place_proteins(b, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs, nsteps/2-post[Q-1]+psteps[typet[Q-1]]);
      } else
		locs = place_proteins(b, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs);

      for (int ff = nsteps/2; ff < nsteps; ff++) {
	    WW = WW*calculateW(b->v[ff]);
      }


      for (int j = 0; j < N; j++) {
          int pos = 0;
          if ((WA[j].nP != 0)) {
	      int lpos = nsteps/2-WA[j].locs[WA[j].nP-1]/(2*nstructs);
	      int stx;
	      int lnl = proteins[stx=(WA[j].locs[WA[j].nP-1]%(2*nstructs)%(nstructs))].nsteps;
	      if (lnl - lpos > 0) {
            int tpx = 0;
	        for (int xx = 0; xx < ntypes; xx++) {
	  	      if (pstartindex[xx] <= stx) tpx = xx;
	        }
	        int reverse = WA[j].locs[WA[j].nP-1]%(2*nstructs)/nstructs;
	        pos += 1;
	        for (int xx = 0; xx < tpx; xx++) {
	          pos += 2*(psteps[xx]-1);
	        }
	        pos += 2*(lpos-1)+reverse;
	      } else pos = 0;
	    }
		if (Q == pos) {
			matrix Wt = calculateW(WA[j].W)*WW;
			rdist->add_data(sqrt(Wt(1,4)*Wt(1,4)+Wt(2,4)*Wt(2,4)+Wt(3,4)*Wt(3,4)));
			eh.add_data(nPs+WA[j].nP, sqrt(Wt(1,4)*Wt(1,4)+Wt(2,4)*Wt(2,4)+Wt(3,4)*Wt(3,4)));
			eh2.add_data(nPs+WA[j].nP, 1.0+Wt(3,3));
			eh3.add_data(nPs+WA[j].nP, fabs(calculatetp(Wt)(3,1)));
			ncomb++;
		}
      }

	}
  }

  printf("%d total combinations\n", ncomb);
  rdist->printhistogram_pm("radial-HU");
  eh.printhistogram_pm("radial-allHU");
  eh2.printhistogram_pm("gamma-allHU");
  eh3.printhistogram_pm("costau-allHU");
}



void closure::reducedcombination() {

  nbounds = 2;

  matrix *W0 = new matrix[nbounds];

  matrix Wr(4,4);

  int *Ntw = new int[nbounds];
  FLOAT *avgRoll = new FLOAT[nbounds];
  read_reference("ref_lac3", &Wr);
  W0[0].setsize(4,4);
  W0[0] = invert(Wr);
  read_reference("ref_lac4", &Wr);
  W0[1].setsize(4,4);
  W0[1] = invert(Wr);

  for (int q = 0; q < nbounds; q++) {
    Ntw[q] = 0;
    avgRoll[q] = 0.0;
  }

  FLOAT Lk0 = 0.0;

  for (int i = 0; i < nsteps; i++) {
    Lk0 += I[bpstep[i]](3,1)/(360.0);
  }

  Lk0 = round(Lk0);

  matrix W1;

  histogram *hwrithe[nbounds];
  for (int i = 0; i < nbounds; i++) hwrithe[i] = new histogram(nsteps/2, 0.05);
  histogram *htwist[nbounds];
  for (int i = 0; i < nbounds; i++) htwist[i] = new histogram(nsteps*3, 0.05);
  histogram *hlink[nbounds];
  for (int i = 0; i < nbounds; i++) hlink[i] = new histogram(nsteps/5, 1.0);
  histogram *hRg[nbounds];
  for (int i = 0; i < nbounds; i++) hRg[i] = new histogram(nsteps, 1.0);

  FILE *df;

  if (make_pics || make_pics_sc) {
    char descfile[80];
    sprintf(descfile, "structures-%dbp-ID%d.info", nsteps, key);
    df = fopen(descfile, "w");
  }

  int *Noverlap = new int[nbounds];
  for (int i = 0; i < nbounds; i++) Noverlap[i] = 0;

  ncombinations = 0;

  int nZ = (int)(((FLOAT)(nsteps)*3.4+200.0)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4+200.0)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4+200.0)/epsilon)+2;

  int i;

  int dis = (int)ceil(radi/epsilon);

  element_top *E = new element_top[nZ*nY*nX];
  for (i = 0; i < nZ*nY*nX; i++) {
    E[i].n = 0;
  }

  printf("step 2: binning configurations, %ld bytes in grid (eps = %lf)\n", 4*nZ*nY*nX, epsilon);

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
  FLOAT trff = 180.0*acos(gam)/M_PI;
  FLOAT twfactor = twi;
  long long Pr = 0;
  int Ptr = 0;
  int Ptw = 0;

  int twsum = 0;

  matrix W2;
  matrix newW;
  matrix tpt;
  matrix tptp;
  matrix W;

  matrix WM;

  FLOAT E1, E2, Et;

  Et = 0.0;

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
  long *locs[1000];

  bdna **Qb = new bdna*[Qtotal+1];
  for (int cc = 0; cc <= Qtotal; cc++) Qb[cc] = new bdna(nsteps);

  extended_histogram *posh[nbounds];
  for (int i = 0; i < nbounds; i++) {
	posh[i] = new extended_histogram(nsteps, 1.0, nsteps/10);
  }
  extended_histogram *posh2[nbounds];
  for (int i = 0; i < nbounds; i++) {
	posh2[i] = new extended_histogram(nsteps, 1.0, nsteps/10);
  }

  histogram *nh[nbounds];
  for (int i = 0; i < nbounds; i++) nh[i] = new histogram(lround(nsteps/8.0), 1.0);

  histogram *hroll[nbounds];
  for (int i = 0; i < nbounds; i++) hroll[i] = new histogram(151.0, 1.0);

  printf("\nstep 3: generating second half segments (%d for each first half) and combining\n", Qtotal+1);

  FLOAT avgP = 0.0;

  init_positions(nsteps/2-30, nsteps+30);

  int pT = 0;
  int N2 = 0;
  for (i = 0; i < N; i++) {


    E2 = 0.0;
    generate_config2();

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
      v[ff] = tptp;
      E2 += calculate_elenergy(v[ff], ff);
	  b->v[ff] = v[ff];
    }

    for (int Q = 0; Q <= Qtotal; Q++) {
      if (Q != 0) {
		place_protein(b, nsteps/2-post[Q-1], typet[Q-1], directt[Q-1]);
		locs[Q] = place_proteins(b, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs[Q], nsteps/2-post[Q-1]+psteps[typet[Q-1]]);
      } else
		locs[Q] = place_proteins(b, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs[Q]);

      if (Q == 0) avgP += (FLOAT)nPs[Q];

      WW[Q] = identity(4);
      for (int ff = nsteps/2; ff < nsteps; ff++) {
	    WW[Q] = WW[Q]*calculateW(b->v[ff]);
	    Qb[Q]->v[ff] = b->v[ff];	
		b->v[ff] = v[ff];
      }
    }

	if (i < randomcombinations) {

	  int elementID = int((FLOAT)N*erand48(ckey));

	  pushkey();

	  ckey[0] = WA[elementID].key0;
	  ckey[1] = WA[elementID].key1;
	  ckey[2] = WA[elementID].key2;

	  generate_config1();
	  for (int j = 0; j < nsteps/2; j++) {
	       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
	       v[j] = tpt;
	  }
   	  place_proteins(this, WA[elementID].locs, WA[elementID].nP);

	  popkey();

      int pos = 0;
      if ((WA[elementID].nP != 0)) {
	  int lpos = nsteps/2-WA[elementID].locs[WA[elementID].nP-1]/(2*nstructs);
	  int stx;
	  int lnl = proteins[stx=(WA[elementID].locs[WA[elementID].nP-1]%(2*nstructs)%(nstructs))].nsteps;
	  if (lnl - lpos > 0) {
        int tpx = 0;
	    for (int xx = 0; xx < ntypes; xx++) {
	      if (pstartindex[xx] <= stx) tpx = xx;
	    }
	    int reverse = WA[elementID].locs[WA[elementID].nP-1]%(2*nstructs)/nstructs;
	    pos += 1;
	    for (int xx = 0; xx < tpx; xx++) {
	      pos += 2*(psteps[xx]-1);
	    }
	    pos += 2*(lpos-1)+reverse;
	  } else pos = 0;
	  }

	  for (int ii = nsteps/2; ii < nsteps; ii++) {
		v[ii] = Qb[pos]->v[ii];
	  }

	  int nPtotal;

	  long int *lls = new long int[nPtotal = WA[elementID].nP+nPs[pos]];
	  for (int ii = 0; ii < WA[elementID].nP; ii++) lls[ii] = WA[elementID].locs[ii];
	  for (int ii = 0; ii < nPs[pos]; ii++) lls[ii+WA[elementID].nP] = locs[pos][ii];

	  if (!overlap(lls, nPtotal)) pT++;

	  delete [] lls;

	}

    for (int q = 0; q < nbounds; q++) {

      for (int Q = 0; Q <= Qtotal; Q++) {
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
          if ((WA[elementID].nP != 0)) {
	      int lpos = nsteps/2-WA[elementID].locs[WA[elementID].nP-1]/(2*nstructs);
	      int stx;
	      int lnl = proteins[stx=(WA[elementID].locs[WA[elementID].nP-1]%(2*nstructs)%(nstructs))].nsteps;
	      if (lnl - lpos > 0) {
            int tpx = 0;
	        for (int xx = 0; xx < ntypes; xx++) {
	  	      if (pstartindex[xx] <= stx) tpx = xx;
	        }
	        int reverse = WA[elementID].locs[WA[elementID].nP-1]%(2*nstructs)/nstructs;
	        pos += 1;
	        for (int xx = 0; xx < tpx; xx++) {
	          pos += 2*(psteps[xx]-1);
	        }
	        pos += 2*(lpos-1)+reverse;
	      } else pos = 0;
	    }

	    if (pos == Q) {

	    newW = calculateW(W)*WM;
	    tpt = calculatetp(newW);
            rc = tpt(4,1)*tpt(4,1)+tpt(5,1)*tpt(5,1)+tpt(6,1)*tpt(6,1);
            if ((rc < factor)) {
	      Pr++;
	      if ((fabsf(tpt(1,1)) < trff) && (tpt(2,1) > 0.0) && (tpt(2,1) < 140.0)) {
	        Ptr++;
	        if (fabsf(tpt(3,1)) < twfactor) {
//			  writematrix(stdout, calculateW(W)*WW[pos]);
			  FLOAT ro = tpt(2,1);
			  int nPtotal;
			  int nspace;
			  twsum++;
			  pushkey();

			  ckey[0] = WA[elementID].key0;
			  ckey[1] = WA[elementID].key1;
			  ckey[2] = WA[elementID].key2;

		      E1 = 0.0;
		      generate_config1();
			  for (int j = 0; j < nsteps/2; j++) {
			       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
			       v[j] = tpt;
				E1 += calculate_elenergy(v[j], j);
			  }
    		  place_proteins(this, WA[elementID].locs, WA[elementID].nP);

			  popkey();
			  matrix Wi = identity(4);
			  for (int ii = 0; ii < nsteps/2; ii++) {
				Wi = Wi * calculateW(v[ii]);
			  }
			  for (int ii = nsteps/2; ii < nsteps; ii++) {
				v[ii] = Qb[pos]->v[ii];
				Wi = Wi * calculateW(v[ii]);
			  }
//			  writematrix(stdout, Wi);

	          long *lls = new long[nPtotal = WA[elementID].nP+nPs[Q]];
		
	          for (int i = 0 ; i < WA[elementID].nP; i++) lls[i] = WA[elementID].locs[i];
              for (int i = 0 ; i < nPs[Q]; i++) lls[i+WA[elementID].nP] = locs[Q][i];

			  if (neglect || (!overlap(lls, nPtotal))) {

		        Ntw[q]++;
			avgRoll[q] += ro;
	  		    long tspace = 0;

	            for (int i = 0; i < nPs[Q]+WA[elementID].nP-1; i++) {
				  posh[q]->add_data(nPtotal, nspace = (lls[i+1]/(2*nstructs))-(lls[i]/(2*nstructs)));
				  posh2[q]->add_data(nPtotal, lls[i]/(2*nstructs));
				  tspace += nspace;
                }
				if (nPtotal > 0) posh2[q]->add_data(nPtotal, lls[nPtotal-1]/(2*nstructs));
		  	  //  if ((nPtotal > 1) && (equals(boundary[q], ringclosureboundary())))
			//	  posh[q]->add_data(nPtotal, nsteps-tspace);
		  	    nh[q]->add_data(WA[elementID].nP+nPs[Q]);
				FLOAT writhe = calculate_writhe();
				FLOAT twist = calculate_twist();
				FLOAT Rg = calculateR();
				hroll[q]->add_data(ro);
				hwrithe[q]->add_data(writhe);
				htwist[q]->add_data(twist);
				hlink[q]->add_data(round(writhe + twist));
				hRg[q]->add_data(Rg);
				if (make_pics || make_pics_sc) {
				  fprintf(df, "%d(B%d) writhe = %lf, twist = %lf, Rg = %lf Ang, opening = %lf deg, %d proteins: ", Ntw[q], q, calculate_writhe(), calculate_twist(), Rg, ro, nPtotal);
				  for (int pp = 0; pp < nPtotal; pp++) {
				    fprintf(df, "%d(%d%c) ", lls[pp]/(2*nstructs), lls[pp]%(2*nstructs)%nstructs, (lls[pp]%(2*nstructs)/nstructs == 0 ? '+' : '-'));
				  }
				  fprintf(df, "\n");
				}
				if ((make_pics) || ((make_pics_sc) && (fabs(Lk0-writhe-twist) > 0.01))) {
					char s1[40], s2[40];
					sprintf(s1, "structure-%dbp-ID%d_B%d-%d.dat", nsteps, key, q, Ntw[q]);
					sprintf(s2, "proteins-%dbp-ID%d_B%d-%d.pdb", nsteps, key, q, Ntw[q]);
					illustrate_config(s1, s2, lls, nPtotal);
				}
			  } else {
				Noverlap[q]++;
			  }
			  delete [] lls;
	    }
	    } 
        }       
            ncombinations++;
	}
  }
  }
 } 
}
  for (int Q = 0; Q <= Qtotal; Q++) delete [] locs[Q];

}

  FLOAT randomfraction = 1.0;
  if (randomcombinations > 0) printf("%lf random accepted fraction\n", randomfraction = (FLOAT)pT/(FLOAT)randomcombinations);

  printf("average protein (2nd half) = %f\n", avgP/((FLOAT)N));

  FLOAT wr, wtr, wtw;

  printf("\nN_tw (by boundary, number of overlaps in parentheses) = ");


  for (int q = 0; q < nbounds; q++) {
    printf("%d (and %d, %lf)  ", Ntw[q], Noverlap[q], avgRoll[q]/(FLOAT)Ntw[q]);
    Ptw += Ntw[q];
  }

  printf("\n%lf circular accepted fraction\n", (FLOAT)Ptw/(FLOAT)twsum);

  printf("\nFractional occupancy (by boundary): ");
  for (int q = 0; q < nbounds; q++) {
    printf("%f ", (FLOAT)Ntw[q]/(FLOAT)Ptw);
  }
  printf("\n\n");

  printf("N_r   = %lld %le\n", Pr, wr = (FLOAT)Pr/((FLOAT)N*(FLOAT)N*4.0*M_PI*sqrt(factor)*sqrt(factor)*sqrt(factor)/3.0));
  printf("N_gam = %d %le\n", Ptr, wtr = (FLOAT)Ptr/(2.0*sin(acos(gam))*140.0/360.0*(FLOAT)Pr));
  printf("N_tw  = %d %le\n", Ptw, wtw = (FLOAT)Ptw*360.0/(4.0*twfactor*M_PI*(FLOAT)Ptr));

  FLOAT J;
  printf("J[M] = %le\n", J = 4.0*M_PI*10000.0*wr*wtr*wtw/6.022);
  if (randomcombinations > 0) {
	printf("J_overlap_corrected[M] = %le\n", J/randomfraction);
  }
  printf("<U> (2 proteins) = %lf kT (%d samples)\n", Et/(FLOAT)N2, N2);
  printf("\n\n%lld combinations done. max bin equal %d\n", ncombinations, max);


  char sJ[50];
  sprintf(sJ, "Jfactor-%dbp_ID%d", nsteps, key); 
  FILE *fJ = openfwrite(sJ);

  fprintf(fJ, "J = %le M, total closed configurations = %d\n", J, Ptw);
  fprintf(fJ, "sample size = %d squared, %lld total configurations\n\n", N, (long long)N*(long long)N);

  fprintf(fJ, "Epsilon values: R < %lf Ang, cos(y) > %lf, tau < %lf\n", radi, gam, twi);
  fprintf(fJ, "Proteins:\n");
  for (int q = 0; q < ntypes; q++) {
    fprintf(fJ, "%s (%d steps), w=%lf\n", pnames[q], psteps[q], Pp[q]);
  }
  for (int q = 0; q < nbounds; q++) {
    fprintf(fJ, "Boundary %d:\n", q);
    writematrix(fJ, W0[q]);
  }
  fprintf(fJ, "\nJ by boundary = ");
  for (int q = 0; q < nbounds; q++) {
    fprintf(fJ, "%le ", J*(FLOAT)Ntw[q]/(FLOAT)Ptw);
  }
  fprintf(fJ, "\n");
  fclose(fJ);


  for (int i = 0; i < nbounds; i++) {
    char s[50];
    if (Ntw[i] > 0) {
      sprintf(s, "writhe_distribution-%dbp_ID%d_B%d", nsteps, key, i);
      hwrithe[i]->printhistogram_mm(s);
      sprintf(s, "twist_distribution-%dbp_ID%d_B%d", nsteps, key, i);
      htwist[i]->printhistogram_pm(s);
      sprintf(s, "link_distribution-%dbp_ID%d_B%d", nsteps, key, i);
      hlink[i]->printhistogram_pi(s);
      sprintf(s, "R_g_distribution-%dbp_ID%d_B%d", nsteps, key, i);
      hRg[i]->printhistogram_pm(s);

      sprintf(s, "roll_distribution-%dbp_ID%d_B%d", nsteps, key, i);
      hroll[i]->printhistogram_pm(s);

      char posf[60];
      sprintf(posf, "protein_distances-%dbp_ID%d_B%d", nsteps, key, i);
      posh[i]->printhistogram_pi(posf, 2);
      sprintf(posf, "protein_positions-%dbp_ID%d_B%d", nsteps, key, i);
      posh2[i]->printhistogram_pi(posf, 1);
      sprintf(posf, "protein_number-%dbp_ID%d_B%d", nsteps, key, i);
      nh[i]->printhistogram(posf);
    }
  }
  if (make_pics || make_pics_sc) fclose(df);


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
*/
  P = 0.0;

  short nPs = 0;

  init_positions(-30, nsteps/2+30);
  for (int i = 0; i < N; i++) {

    long int *locations = place_proteins(this, 0, nsteps/2, -30, nsteps/2+30, 0, nPs);

    for (int k = 0; k < nPs; k++) {
      alocations[locations[k]/(nstructs*2)]++;
    }

    delete [] locations;

  } 

  init_positions(nsteps/2-30, nsteps+30);
  for (int i = 0; i < N; i++) {
    long int *locations = place_proteins(this, nsteps/2, nsteps, nsteps/2-30, nsteps+30, 1, nPs);

    for (int k = 0; k < nPs; k++) {
	  alocations[locations[k]/(nstructs*2)]++;
    }

    delete [] locations;

  }

  long int Nt = 0;
  int Ns = 0;

  for (int i = 0; i < nsteps; i++) {
    printf("%d %d\n", i, alocations[i]);
    Nt += alocations[i];
    if (alocations[i] != 0) Ns++;
  }
  printf("Average spacing = %f\n", (FLOAT)N/((FLOAT)Nt/(FLOAT)Ns));

  init_positions(-30, nsteps+30);

  histogram *rdist = new histogram(lround(3.4*(FLOAT)nsteps/10.0+5.0), 10.0);
  histogram *rpdist = new histogram(lround(3.4*(FLOAT)nsteps/10.0+5.0), 10.0);

  matrix avgA = identity(4);
  avgA.setv(1,1, 0.0);
  avgA.setv(2,2, 0.0);
  avgA.setv(3,3, 0.0);
  avgA.setv(4,4, 0.0);
  int max = 0;
  int Nn = 0;
  FLOAT avgX = 0.0;
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

   locations = place_proteins(this, 0, nsteps, -30, nsteps+30, 1, nPs);

    if (nPs > max) max = nPs;

    NPS += nPs;

    temp = identity(4);
    temp2 = identity(4);

    if (!overlap()) { 
      Nn++;

      for (int j = 0; j < nsteps; j++) {
        temp = temp * calculateW(v[j]);
        P += temp(3,4)-temp2(3,4);
        temp2 = temp2 * calculateW(v[j]);
      }
      avgX += sqrt(temp(1,4)*temp(1,4)+temp(2,4)*temp(2,4)+temp(3,4)*temp(3,4));
      avgA = avgA + temp;
      if ((make_pics) && (i < 30)) {
	printf("Making %d(%d) of 30 pics, %d proteins\n", i+1, i, nPs);
	char s1[40], s2[40];
	sprintf(s1, "structure-%dbp-ID%d-%d.dat", nsteps, key, i);
	sprintf(s2, "proteins-%dbp-ID%d-%d.pdb", nsteps, key, i);
	illustrate_config(s1, s2, locations, nPs);
      }
    }
    rdist->add_data(sqrt(temp(1,4)*temp(1,4)+temp(2,4)*temp(2,4)+temp(3,4)*temp(3,4))); 
    rpdist->add_data((sqrt(temp(1,4)*temp(1,4)+temp(2,4)*temp(2,4)+temp(3,4)*temp(3,4))), nPs);
    delete [] locations;
  }

  printf("<x>(0) = %lf\n", avgX/(FLOAT)N);
  printf("Overlap factor = %lf\n (N_actual = %d)\n", (FLOAT)Nn/(FLOAT)N, Nn);

  avgA = (1.0/(FLOAT)Nn)*avgA;
  writematrix(stdout, avgA);

  printf("Average A1:N = \n");
  writematrix(stdout, avgA);

  matrix accum = identity(4);
  for (int i = 0; i < 10000; i++) accum = accum * avgA;
  printf("A1:N^10000 = \n");
  writematrix(stdout, accum);

  printf("Persistence length (Ang) = %f, Avg proteins = %f, Avg elastic energy = %f kT\n", P/(FLOAT)Nn, NPS/(FLOAT)N, energy/(FLOAT)N);

//  rdist->printhistogram_pm("radial-HU");
//  rpdist->printhistogram_pm("radial-proteins");


//  for (int i = 0; i < nsteps; i++) {
//    printf("N_HU[%d] = %d\n", i, alocations[i]);
//  } 

}




void closure::solveclosureproblem(int N1, int M1, int reduced) {

  matrix W;

  matrix Wstart(4,4);

  read_reference("ref_lac1", &Wstart);  

  key = ckey[0];

  printf("Initializing proteins\n");
  init_proteins();
  printf("Done\n");

  WA = new conf[N1];

  printf("Key = %d %d %d\n", ckey[0], ckey[1], ckey[2]);

  //  FLOAT *P = new FLOAT[3*N1*nsteps];  // positions of 2nd -> last base (1st at 0 0 0)
  
  N = N1;

  matrix temp;
  matrix Wall;
  matrix tpt;
  matrix endth;

  stress_free_state(1);

  int sum = 0;
 
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
  init_positions(-30, nsteps/2+30); 

  int olap = 0;

  for (int i = 0; i < N; i++) {

    do {
	WA[i].key0 = ckey[0];
	WA[i].key1 = ckey[1];
	WA[i].key2 = ckey[2];

    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

    }

    WA[i].locs = place_proteins(this, 0, nsteps/2, -30, nsteps/2+30, 0, WA[i].nP);
//    sortlocs(WA[i].locs, WA[i].nP);

    temp = identity(4);
    for (int j = 0; j < nsteps/2; j++) {
       temp = temp * calculateW(v[j]);
    }

    endth = calculatetp(Wstart*temp);

    W = endth;

    WA[i].W = W;
	sum++;
	if ((checkoverlap) && (olap = overlap1(WA[i].locs, WA[i].nP))) {
/*	  for (int f = nsteps/2; f < nsteps; f++) v[f] = I[bpstep[f]];
	  for (int f = 0; f < WA[i].nP; f++) printf("%d ", WA[i].locs[f]/(2*nstructs));
	  printf("\n");
	  if (WA[i].nP == 1) { illustrate_config("bad.dat", "pro.pdb", WA[i].locs, WA[i].nP);
	  return; } */
	}
    } while ((checkoverlap) && (olap));

    avgP += WA[i].nP;

  }

  avgP /= (FLOAT)N;

  printf("average protein (first half) = %f, %d halves discarded\n", avgP, sum-N);

  if (reduced) {

    reducedcombination();

  } else combination();

  printf("done.\n");

  delete [] WA;
//  delete [] HU;


}


