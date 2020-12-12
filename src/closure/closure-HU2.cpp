
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "closure-HU2.h"
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


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
  printf("srand48 key = %ld\n", A);
  srand48(A);
}

closure::closure(long int key) {
  srand48(key);
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
    HU[k].nsteps = 16;
    HU[k].v = new matrix[HU[k].nsteps];
    for (int i = 0; i < 16; i++) HU[k].v[i].setsize(6,1);
    sprintf(s, "HU-%d.dat", k+1);
    HU[k].readconfig(s);
  }
}

void closure::set_HU(int x) {
  NHU = x;
}

void closure::place_HU(bdna *b, int a) {
  int c = (int)(4*drand48());
  for (int i = 0; i < 16; i++) 
    b->v[i+a] = HU[c].v[i];
}

short *closure::place_HU_half(bdna *b, short a, short a2, short pref, short &nhu, int nT) {

}

short *closure::place_HU(bdna *b, short a, short a2, short pref, short &nhu) {

  short s[1000];
  short *locations;
  nhu = 0;

  for (short x = a; x < a2; x++) {
    if ((drand48() < PHU) && ((pref == 0) || (x-pref < 20) || (x-pref > 22)) && ((nhu == 0) || (x-s[nhu-1] < 20) || (x-s[nhu-1] > 22))) {
      s[nhu++] = x;
      int c = (int)(4*drand48());
      for (int j = 0; j < 16; j++) {
        b->v[j+x] = HU[c].v[j];
      }
      x += 15;
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

  if (strlen(nhuH) > 0) nh = new histogram(20, 1.0);
  if (strlen(positionH) > 0) position = new histogram(nsteps, 1.0);
  if (strlen(spacingH) > 0) spacing = new histogram(nsteps, 1.0);


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
	b->v[Q] = v[Q];
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
	      pos = 22-(WA[elementID].locs[WA[elementID].nhu-1]-(nsteps/2-22));

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
                  short *lls = new short[WA[elementID].nhu+nhus[Q]];
		
	          for (int i = 0 ; i < WA[elementID].nhu; i++) lls[i] = WA[elementID].locs[i];
                  for (int i = 0 ; i < nhus[Q]; i++) lls[i+WA[elementID].nhu] = locs[Q][i];
	          for (int i = 0; i < nhus[Q]+WA[elementID].nhu; i++) {
		    if (strlen(positionH) > 0) position->add_data((FLOAT)(lls[i]+1));
		    if ((strlen(spacingH) > 0) && (i != 0)) spacing->add_data((FLOAT)(lls[i]-lls[i-1]));
		  }
		  if (strlen(nhuH) > 0) nh->add_data((FLOAT)(nhus[Q]+WA[elementID].nhu));
		 
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

  if (strlen(nhuH) > 0) nh->printhistogram(nhuH);
  if (strlen(positionH) > 0) position->printhistogram(positionH);
  if (strlen(spacingH) > 0) spacing->printhistogram(spacingH);

  delete nh;
  delete position;
  delete spacing; 
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

    locations2 = place_HU(this, ((nhu != 0) && (locations[nhu-1] > nsteps/2-16) ? locations[nhu-1]+16 : nsteps/2), nsteps-16, (nhu != 0 ? locations[nhu-1] : 0), nhu2);

    for (int k = 0; k < nhu2; k++) {
      alocations[locations2[k]]++;
    }

    Nhu += nhu+nhu2;

    temp = identity(4);
    temp2 = identity(4);

    for (int j = 0; j < nsteps; j++) {
      temp = temp * calculateW(v[j]);
      P += temp(3,4)-temp2(3,4);
      temp2 = temp2 * calculateW(v[j]);
    }

    delete [] locations;
    delete [] locations2;

  }

  printf("Persistence length (Ang) = %f, Avg HU = %f, Avg elastic energy = %f kT\n", P/(FLOAT)N, Nhu/(FLOAT)N, energy/(FLOAT)N);

  for (int i = 0; i < nsteps; i++) {
    printf("N_HU[%d] = %d\n", i, alocations[i]);
  } 

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


