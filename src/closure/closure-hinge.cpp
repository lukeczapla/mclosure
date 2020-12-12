
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "closure-hinge.h"
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
  ckey = A;
  srand48(A);
}

closure::closure(long int key) {
  srand48(key);
  ckey = key;
  printf("srand48 key = %ld\n", key);
}

closure::~closure() {
  delete [] boundary;
  delete [] z;
  delete [] lambda;
  delete [] factor;
}


void closure::init_HU() {
  if (PHU == 0.0) return;
  HU.nsteps = 14;
  HU.v = new matrix[HU.nsteps];
  for (int i = 0; i < 14; i++) HU.v[i].setsize(6,1);
  HU.readconfig("HU.dat");
}

void closure::set_HU(FLOAT x) {
  PHU = x;
}

void closure::place_HU(int a, int b) {

  if (PHU == 0.0) return;

  int x = a+(int)((b-a)*drand48());

  int p = x;

  int q = 0;

  for (; x < b; x++) {
    if (drand48() < PHU) {
      if (q == 0) q = x;
      for (int j = 0; j < 16; j++) {
        v[j+x] = HU.v[j];
      }
      x += 15;
    }
  }

  if ((q == 0) || (q > p+16)) q = p+16;

  for (int i = a; i < q-16; i++) {
    if (drand48() < PHU) {
      for (int j = 0; j < 16; j++) {
        v[j+i] = HU.v[j];
      }
      i += 15;
    }
  }

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
    writematrix(stdout, W0[q]);
    Ntw[q] = 0;
  }

  matrix W1;

  matrix pretp(6,1);
  matrix posttp(6,1);

  pretp.setv(1,1,  4.8158);
  pretp.setv(2,1,  58.0379);
  pretp.setv(3,1,  64.6256);
  pretp.setv(4,1,  -17.8666);
  pretp.setv(5,1,  -2.2604);
  pretp.setv(6,1,  74.5757);




  posttp.setv(1,1,  -4.8233);
  posttp.setv(2,1,  58.0204);
  posttp.setv(3,1,  64.6245);
  posttp.setv(4,1,  17.8621);
  posttp.setv(5,1,  -2.2528);
  posttp.setv(6,1,  74.5819);


  matrix preW, postW;

  preW = calculateW(pretp);
  postW = calculateW(posttp);


  ncombinations = 0;

  int nZ = (int)(((FLOAT)(nsteps)*3.4+160.0)/epsilon);
  int nY = (int)(((FLOAT)(nsteps)*3.4+160.0)/epsilon);
  int nX = (int)(((FLOAT)(nsteps)*3.4+160.0)/epsilon);

  int i;

  int dis = (int)ceil(radi/epsilon);

  element_top *E = new element_top[nZ*nY*nX];
  for (i = 0; i < nZ*nY*nX; i++) {
    E[i].n = 0;
  }

  printf("step 2: binning configurations\n");

  for (i = 0; i < N; i++) {
    W1 = preW*calculateW(WA[i]);
    FLOAT z = -(W1(1,3)*W1(1,4)+W1(2,3)*W1(2,4)+W1(3,3)*W1(3,4));
    FLOAT y = -(W1(1,2)*W1(1,4)+W1(2,2)*W1(2,4)+W1(3,2)*W1(3,4));
    FLOAT x = -(W1(1,1)*W1(1,4)+W1(2,1)*W1(2,4)+W1(3,1)*W1(3,4));

    addelement(E[int(z/epsilon+nZ/2)+nZ*int(y/epsilon+nY/2)+nZ*nY*int(x/epsilon+nX/2)], i);
  }


  FLOAT rc;

  long long Pr = 0;

  matrix WW;
  matrix W2;
  matrix newW;
  matrix tpt;
  matrix tptp;
  matrix W;

  matrix WM;

  char ss[40];
  sprintf(ss, "hinge-coords-%d", nsteps);

  histogram *rad = new histogram(300, 1.0);
  trajectory *TT = new trajectory(ss, 1, WRITE_TRAJECTORY);
  bdna btemp(1);

  printf("\nstep 3: generating second half segments and combining\n");

  for (i = 0; i < N; i++) {

    generate_config2();

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
      v[ff] = tptp;
    }

    WW = identity(4);
    for (int ff = nsteps/2; ff < nsteps; ff++) {
      WW = WW * calculateW(v[ff]);
    }

    for (int q = 0; q < nbounds; q++) {

      WM = WW * postW * W0[q];

      FLOAT Z = WM(3,4);
      FLOAT Y = WM(2,4);
      FLOAT X = WM(1,4);
      int Zindex = (int)(Z/epsilon+nZ/2);
      int Yindex = (int)(Y/epsilon+nY/2);
      int Xindex = (int)(X/epsilon+nX/2);
   
      for (int m = (Zindex-dis < 0 ? 0: Zindex-dis); m <= (Zindex+dis > nZ-1 ? nZ-1 : Zindex+dis); m++)
      for (int n = (Yindex-dis < 0 ? 0: Yindex-dis); n <= (Yindex+dis > nY-1 ? nY-1 : Yindex+dis); n++)
      for (int o = (Xindex-dis < 0 ? 0: Xindex-dis); o <= (Xindex+dis > nX-1 ? nX-1 : Xindex+dis); o++) {
      if (30.0*30.0 < epsilon*epsilon*((Zindex-m)*(Zindex-m)+(Yindex-n)*(Yindex-n)+(Xindex-o)*(Xindex-o)) < 44.0*44.0)          
	for (int l = 0; l < E[m+nZ*n+nZ*nY*o].n; l++) {

            W = WA[getelement(E[m+nZ*n+nZ*nY*o], l)];

	    newW = preW*calculateW(W)*WM;
	    tpt = calculatetp(newW);
            rc = tpt(4,1)*tpt(4,1)+tpt(5,1)*tpt(5,1)+tpt(6,1)*tpt(6,1);
            if (((rc > 34.0*34.0) && (rc < 40.0*40.0)) && ((tpt(6,1) < -34.0) && (tpt(6,1) > -40.0))) 
		  if ((tpt(1,1) < 30.0) && (tpt(1,1) > -30.0)) if ((tpt(3,1) < 33.0) && (tpt(3,1) > -33.0)) {
	      Pr++;

		  matrix Wref = preW*calculateW(W)*WW;
//		  matrix Wi = identity(4);
//		  Wi.setv(2,2, -Wi(2,2));
//		  Wi.setv(3,3, -Wi(3,3));
		  
		  matrix postrev = Wref;//*postW;

	//	  postrev.setv(1,2, -postrev(1,2));
	//	  postrev.setv(2,2, -postrev(2,2));
	//	  postrev.setv(3,2, -postrev(3,2));
	//	  postrev.setv(1,3, -postrev(1,3));
	//	  postrev.setv(2,3, -postrev(2,3));
	//	  postrev.setv(3,3, -postrev(3,3));


		  char s[40];
		  sprintf(s, "lac_ref_%lld", Pr);
		  FILE *fref = fopen(s, "w");
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,4), postrev(2,4), postrev(3,4));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,1), postrev(2,1), postrev(3,1));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,2), postrev(2,2), postrev(3,2));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,3), postrev(2,3), postrev(3,3));

		  postrev = postrev * postW;
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,4), postrev(2,4), postrev(3,4));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,1), postrev(2,1), postrev(3,1));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,2), postrev(2,2), postrev(3,2));
		  fprintf(fref, " %lf %lf %lf \n", postrev(1,3), postrev(2,3), postrev(3,3));


		  fclose(fref);

		  

		  strcat(s, ".dat");
		  fref = fopen(s, "w");
		  writematrix(fref, newW);
		  writematrix(fref, calculatetp(newW));
		  fclose(fref);

	      btemp.v[0] = tpt;
	      TT->add_coordinates(&btemp);
		  FLOAT myR = 0.0;
	      matrix R = invert(preW)*newW*invert(postW);
	      matrix tpR = calculatetp(R);
          rad->add_data(myR=sqrt(tpR(4,1)*tpR(4,1)+tpR(5,1)*tpR(5,1)+tpR(6,1)*tpR(6,1)));
	    }

            ncombinations++;

          }
        }
    }
  }

  char s[40];
  sprintf(s, "radial-hinge-%d", nsteps);
  rad->printhistogram(s);

  delete rad;
  delete TT;

  FLOAT wr;

  printf("N_r     = %lld\n", Pr);
  printf("W(r=37) = %e\n", wr = (FLOAT)Pr/((FLOAT)N*(FLOAT)N*4.0*M_PI*(pow(40.0,3.0)-pow(34.0,3.0))/3.0));

  printf("\nJ[M] (W(0)/N_A) = %e\n", 10000*wr/6.022);

  printf("\n%lld combinations done. max bin equal %d\n", ncombinations, max);

 
  delete [] E;

  delete [] W0;

}


void closure::persistence(int N) {

  printf("N = %d, calculating persistence length, averages...\n", N);


  FLOAT energy = 0.0;
  FLOAT P = 0.0;

  FLOAT tw = 0.0;

  FLOAT Ptw = 0.0;

  matrix tpt;

  matrix temp, temp2; 

  for (int i = 0; i < N; i++) {

    generate_config1();
    generate_config2();

    for (int j = 0; j < nsteps; j++) {
      tpt = (z[bpstep[j]]*v[j]);
      tw += tpt(3,1);
      tpt = tpt + I[bpstep[j]];
      v[j] = tpt;
      Ptw += cos(M_PI*tw/180.0);
      energy += calculate_elenergy(v[j], j);
    }

    temp = identity(4);
    temp2 = identity(4);

    for (int j = 0; j < nsteps; j++) {
      temp = temp * calculateW(v[j]);
      P += temp(3,4)-temp2(3,4);
      temp2 = temp2 * calculateW(v[j]);
    }

  }

  printf("Persistence length (Ang) = %f, Twist persistence length (bp) = %f, Avg elastic energy = %f kT\n", P/(FLOAT)N, Ptw/(FLOAT)N, energy/(FLOAT)N);

}



void closure::solveclosureproblem(int N1, int M1, int reduced) {

  matrix W;


  WA = new matrix[N1];

  //  FLOAT *P = new FLOAT[3*N1*nsteps];  // positions of 2nd -> last base (1st at 0 0 0)
  
  N = N1;

  matrix temp;
  matrix tempW, Wall;
  matrix tpt;
  matrix endth;

  stress_free_state(1);


//  for (int i = 0; i < nbounds; i++) 
//    for (int j = 0; j <= 100; j++) p[i][j] = 0;

  printf("\nstep 1: generating %d half-configurations\n", N);

  FLOAT avgenergy = 0.0;
 


  for (int i = 0; i < N; i++) {

    FLOAT energy = 0.0;
    temp = identity(4);
    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

    }

    for (int j = 0; j < nsteps/2; j++) {
       energy += calculate_elenergy(v[j], j);
       tempW = calculateW(v[j]);
       temp = temp * tempW;
    }

    avgenergy += energy;

    endth = calculatetp(temp);

    W = endth;

    WA[i] = W;

  }

  avgenergy /= (FLOAT)N;

  printf("average energy ~= %e\n", 2.0*avgenergy);

  if (reduced) {

    reducedcombination();

  } 

  printf("done.\n");

}


