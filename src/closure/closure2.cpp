
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../trajectory/trajectory.h"
#include "closure2.h"
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
  HU.nsteps = 16;
  HU.v = new matrix[HU.nsteps];
  for (int i = 0; i < 16; i++) HU.v[i].setsize(6,1);
  HU.readconfig("HU.dat");
}

void closure::set_HU(FLOAT x) {
  PHU = x;
}

void closure::place_HU(int a) {
  for (int i = a; i < a+16; i++)
    v[i] = HU.v[i-a];
}

void closure::place_HU(int a, int b) {

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
    for (j = 1; j <= 6; j++) 
      lambda[i].setv(j, 1, d[j-1]);

//    writematrix(stdout, lambda[i]);

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



void closure::reducedcombination(int hu) {

  int dis = (int)ceil(radi/epsilon);

  matrix *W0 = new matrix[nbounds];
  for (int q = 0; q < nbounds; q++) W0[q] = invert(calculateW(boundary[q]));

  matrix W1;

  int nZ = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nY = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;
  int nX = (int)(((FLOAT)(nsteps)*3.4)/epsilon)+2;  

  int i;

  element_top *E = new element_top[nZ*nY*nX];
  for (i = 0; i < nZ*nY*nX; i++) {
    E[i].n = 0;
  }

  printf("step 2: binning configurations\n");

  for (i = 0; i < N; i++) {
    W1 = calculateW(WA[i]);
    FLOAT z = -(W1(1,3)*W1(1,4)+W1(2,3)*W1(2,4)+W1(3,3)*W1(3,4));
    FLOAT y = -(W1(1,2)*W1(1,4)+W1(2,2)*W1(2,4)+W1(3,2)*W1(3,4));
    FLOAT x = -(W1(1,1)*W1(1,4)+W1(2,1)*W1(2,4)+W1(3,1)*W1(3,4));

    addelement(E[int(z/epsilon+nZ/2)+nZ*int(y/epsilon+nY/2)+nZ*nY*int(x/epsilon+nX/2)], i);
  }

  FLOAT rc;

  FLOAT factor = radi*radi; //now 10 ang, was 4% cl or FLOAT factor = 0.0016*cl*cl;
  FLOAT trfactor = gam;
  FLOAT twfactor = twi;

  matrix WW;
  matrix WM;
  matrix W2;
  matrix tptp;
  matrix newW;
  matrix tpt;

  printf("step 3: generating second half segments and combining\n");
  for (i = 0; i < N; i++) {
    generate_config2();
    WW = identity(4);
    for (int ff = nsteps/2; ff < nsteps; ff++) {
      tptp = (z[bpstep[ff]]*v[ff])+I[bpstep[ff]];
      v[ff] = tptp;
    }

    place_HU(nsteps/2-hu-1);
    place_HU(nsteps/2-hu+15, nsteps-17);
    

    for (int ff = nsteps/2; ff < nsteps; ff++) {
      WW = WW * calculateW(v[ff]);
    }

    for (int q = 0; q < nbounds; q++) {

    WM = WW * W0[q];

    FLOAT Z = WM(3,4);
    FLOAT Y = WM(2,4);
    FLOAT X = WM(1,4);
    int Zindex = (int)(Z/epsilon+nZ/2);
    int Yindex = (int)(Y/epsilon+nY/2);
    int Xindex = (int)(X/epsilon+nX/2);
    
    matrix W;

    for (int m = (Zindex-dis < 0 ? 0: Zindex-dis); m <= (Zindex+dis > nZ-1 ? nZ-1 : Zindex+dis); m++)
    for (int n = (Yindex-dis < 0 ? 0: Yindex-dis); n <= (Yindex+dis > nY-1 ? nY-1 : Yindex+dis); n++)
    for (int o = (Xindex-dis < 0 ? 0: Xindex-dis); o <= (Xindex+dis > nX-1 ? nX-1 : Xindex+dis); o++) {
      if (epsilon*epsilon*((Zindex-m)*(Zindex-m)+(Yindex-n)*(Yindex-n)+(Xindex-o)*(Xindex-o)) < radi*radi+2.0*epsilon*epsilon)
        for (int l = 0; l < E[m+nZ*n+nZ*nY*o].n; l++) {

          W = WA[getelement(E[m+nZ*n+nZ*nY*o], l)];

	  newW = calculateW(W)*WM;
	  tpt = calculatetp(newW);
          rc = tpt(4,1)*tpt(4,1)+tpt(5,1)*tpt(5,1)+tpt(6,1)*tpt(6,1);
          if (rc < factor) {
	    Pr++;
	    if (newW(3,3) > trfactor) {
	      Ptr++;
	      if (fabsf(tpt(3,1)) < twfactor) {
		Ntw[q]++;
	      }
	    };
	  }

          ncombinations++;

        }
    }
    }
  }

  if (hu == 15) {

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

  printf("N_r   = %lld %e\n", Pr, wr = (FLOAT)Pr/(16.0*(FLOAT)N*(FLOAT)N*4.0*M_PI*sqrt(factor)*sqrt(factor)*sqrt(factor)/3.0));
  printf("N_gam = %d %e\n", Ptr, wtr = (FLOAT)Ptr/((1.0-trfactor)*(FLOAT)Pr));
  printf("N_tw  = %d %e\n", Ptw, wtw = (FLOAT)Ptw*360.0/(4.0*twfactor*M_PI*(FLOAT)Ptr));

  printf("J[M] = %le\n", 4.0*M_PI*10000*wr*wtr*wtw/6.022);

  printf("\n\n%lld combinations done. max bin equal %d\n", ncombinations, max);

  }

  for (i = 0; i < nZ*nY*nX; i++) {
    if (E[i].n > 0) delete [] E[i].list;
  }
                      
  delete [] E;

}


void closure::solveclosureproblem(int N1, int M1, int reduced) {

  matrix W;

  init_HU();

  ncombinations = 0;

  N = N1;

  WA = new matrix[N];

  Ntw = new int[nbounds];
  for (int q = 0; q < nbounds; q++) Ntw[q] = 0;

  matrix temp;
  matrix tempW, Wall;
  matrix tpt;
  matrix endth;

  stress_free_state(1);

  for (int hu = 0; hu < 16; hu++) {

  printf("\ncalculation %d of 16\n", hu+1);
  printf("step 1: generating %d half configurations\n", N);

  FLOAT avgenergy = 0.0;
 
  for (int i = 0; i < N; i++) {

    FLOAT energy = 0.0;
    temp = identity(4);
    generate_config1();

    for (int j = 0; j < nsteps/2; j++) {

       tpt = (z[bpstep[j]]*v[j])+I[bpstep[j]];
       v[j] = tpt;

    }

    place_HU(nsteps/2-hu-1);

    place_HU(1, nsteps/2-hu-17);

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

  if (reduced) {

    reducedcombination(hu);

  }

  }

  delete [] WA;

  printf("done.\n");

}


