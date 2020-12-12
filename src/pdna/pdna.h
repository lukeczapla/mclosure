#ifndef _PDNA_H
#define _PDNA_H

#include "../bdna/bdna.h"
#include "../proteinPDB/proteinPDB.h"

FLOAT randg();

struct protein {
	int pos;
	int type;
	int str;
	int dir;
};

class pdna : public bdna {

  matrix *z;
  matrix *lambda;

  proteinPDB *proteinPDBs;
  bdna *proteins;

  matrix *factor;

  int ntypes;
  int nstructs;
  int *nproteins;
  int *structtype;
  char pfilenames[100][20];
  char pnames[100][20];
  int *psteps;
  int *pstartindex;
  FLOAT *Pp;
FLOAT *Pp2;
  matrix *vchange;

  matrix *v2;
  matrix *vchange2;
bdna *q;
  protein *pos;
  protein *poschange;
  int nP, nPchange;

public:

  	pdna();
  	~pdna();


  	void printeigenvalues();

  	void initialize_params(double b);

  	void generate_config();
  	void convert_steps();

  	void set_protein(int n, FLOAT x);
  	void init_proteins();

  	double move(double beta);

	void printnew3dna(char *filename);
  	void place_protein(bdna *b, short int a, short int type, char reverse);

  	void illustrate_config(char *f1, char *f2);

  	int overlap();
  	int overlap(long int *locs, int nP);

	void revert();
	void accept();
	int isvalid();
	double omega(int N, int N2, int m);

	inline int Pnum() { return nP; }
	inline int Pnumnew() { return nPchange; }
	void printP();
	
	int oldEEp();
	int newEEp();
        int oldEEbin();
        int newEEbin();
	double oldsigma();
	double newsigma();
	double oldEEfull();
	double newEEfull();

	
	int oldEEt();
	int newEEt();
	double chainS(int x1, int x2);
	double twist();
	int newEE();
	int oldEE();
	
	double dE(double &elOld, double &elNew);

  	FLOAT *Pproteins;

};

#endif

