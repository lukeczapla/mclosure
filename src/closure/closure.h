#ifndef _CLOSURE_H_
#define _CLOSURE_H_

#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../proteinPDB/proteinPDB.h"

#define MAX_PROTEINS 500
#define DNA_DNA_DIST 20.0
#define NEAREST_CHECK 50

FLOAT randg();
matrix ringclosureboundary();

struct element_item {
  unsigned short list0;
  unsigned short list1;
  unsigned short list2;
  int index;
};

struct element_top_himem {
  int n;
  element_item *data;
};

struct element_top {
  int n;
  unsigned short *list0;
  unsigned short *list1;
  unsigned short *list2;
};

struct pair {
  int a;
  int b;
  int bound;
};

class closure : public bdna {

  matrix *z;
  matrix *lambda;

  int pairfile;

  matrix *WA;

  matrix *boundary;
  FLOAT epsilon;

  int himem;
  int N,M;
  int neglect;


  FLOAT radi, gam, twi;
  FLOAT bin_width;

  matrix *factor;
  int make_pics;
  int make_pics_sc;
  FLOAT sc_factor;

  int nbounds;

  int nproteins;
  char pnames[MAX_PROTEINS][20];
  bdna *proteins;
  int *pposition;

  proteinPDB *proteinPDBs;

  int restrictends;
  int checkoverlap;
  int randomcombinations;

  element_item *corr;

  unsigned short tempkey[3];
  unsigned short ckey[3];

public:

  int skey;

  long long ncombinations;

  closure();
  closure(unsigned short key[3]);

  ~closure();

  int overlap1();
  int overlap2();

  int overlap();

  void pushkey();
  void popkey();

  void setrestrictions(int r);

  void printeigenvalues();

  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);
  void setpics(int pics, int sc_pics, FLOAT sc_width);
  void setmem(int himemory);
  void setoverlap(int o);
  void setrandomsamples(int o);
  void setbinsize(FLOAT sz);

  inline void removeoverlap(int ng) { neglect = ng; };


  void read_constraints();
  void place_constraints();

  void printproteins(char *s);

  void generate_config1();
  void generate_config2();

  void persistence(int n);

  void reducedcombination(element_top *E, element_top_himem *EM);

  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);

  void solveextension(int N1, FLOAT fraction);
  void solveclosureproblem(int N1, int M1, int reduced);

};

#endif

