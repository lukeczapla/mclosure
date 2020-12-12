#ifndef _CLOSURE_H_
#define _CLOSURE_H_

#include <pthread.h>
#include <stdlib.h>

#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../proteinPDB/proteinPDB.h"

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

class closure : public bdna {

  matrix *z;
  matrix *lambda;

  element_top *E;
  element_top_himem *EM;
  float *Wmem;

  pthread_mutex_t mutex;

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
  int nnodes; 

  int nproteins;
  char pnames[10][20];
  bdna *proteins;
  int *pposition;

  proteinPDB *proteinPDBs;

  int restrictends;
  int checkoverlap;
  int randomcombinations;

  element_item *corr;

  int *Nrgt;

  unsigned short ckey[3];

public:

  closure();
  closure(unsigned short key[3]);

  ~closure();

  int overlap1();
  int overlap2();

  int overlap();
  int overlap(bdna *b);

  inline void cleanup() { pthread_mutex_destroy(&mutex); };

  void pushkey(unsigned short *tempkey, unsigned short *key);
  void popkey(unsigned short *tempkey, unsigned short *key);

  void setrestrictions(int r);


  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);
  void setpics(int pics, int sc_pics, FLOAT sc_width);
  void setmem(int himemory);
  void setoverlap(int o);
  void setrandomsamples(int o);
  void setbinsize(FLOAT sz);

  inline void removeoverlap(int ng) { neglect = ng; };

  matrix memW(int index);
  void setmemW(matrix tpt, int index);
  element_item *getelement(struct element_top_himem e, int n);
  unsigned short *getelement(struct element_top e, int n);
  void addelement(struct element_top &e, unsigned short *index);
  void addelement(struct element_top_himem &e, unsigned short *keyi, int index);

  void read_constraints();

  void printproteins(char *s);

  void m_generate_config1(bdna *b, unsigned short *key, drand48_data *buffer);
  void m_generate_config2(bdna *b, unsigned short *key, drand48_data *buffer);
  void m_place_constraints(bdna *b);

  void generate_config1(unsigned short *key);
  void generate_config2(unsigned short *key);
  void place_constraints();

  void persistence(int n);

  void setnodes(int n);

  void reducedcombination(int bkey);

  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);

  void writeJ(FILE *out);
  void initE(int N1);

  void solveclosureproblem(int node, int offset);

};

#endif

