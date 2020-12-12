#ifndef _CLOSURE_RAD_H_
#define _CLOSURE_RAD_H_

#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../proteinPDB/proteinPDB.h"

FLOAT randg();

class closure : public bdna {

  matrix *z;
  matrix *lambda;

  matrix *WA;


  int N,M;
  FLOAT bin_width;

  matrix *factor;

  int nproteins;
  char pnames[10][20];
  bdna *proteins;
  int *pposition;

  proteinPDB *proteinPDBs;

  unsigned short tempkey[3];
  unsigned short ckey[3];

public:

  int skey;

  long long ncombinations;

  closure();
  closure(unsigned short key[3]);

  ~closure();

  void printeigenvalues();
  void initialize_params();
  void setbinsize(double sz);

  void read_constraints();
  void place_constraints();

  void generate_config1();
  void generate_config2();

  int overlap1();
  int overlap2();
  int overlap();

  void combination(int N, char *rfile);

};

#endif

