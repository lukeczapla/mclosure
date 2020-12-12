
#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"

FLOAT randg();
matrix ringclosureboundary();

class closure : public bdna {

  matrix *z;
  matrix *lambda;

  bdna HU;

  matrix *WA;

  matrix *boundary;
  FLOAT epsilon;

  int N,M;

  FLOAT radi, gam, twi;

  matrix *factor;

  long int ckey;

  int nbounds;
  FLOAT PHU;


public:


  long long ncombinations;

  closure();
  closure(long int key);

  ~closure();

  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);

  void generate_config1();
  void generate_config2();

  void init_HU();
  void place_HU(int a, int b);

  void set_HU(FLOAT x);

  void persistence(int n);

  void reducedcombination();

  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int N1, int M1, int reduced);

};
