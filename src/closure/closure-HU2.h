
#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"

FLOAT randg();
matrix ringclosureboundary();

struct conf {
  short nhu;
  short *locs;
  matrix W;
};

class closure : public bdna {

  matrix *z;
  matrix *lambda;

  bdna *HU;

  conf *WA;

  char positionH[30];
  char nhuH[30];
  char spacingH[30];

  histogram *position;
  histogram *nh;
  histogram *spacing;

  matrix *boundary;
  FLOAT epsilon;

  int N,M;

  FLOAT radi, gam, twi;

  matrix *factor;

  int nbounds;
  int NHU;

public:


  long long ncombinations;

  closure();
  closure(long int key);

  ~closure();


  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);

  inline void set_histograms(char *pH, char *sH, char *nH) {
    strcpy(positionH, pH);
    strcpy(spacingH, sH);
    strcpy(nhuH, nH);
  };


  void persistence(int N);

  void generate_config1();
  void generate_config2();

  void set_HU(int x);
  void init_HU();

  void place_HU(bdna *b, int a);
  short *place_HU(bdna *b, short a, short a2, short pref, short &nhu);
  short *place_HU_half(bdna *b, short a, short a2, short pref, short &nhu, int nT);
  void reducedcombination();
  
  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int N1, int M1, int reduced);

};
