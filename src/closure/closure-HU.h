
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

struct pair {
  int a;
  int b;
  int bound;
};


class closure : public bdna {

  matrix *z;
  matrix *lambda;


  long int ckey;

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
  FLOAT PHU;

public:


  long long ncombinations;

  closure();
  closure(long int key);

  ~closure();

  void init_pairfile();
  int open_pairfile();
  void write_pair(int a, int b, int bound);
  void close_pairfile();
  int read_pair(int &a, int &b, int &bound);

  void regeneratesystem(int);

  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);

  inline void fast_translate(int a, int b) {
    for (int i = a; i < b; i++)
      v[i] = (z[bpstep[i]]*v[i])+I[bpstep[i]];
  };

  inline void set_histograms(char *pH, char *sH, char *nH) {
    strcpy(positionH, pH);
    strcpy(spacingH, sH);
    strcpy(nhuH, nH);
  };

  inline bdna get_HU(int i) {return HU[i];};

  void persistence(int N);

  void generate_config1();
  void generate_config2();

  void set_HU(FLOAT x);
  void init_HU();

  void place_HU(bdna *b, int a);
  short *place_HU(bdna *b, short a, short a2, short pref, short &nhu);
  short *place_HU_forward(bdna *b, short a, short a2, short pref, short &nhu);
  short *place_HU_reverse(bdna *b, short a, short a2, short pref, short &nhu);

  void reducedcombination();
  
  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int N1, int M1, int reduced);

};
