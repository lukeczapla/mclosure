
#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"

FLOAT randg();
matrix ringclosureboundary();

struct conf {
  short int nP;
  long int *locs;
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


  int make_pics;
  int make_pics_sc;

  long int ckey;

  bdna *proteins;

  conf *WA;


  histogram *position;
  histogram *nh;
  histogram *spacing;

  matrix *boundary;
  FLOAT epsilon;

  int N,M;

  FLOAT radi, gam, twi;

  matrix *factor;

  int ntypes;
  int nstructs;
  int *nproteins;
  int *structtype;
  char pfilenames[10][20];
  char pnames[10][20];
  int *psteps;
  int *pstartindex;
  FLOAT *Pp;
  int nbounds;

public:


  long long ncombinations;

  closure();
  closure(long int key);

  ~closure();

  inline void set_pics(int p) { make_pics = p; };
  inline void set_pics_sc(int p) { make_pics_sc = p; };

  void init_pairfile();
  int open_pairfile();
  void write_pair(int a, int b, int bound);
  void close_pairfile();
  int read_pair(int &a, int &b, int &bound);
  void swap_positions(short a, short a2); 
  void init_positions(short a, short a2);
  void regeneratesystem(int);

  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);

  inline void fast_translate(int a, int b) {
    for (int i = a; i < b; i++)
      v[i] = (z[bpstep[i]]*v[i])+I[bpstep[i]];
  };

//  inline void set_histograms(char *pH, char *sH, char *nH) {
//    strcpy(positionH, pH);
//    strcpy(spacingH, sH);
//    strcpy(nhuH, nH);
//  };


  void persistence(int N);

  void generate_config1();
  void generate_config2();

  void set_protein(int n, FLOAT x);
  void init_proteins();

  void place_protein(bdna *b, short int a, short int type, char reverse);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np, short exclude);

  FLOAT *Pproteins;

  void add_sort(long int *locs, int nlocs, long int X);
  void reducedcombination();
  
  void setepsilon(FLOAT e);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int N1, int M1, int reduced);

};
