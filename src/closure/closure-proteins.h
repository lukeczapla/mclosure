
#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../proteinPDB/proteinPDB.h"

#define MAX_PROTEINS 500
#define DNA_DNA_DIST 20.0

FLOAT randg();
matrix ringclosureboundary();

struct conf {
  short int nP;
  long int *locs;
  matrix W;
  unsigned short key0;
  unsigned short key1;
  unsigned short key2;
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


  proteinPDB *proteinPDBs;
  bdna *proteins;

  conf *WA;


  histogram *position;
  histogram *nh;
  histogram *spacing;

  matrix *boundary;
  FLOAT epsilon;

  int N,M;

  FLOAT radi, gam, twi;
  int checkoverlap;

  matrix *factor;

  int neglect;

  int ntypes;
  int nstructs;
  int *nproteins;
  int *structtype;
  char pfilenames[MAX_PROTEINS][20];
  char pnames[MAX_PROTEINS][20];
  int *psteps;
  int *pstartindex;
  FLOAT *Pp;
  int nbounds;
  int randomcombinations;

  unsigned short key;
  unsigned short tempkey[3];
  unsigned short ckey[3];


public:


  long long ncombinations;

  closure();
  closure(unsigned short *key);

  ~closure();

  inline void printparams() {
    for (int i = 0; i < nsteps; i++) writematrix(stdout, factor[i]);
  }


  inline void set_pics(int p) { make_pics = p; };
  inline void set_pics_sc(int p) { make_pics_sc = p; };

  void printeigenvalues();

  void pushkey();
  void popkey();
  void pushq();
  void popq();

  void swap_positions(short a, short a2); 
  void init_positions(short a, short a2);

  void initialize_params();

  void setparams(FLOAT a, FLOAT b, FLOAT c);
  void setoverlap(int o);

  inline void removeoverlap(int ng) { neglect = ng; };

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
  void convert_steps();

  void set_protein(int n, FLOAT x);
  void init_proteins();


  void place_protein(bdna *b, short int a, short int type, char reverse);
  void place_proteins(bdna *b, long int *locs, int nP);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np, short exclude);

  void illustrate_config(char *f1, char *f2, long int *locs, int nP);

  int overlap1(long int *locs, int nP);
  int overlap2(long int *locs, int nP);

  int overlap();
  int overlap(long int *locs, int nP);

  FLOAT *Pproteins;

  void add_sort(long int *locs, int nlocs, long int X);
  void combination();
  void reducedcombination();
  
  void setepsilon(FLOAT e);
  void setrandomsamples(int o);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int N1, int M1, int reduced);

};
