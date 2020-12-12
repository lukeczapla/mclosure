
#include <pthread.h>

#include "../bdna/bdna.h"
#include "../spatial/spatial.h"
#include "../histogram/histogram.h"
#include "../proteinPDB/proteinPDB.h"

matrix ringclosureboundary();

struct conf {
  float W[6];
  short int nP;
  long int *locs;
  unsigned short key0;
  unsigned short key1;
  unsigned short key2;
};

struct element_top {
  int n;
  int *list;
};


class closure : public bdna {


  element_top *E;
  conf *WA;

  pthread_mutex_t mutex;

  matrix *z;
  matrix *lambda;

  int make_pics;
  int make_pics_sc;

  proteinPDB *proteinPDBs;
  bdna *proteins;

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
  char pfilenames[10][20];
  char pnames[10][20];
  int *psteps;
  int *pstartindex;
  int nnodes;
  FLOAT *Pp;

  int *Nrgt;

  int nbounds;
  int randomcombinations;

  unsigned short skey;
  unsigned short ckey[3];


public:


  long long ncombinations;

  closure();
  closure(unsigned short *key);

  ~closure();

  void addelement(struct element_top &e, int index);

  inline void set_pics(int p) { make_pics = p; };
  inline void set_pics_sc(int p) { make_pics_sc = p; };

  void setnodes(int nn);

  void pushkey(unsigned short *k1, unsigned short *k2);
  void popkey(unsigned short *k1, unsigned short *k2);
  void pushq();
  void popq();

  void swap_positions(short a, short a2); 
  void init_positions(short a, short a2);

  void m_swap_positions(short *q1, short a, short a2, unsigned short *key, struct drand48_data *buffer);
  void m_init_positions(short *q1, short a, short a2, unsigned short *key, struct drand48_data *buffer);

  void initialize_params();

  inline void cleanup() { pthread_mutex_destroy( &mutex ); };

  void setparams(FLOAT a, FLOAT b, FLOAT c);
  void setoverlap(int o);

  inline void removeoverlap(int ng) { neglect = ng; };

  inline void fast_translate(int a, int b) {
    for (int i = a; i < b; i++)
      v[i] = (z[bpstep[i]]*v[i])+I[bpstep[i]];
  };

  inline unsigned short MCkey() {
    return skey;
  }

//  inline void set_histograms(char *pH, char *sH, char *nH) {
//    strcpy(positionH, pH);
//    strcpy(spacingH, sH);
//    strcpy(nhuH, nH);
//  };


  void m_generate_config1(bdna *b, unsigned short *key, struct drand48_data *buffer);
  void m_generate_config2(bdna *b, unsigned short *key, struct drand48_data *buffer);


  void convert_steps();

  void set_protein(int n, FLOAT x);
  void init_proteins();


  long int *m_place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np, unsigned short *key, struct drand48_data *buffer, short *q1);
  long int *m_place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np, short exclude, unsigned short *key, struct drand48_data *buffer, short *q1);


  void place_protein(bdna *b, short int a, short int type, char reverse);
  void place_proteins(bdna *b, long int *locs, int nP);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np);
  long int *place_proteins(bdna *b, short a, short a2, short b1, short b2, char overhang, short &np, short exclude);

  void illustrate_config(bdna *b, char *f1, char *f2, long int *locs, int nP);

  int overlap1(long int *locs, int nP);
  int overlap2(long int *locs, int nP);

  int overlap();
  int overlap(long int *locs, int nP);
  int overlap(bdna *b, long int *locs, int nP);

  FLOAT *Pproteins;

  void writeJ(FILE *out);
  void initE(int N1);

  void add_sort(long int *locs, int nlocs, long int X);
  void reducedcombination(int node);
  
  void setepsilon(FLOAT e);
  void setrandomsamples(int o);
  void setboundarynum(int n);
  void addboundary(matrix bc);
  void solveclosureproblem(int node, int offset);

};
