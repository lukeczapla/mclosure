#ifndef _PROTEINPDB_H
#define _PROTEINPDB_H

#include "../defines.h"
#include "../matrix/matrix.h"
#include "../bdna/bdna.h"

class proteinPDB {
protected:
  int natoms;
  FLOAT *x; FLOAT *y; FLOAT *z;
  FLOAT *r;

public:

  proteinPDB();
  proteinPDB(char *fname);
  proteinPDB(const proteinPDB &p);

  ~proteinPDB();

  int readPDB(char *fname);

  friend int proteinoverlap(proteinPDB p, proteinPDB q, matrix Wp, matrix Wq);
  friend int pDNAoverlap(proteinPDB p, matrix Wp, int pos, int plen, bdna *b, FLOAT *x, FLOAT *y, FLOAT *z);
  friend int pDNAoverlap(proteinPDB p, matrix Wp, int pos, int plen, FLOAT *x, FLOAT *y, FLOAT *z, int x1, int x2);

  friend proteinPDB *translateprotein(proteinPDB *p, matrix W);
  
};


#endif
