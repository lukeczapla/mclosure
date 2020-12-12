#ifndef TRANSLATE_TOOLS
#define TRANSLATE_TOOLS

#include "matrix/matrix.h"
#include "bdna/bdna.h"
#include "defines.h"


void read_reference(char *fname, matrix *W);
void rewrite_pdb(char *in, char *out, matrix W, matrix new1, char append);

void translate_position(matrix W, float &x, float &y, float &z);

FLOAT det3(matrix A);

#endif
