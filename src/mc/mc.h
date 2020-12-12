#ifndef _MC_H
#define _MC_H

#include "../pdna/pdna.h"
#include "../parameters/parameters.h"
#include "../socket/socket.h"

double mold(pdna *b, parameter_type *p);
double mnew(pdna *b, parameter_type *p);

void runsimulation(pdna *b, parameter_type *p, socket_type *c, double beta, char *host, int node, int Neq, char *savefile, int Nc);

#endif

