#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mc.h"


double dUstar;

void clear(data_type *d) {
	for (int i = 0; i < d->size; i++) d->Ustar[i] = 0;
}



int bin(parameter_type *p, double v) {
	if ((v >= p->start) && (v < p->end)) return (int)((v-p->start)*(double)(p->Nbins)/(p->end-p->start));
	else if (v <= p->start) return 0;
	else return (p->Nbins - 1);
}


double penalty(parameter_type *p, data_type *d, double v) {
	if ((v >= p->start) && (v < p->end)) return (dUstar * d->Ustar[bin(p, v)]);
// if outside the bins, return "infinity"
	return 1e20;
}


int accept(double dE) {

	// if the new energy is lower, accept
	if (dE < 0.0) return 1;
	
	// if the random number between 0 and 1 is less than exp(-dE), where dE is positive, accept
	if (drand48() < exp(-dE)) return 1;
	
	// otherwise, reject
	return 0;
	
}

double mold(pdna *b, parameter_type *p) {
	if (p->measure == 1) return b->oldEE();
	if (p->measure == 2) return b->oldEEt();
	if (p->measure == 3) return b->Pnum();
	if (p->measure == 5) return b->oldEEfull();
	if (p->measure == 6) return b->oldsigma();
        if (p->measure == 7) return b->oldEEbin();
	if (p->measure == 8) return b->calculate_twist_open() + b->calculate_writhe_open() - 34.28*(float)(b->nsteps)/360.0;
}

double mnew(pdna *b, parameter_type *p) {
	if (p->measure == 1) return b->newEE();
	if (p->measure == 2) return b->newEEt();
	if (p->measure == 3) return b->Pnumnew();
	if (p->measure == 5) return b->newEEfull();
	if (p->measure == 6) return b->newsigma();
        if (p->measure == 7) return b->newEEbin();
        if (p->measure == 8) return b->calculate_twist_open() + b->calculate_writhe_open() - 34.28*(float)(b->nsteps)/360.0;
}


void runsimulation(pdna *b, parameter_type *p, socket_type *c, double beta, char *host, int node, int Neq, char *savefile, int Nc) {
	data_type *d = new data_type;
	data_type *upd;
	
// request the parameters that start the simulation
	d->size = p->Nbins;
	d->dUstar = p->dUstar;
	clear(d);
	upd = requestupdateparameters(c, d);
	d->key = upd->key;
	double measureold, measurenew;
	double Uold, Unew;

	measureold = mold(b, p);
	measurenew = measureold;
	double elOld, elNew;
	
	printf("starting measure = %lf, minimum = %d, equilibrating (%d moves)...\n", measureold, (int)p->start, Neq);
	
	double dE;
	int m;
	
	double S = 0.0;
	
	Uold = penalty(p, upd, measureold);
	
	for (int i = 0; i < Neq; i++) {
		S = b->move(beta);
		dE = b->dE(elOld, elNew);
		measurenew = mnew(b, p);
		Unew = penalty(p, upd, measurenew);
		if (accept(beta*(dE + S + Unew - Uold)) && (((measurenew < p->end) && (measurenew > p->start)) || (p->measure == 3))) {
			b->accept();
			measureold = measurenew;
			Uold = Unew;
		} else b->revert();
	}
	
	measureold = mold(b, p);
	
	printf("starting... measure = %lf\n", measureold);
	
	int i = 0;
	int j = 0;
	
	long int Nt = 0;
	long int Naccept = 0;
	
	long int avgp = 0;
	long int *M = new long int[2050];
	long int *Mt = new long int[2050];
	
	long int Mp[2050][100];
	for (int i = 0; i < 2050; i++) for (int j = 0; j < 100; j++) Mp[i][j] = 0;
	double avgE[2050];
	for (int i = 0; i < 2050; i++) {
		M[i] = 0;
		Mt[i] = 0;
		avgE[i] = 0.0;
	}
	
	int pmax = 0;
	S = 0.0;
	while (1) {
		S = b->move(beta);
		dE = b->dE(elOld, elNew);
		Nt++;
		measurenew = mnew(b, p);
		Unew = penalty(p, upd, measurenew);
		if (accept(beta*(dE + S + Unew - Uold))) {
			b->accept();
			Naccept++;
			measureold = measurenew;
			upd->Ustar[bin(p, measurenew)]++;
			d->Ustar[bin(p, measurenew)]++;
			Uold = Unew + dUstar;
			avgp += b->Pnum();
			M[bin(p, measurenew)] += b->Pnum();
			Mp[bin(p, measurenew)][b->Pnum()]++;
			Mt[bin(p, measurenew)]++;
			avgE[bin(p, measurenew)] += elNew;
			if (b->Pnum() > pmax) pmax = b->Pnum();
		} else {
			b->revert();
			upd->Ustar[bin(p, measureold)]++;
			d->Ustar[bin(p, measureold)]++;
			Uold += dUstar;
			avgp += b->Pnum();
			M[bin(p, measureold)] += b->Pnum();
			Mp[bin(p, measureold)][b->Pnum()]++;
			Mt[bin(p, measureold)]++;
			avgE[bin(p, measureold)] += elOld;

			if (b->Pnum() > pmax) pmax = b->Pnum();
		}
		i++;
		if (i == Nc) {
			delete upd;
			upd = requestupdateparameters(c, d);
			if (upd->key != d->key) {
				printf("keys do not match! requesting parameters\n");
				p = requestreturnparameters(c);
				d->key = upd->key;
				d->size = p->Nbins;
			}
			clear(d);
			d->dUstar = upd->dUstar;
			dUstar = upd->dUstar;
			Uold = penalty(p, upd, measureold);
			i = 0;
			j++;
		}
		if (j == 1000) {
			printf("Acceptance rate = %lf, nP = %d, avg p = %lf avg p (closed) = %lf\n", (double)Naccept/(double)Nt, b->Pnum(), (double)(avgp) / (double)(Nt), (double)(M[0])/(double)(Mt[0]));
			for (int i = 0; i < p->Nbins; i++) {
				printf("%d %lf %lf %ld\n", i+1, (double)M[i]/(double)Mt[i], avgE[i]/(double)Mt[i], Mt[i]);
			}
			printf("\n");
			j = 0;
		}
	}

}

