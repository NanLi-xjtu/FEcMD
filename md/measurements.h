#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "main.h"

#define PropZero(v)			\
	v.sum = 0.,			\
	v.sum2 = 0.
#define PropAccum(v)			\
	v.sum += v.val,			\
	v.sum2 += Sqr (v.val)
#define PropAvg(v,n)			\
	v.sum /= n,			\
	v.sum2 = sqrt ( Max (v.sum2 / n - Sqr (v.sum), 0.))
#define PropEst(v)			\
	v.sum, v.sum2

extern double *histVel, rangeVel;
extern int countVel, limitVel, sizeHistVel, stepVel;

void EvalProps ();
void AccumProps(int icode);
void EvalVelDist ();

#endif
