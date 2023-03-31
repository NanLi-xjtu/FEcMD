#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "main.h"

#define PCR4(r, ro, v, a, a1, a2, t) \
	r.t = ro.t + deltaT * v.t + wr * (cr[0] * a.t + \
	      cr[1] * a1.t + cr[2] * a2.t)
#define PCV4(r, ro, v, a, a1, a2, t) \
	v.t = (r.t - ro.t) / deltaT + wv * (cv[0] * a.t + \
	      cv[1] * a1.t + cv[2] * a2.t)

#define PR(t) \
	PCR4 (mol[n].r, mol[n].r, mol[n].rv, mol[n].ra, \
	      mol[n].ra1, mol[n].ra2, t)
#define PRV(t) \
	PCV4 (mol[n].r, mol[n].ro, mol[n].rv, mol[n].ra, \
	      mol[n].ra1, mol[n].ra2, t)
#define CR(t) \
	PCR4 (mol[n].r, mol[n].ro, mol[n].rvo, mol[n].ra, \
	      mol[n].ra1, mol[n].ra2, t)
#define CRV(t) \
	PCV4 (mol[n].r, mol[n].ro, mol[n].rv, mol[n].ra, \
	      mol[n].ra1, mol[n].ra2, t)

extern double vMax_set, vMin_set;

void LeapfrogStep (int part);
void PredictorStep ();
void CorrectorStep ();

#endif
