#include "integral.h"

double vMax_set, vMin_set;

void LeapfrogStep (int part)
{
	int n, i, flag;

	i = 0;
	if (part == 1) {
//		#pragma omp parallel for num_threads(Nthreads)
		DO_MOL{
			if (mol[n].pedest == 0) {
				VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
				if ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv) > vMax_set) {
					VScale(mol[n].rv, sqrt(vMax_set / ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv))));
				} else if ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv) < vMin_set) {
					VScale(mol[n].rv, sqrt(vMin_set / ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv))));
				}
				
				VVSAdd (mol[n].r, deltaT, mol[n].rv);
			}
		}
	} else {
//		#pragma omp parallel for num_threads(Nthreads)
		DO_MOL{
			if (mol[n].pedest == 0) {
				VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
				if ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv) > vMax_set) {
					VScale(mol[n].rv, sqrt(vMax_set / ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv))));
				} else if ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv) < vMin_set) {
					VScale(mol[n].rv, sqrt(vMin_set / ((1. / NDIM) * mol[n].mass * VLenSq (mol[n].rv))));
				}
			}
		}
	}
}

void PredictorStep ()
{
	double cr[] = {19., -10., 3.}, cv[] = {27., -22., 7.}, div = 24., wr, wv;
	int n;

	wr = Sqr (deltaT) / div;
	wv = deltaT / div;
	#pragma omp parallel for num_threads(Nthreads)
	DO_MOL {
		if (mol[n].pedest == 1) continue;
		mol[n].ro = mol[n].r;
		mol[n].rvo = mol[n].rv;
		PR (x);
		PRV (x);
		PR (y);
		PRV (y);
		PR (z);
		PRV (z);
		mol[n].ra2 = mol[n].ra1;
		mol[n].ra1 = mol[n].ra;
	}
}

void CorrectorStep ()
{
	double cr[] = {3., 10., -1.}, cv[] = {7., 6., -1.}, div = 24., wr, wv;
	int n;

	wr = Sqr (deltaT) / div;
	wv = deltaT / div;
	#pragma omp parallel for num_threads(Nthreads)
	DO_MOL {
		if (mol[n].pedest == 1) continue;
		CR (x);
		CRV (x);
		CR (y);
		CRV (y);
		CR (z);
		CRV (z);
	}
}

