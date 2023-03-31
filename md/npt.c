#include "npt.h"

double extPressure, g1Sum, g2Sum, massS, massV, varS, varSa,
     varSa1, varSa2, varSo, varSv, varSvo, varV, varVa, varVa1,
     varVa2, varVo, varVv, varVvo;
int maxEdgeCells;

void ComputeDerivsPT ()
{
	double aFac, vFac, EkSum;
	int n;

	EkSum = 0.;
	DO_MOL EkSum += 0.5 * mol[n].mass * VLenSq (mol[n].rv);
	EkSum *= pow (varV, 2./3.);
	g1Sum = 2. * EkSum - 3. * nMol * temperature;
//	g2Sum = 2. * EkSum + virSum - 3. * extPressure * varV;
	g2Sum = 2. * EkSum + virSum - 3. * extPressure * nMol / density;
	aFac = pow (varV, -1./3.);
	vFac = - varSv / varS - 2. * varVv / (3. * varV);
	DO_MOL VSSAdd (mol[n].ra, aFac, mol[n].ra, vFac, mol[n].rv);
	varSa = Sqr (varSv) / varS + g1Sum * varS / massS;
	varVa = varSv * varVv / varS + g2Sum * Sqr (varS) /
	(3. * massV * varV);
}

void PredictorStepPT ()
{
	double cr[] = {19., -10., 3.}, cv[] = {27., -22., 7.}, div = 24., e, wr, wv;

	wr = Sqr (deltaT) / div;
	wv = deltaT / div;

	varSo = varS;
	varSvo = varSv;
	PCR4_NPT (varS, varS, varSv, varSa, varSa1, varSa2);
	PCV4_NPT (varS, varSo, varSv, varSa, varSa1, varSa2);
	varSa2 = varSa1;
	varSa1 = varSa;
	//... (ditto for varV) ...
	varVo = varV;
	varVvo = varVv;
	PCR4_NPT (varV, varV, varVv, varVa, varVa1, varVa2);
	PCV4_NPT (varV, varVo, varVv, varVa, varVa1, varVa2);
	varVa2 = varVa1;
	varVa1 = varVa;

	e = pow (varV, 1. / NDIM);
	latticePara *= (e / region.x);
	density *= VProd (region) / varV;
	VSetAll (region, e);
}

void CorrectorStepPT ()
{
	double cr[] = {3., 10., -1.}, cv[] = {7., 6., -1.}, div = 24., e, wr, wv;

	wr = Sqr (deltaT) / div;
	wv = deltaT / div;
	PCR4_NPT (varS, varSo, varSvo, varSa, varSa1, varSa2);
	PCV4_NPT (varS, varSo, varSvo, varSa, varSa1, varSa2);
	//... (ditto for varV) ...
	PCR4_NPT (varV, varVo, varVvo, varVa, varVa1, varVa2);
	PCV4_NPT (varV, varVo, varVvo, varVa, varVa1, varVa2);

	e = pow (varV, 1. / NDIM);
	latticePara *= (e / region.x);
	density *= VProd (region) / varV;
	VSetAll (region, e);
}

void InitFeedbackVars ()
{
	varS = 1.;
	varV = Cube (region.x);
	varSv = 0.;
	varSa = 0.;
	varSa1 = 0.;
	varSa2 = 0.;
	//... (ditto for varV...) ...
	varVv = 0.;
	varVa = 0.;
	varVa1 = 0.;
	varVa2 = 0.;
}

void ScaleCoords ()
{
	double fac;
	int n;

	fac = pow (varV, -1. / 3.);
	DO_MOL VScale (mol[n].r, fac);
}

void UnScaleCoords ()
{
	double fac;
	int n;

	fac = pow (varV, 1. / 3.);
	DO_MOL VScale (mol[n].r, fac);
}

void ScaleVels ()
{
	double fac;
	int n;

	fac = pow (varV, -1. / 3.);
	DO_MOL VScale (mol[n].rv, fac);
}

void UnScaleVels ()
{
	double fac;
	int n;

	fac = pow (varV, 1. / 3.);
	DO_MOL VScale (mol[n].rv, fac);
}

void UpdateCellSize ()
{
	VSCopy (cells, 1. / rCut, region);
	cells.x = Min (cells.x, maxEdgeCells);
	cells.y = Min (cells.y, maxEdgeCells);
	cells.z = Min (cells.z, maxEdgeCells);
}
