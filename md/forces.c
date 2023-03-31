/****************************************************************************************************
** the different approaches to computing interactions: all pairs cell subdivision and neighbor lists.
*****************************************************************************************************/

#include "forces.h"

VecI cells;
double dispHi, rNebrShell, lj_r, lj_e;
int *cellList, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;

double embedWt, rSwitch, splineA2, splineA3;

//all pairs method for interaction computations
void ComputeForces()
{
	VecR dr;
	double fcVal, uVal, rr, rrCut, rri, rri3, r, o, e;
	int j1, j2, n;

	o = lj_r; //--A
	o /= lUnit * 1.e10; //dimensionless
	e = lj_e * eleChar; //--J
	e /= eUnit; //dimensionless

	rrCut = Sqr (rCut);
	DO_MOL VZero (mol[n].f);
	uSum = 0.;
	virSum = 0.;
	for (j1 = 0; j1 < nMol - 1; j1 ++) {
		for (j2 = j1 + 1; j2 < nMol; j2 ++) {
			VSub (dr, mol[j1].r, mol[j2].r);
			if (strcmp(boundaryCond, "p") == 0) VWrapAll (dr);
			rr = VLenSq (dr);
			if (rr < rrCut){
				r = sqrt (rr);
				rri = 1. / rr;
				rri3 = Cube (rri);
				fcVal = 48. * e / (o * o) * (pow(o/r, 14.) - 0.5 * pow(o/r, 8.));
				uVal = 4. * e * (pow(o/r, 12.) - pow(o/r, 6.)) + e;
				VVSAdd (mol[j1].f, fcVal, dr);
				VVSAdd (mol[j2].f, - fcVal, dr);
				uSum += uVal;
				virSum += fcVal * rr;
			}
		}
	}
}

//cell subdivision
void ComputeForcesCell ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double fcVal, rr, rrCut, rri, rri3, uVal, r, o, e;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

	o = 3.4; //--A
	o /= lUnit * 1.e10; //dimensionless
	e = 1.65672e-21; //--J
	e /= eUnit; //dimensionless

	rrCut = Sqr (rCut);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n++) cellList[n] = -1;
	DO_MOL{
		VSAdd (rs, mol[n].r, 0.5, region);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;

		if (c >= (VProd(cells) + nMol) || c < 0) {
			printf ("error(forces.c): cellList is not big enough,\n \
				cc:(%d %d %d), cells(%d %d %d)\n \
				rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
				rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
				region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
			PrintMovie ();
			exit (1);
		}

		cellList[n] = cellList[c];
		cellList[c] = n;
	}

	DO_MOL VZero (mol[n].f);
	uSum = 0.;
	virSum = 0.;
	for (m1z = 0; m1z < cells.z; m1z ++){
		for (m1y = 0; m1y < cells.y; m1y ++){
			for (m1x = 0; m1x < cells.x; m1x ++){
				VSet (m1v, m1x, m1y, m1z);
				m1 = VLinear (m1v, cells) + nMol;
				for (offset = 0; offset < N_OFFSET; offset ++){
					VAdd (m2v, m1v, vOff[offset]);
					VZero (shift);
					if (strcmp(boundaryCond, "n") == 0) {
						if (m2v.z < 0 || m2v.z >= cells.z) continue;
						if (m2v.y < 0 || m2v.y >= cells.y) continue;
						if (m2v.x < 0 || m2v.x >= cells.x) continue;
					} else if (strcmp(boundaryCond, "p") == 0)
						VCellWrapAll ();

					m2 = VLinear (m2v, cells) + nMol;
					DO_CELL (j1, m1){
						DO_CELL (j2, m2){
							if (m1 != m2 || j2 < j1){
								VSub (dr, mol[j1].r, mol[j2].r);
								if (strcmp(boundaryCond, "p") == 0) VVSub (dr, shift);
								rr = VLenSq (dr);
								if (rr < rrCut){
									r = sqrt (rr);
									rri = 1. / rr;
									rri3 = Cube (rri);
									fcVal = 48. * e / (o * o) * (pow(o/r, 14.) - 0.5 * pow(o/r, 8.));
									uVal = 4. * e * (pow(o/r, 12.) - pow(o/r, 6.)) + e;
									VVSAdd (mol[j1].f, fcVal, dr);
									VVSAdd (mol[j2].f, - fcVal, dr);
									uSum += uVal;
									virSum += fcVal * rr;
								}
							}
						}
					}
				}
			}
		}
	}
}

void ComputeForcesNebr ()
{
	VecR dr;
	double fcVal, rr, rrCut, rri, rri3, uVal, r, o, e;
	int j1, j2, n;

	o = lj_r; //--A
	o /= lUnit * 1.e10; //dimensionless
	e = lj_e * eleChar; //--J
	e /= eUnit; //dimensionless

	rrCut = Sqr (rCut);
	DO_MOL VZero (mol[n].f);
	uSum = 0.;
	virSum = 0.;
	for (n = 0; n < nebrTabLen; n ++){
		j1 = nebrTab[2 * n];
		j2 = nebrTab[2 * n + 1];
		VSub (dr, mol[j1].r, mol[j2].r);
		if (strcmp(boundaryCond, "p") == 0) VWrapAll (dr);
		rr = VLenSq (dr);
		if (rr < rrCut){
			r = sqrt (rr);
			rri = 1. / rr;
			rri3 = Cube (rri);
			fcVal = 48. * e / (o * o) * (pow(o/r, 14.) - 0.5 * pow(o/r, 8.));
			uVal = 4. * e * (pow(o/r, 12.) - pow(o/r, 6.)) + e;
			VVSAdd (mol[j1].f, fcVal, dr);
			VVSAdd (mol[j2].f, -fcVal, dr);
			uSum += uVal;
			virSum += fcVal * rr;
		}
	}
}

//EAM potential
void ComputeForcesEam ()
{
	VecR dr;
	double eDim, fcVal, rCutC, rr, rrCdi, rrCut, rrd, rri, rri3,
	rrSwitch, t, uVal;
	int j1, j2, n;
	eDim = NDIM * (NDIM + 1) * exp (1.);
	rCutC = pow (2., 1./6.);
	rrCut = Sqr (rCut);
	rrCdi = 1. / Sqr (rrCut - Sqr (rCutC));
	rrSwitch = Sqr (rSwitch);
	DO_MOL mol[n].logRho = 0.;
	for (n = 0; n < nebrTabLen; n ++) {
		j1 = nebrTab[2 * n];
		j2 = nebrTab[2 * n + 1];
		VSub (dr, mol[j1].r, mol[j2].r);
		VWrapAll (dr);
		rr = VLenSq (dr);
		if (rr < rrCut) {
			t = Sqr (rrCut - rr);
			mol[j1].logRho += t;
			mol[j2].logRho += t;
		}
	}
	DO_MOL {
		if (mol[n].logRho > 0.) {
			mol[n].logRho = log ((rrCdi / eDim) * mol[n].logRho);
		}
	}
	DO_MOL VZero (mol[n].ra);
	uSum = 0.;
	for (n = 0; n < nebrTabLen; n ++) {
		j1 = nebrTab[2 * n];
		j2 = nebrTab[2 * n + 1];
		VSub (dr, mol[j1].r, mol[j2].r);
		VWrapAll (dr);
		rr = VLenSq (dr);
		if (rr < rrCut) {
			rrd = rrCut - rr;
			if (rr < rrSwitch) {
				rri = 1. / rr;
				rri3 = Cube (rri);
				fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
//				fcVal = 48. * rri3 * (rri3 - 0.5) * rri - 2. * sqrt (rri);
				uVal = 4. * rri3 * (rri3 - 1.);
//				uVal = 4. * rri3 * (rri3 - 1.) + 2. * sqrt (rr);
			} else {
				fcVal = (4. * splineA2 + 6. * splineA3 * rrd) * rrd;
				uVal = (splineA2 + splineA3 * rrd) * Sqr (rrd);
			}
			fcVal = embedWt * fcVal + (1. - embedWt) * 2. * rrCdi * (mol[j1].logRho + mol[j2].logRho + 2.) * rrd;
			VVSAdd (mol[j1].ra, fcVal, dr);
			VVSAdd (mol[j2].ra, - fcVal, dr);
			uSum += uVal;
		}
	}
	t = 0.;
	DO_MOL t += mol[n].logRho * exp (mol[n].logRho);
	uSum = embedWt * uSum + (1. - embedWt) * 0.5 * eDim * t;
}

//EAM params
void EvalEamParams ()
{
	double bb, p, pd, rr, rr3;
	rSwitch = pow (26. / 7., 1. / 6.);
	rr = Sqr (rSwitch);
	rr3 = Cube (rr);
	p = 4. * (1. / rr3 - 1.) / rr3;
	pd = - 48. * (1. / rr3 - 0.5) / (rSwitch * rr3);
//	p = 4. * (1. / rr3 - 1.) / rr3 + 2. * sqrt (rr);
//	pd = - 48. * (1. / rr3 - 0.5) / (rSwitch * rr3) + 2. / sqrt (rr);
	bb = 4. * (1. - sqrt (1. + 3. * p / (2. * rSwitch * pd)));
	splineA2 = (6. * p + bb * rSwitch * pd) / (2. * Sqr (bb * rr));
	splineA3 = - (4. * p + bb * rSwitch * pd) / (2. * Sqr (bb * rr) * bb * rr);
	rCut = rSwitch * sqrt (bb + 1.);
	embedWt = 0.3333;
}

