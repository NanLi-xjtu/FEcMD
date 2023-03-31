#include "heattransport.h"

char nonq_mode[128];
double gravField;
VecI sizeHistGrid;
double **histGrid;
int countGrid, limitGrid, stepGrid;
double *profileT, *profileV;
double enTransSum, wallTempHi, wallTempLo;
Prop thermalCond;

void BuildNebrListHeat ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

	rrNebr = Sqr (rCut + rNebrShell);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
	DO_MOL{
		VSAdd (rs, mol[n].r, 0.5, region);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;
		cellList[n] = cellList[c];
		cellList[c] = n;
	}
	nebrTabLen = 0;
	for (m1z = 0; m1z < cells.z; m1z ++){
		for (m1y = 0; m1y < cells.y; m1y ++){
			for (m1x = 0; m1x < cells.x; m1x ++){
				VSet (m1v, m1x, m1y, m1z);
				m1 = VLinear (m1v, cells) + nMol;
				for (offset = 0; offset < N_OFFSET; offset ++){
					VAdd (m2v, m1v, vOff[offset]);
					VZero (shift);
//					VCellWrapAll ();

					//z boundaries are no longer periodic
					VCellWrap (x);
					VCellWrap (y);
					if (m2v.z < 0 || m2v.z >= cells.z) continue;

					m2 = VLinear (m2v, cells) + nMol;
					DO_CELL (j1, m1){
						DO_CELL (j2,m2){
							if (m1 != m2 || j2 < j1){
								VSub (dr, mol[j1].r, mol[j2].r);
								VVSub (dr, shift);
								if (VLenSq (dr) < rrNebr){
									if (nebrTabLen >= nebrTabMax)
										ErrExit (ERR_TOO_MANY_NEBRS);
									nebrTab[2 * nebrTabLen] = j1;
									nebrTab[2 * nebrTabLen + 1] = j2;
									++ nebrTabLen;
								}
							}
						}
					}
				}
			}
		}
	}
}

void ComputeExternalForce ()
{
	int n;

	DO_MOL mol[n].ra.x += gravField;
}

void ApplyBoundaryCondFlow ()
{
	double vSign;
	int n;

	DO_MOL {
		VWrap (mol[n].r, x);
		VWrap (mol[n].r, y);
		vSign = 0.;
		if (mol[n].r.z >= 0.5 * region.z) vSign = 1.;
		else if (mol[n].r.z < -0.5 * region.z) vSign = -1.;
		if (vSign != 0.) {
			mol[n].r.z = 0.49999 * vSign * region.z;
			VRand3D (&mol[n].rv);
			VScale (mol[n].rv, velMag);
			if (mol[n].rv.z * vSign > 0.) mol[n].rv.z *= -1.;
		}
	}
}

void ApplyBoundaryCondHeat()
{
	double vNew, vSign, vvNew, vvOld;
	int n;
	enTransSum = 0.;
	DO_MOL {
		VWrap (mol[n].r, x);
		VWrap (mol[n].r, y);
		vSign = 0.;
		if (mol[n].r.z >= 0.5 * region.z) vSign = 1.;
		else if (mol[n].r.z < -0.5 * region.z) vSign = -1.;
		if (vSign != 0.) {
			mol[n].r.z = 0.49999 * vSign * region.z;
			vvOld = VLenSq (mol[n].rv); 
			vNew = sqrt (NDIM * ((vSign < 0.) ? wallTempHi : wallTempLo) / mol[n].mass);
			VRand3D (&mol[n].rv);
			VScale (mol[n].rv, vNew);
			vvNew = VLenSq (mol[n].rv);
			enTransSum += 0.5 * vSign * mol[n].mass * (vvNew - vvOld);
			if (mol[n].rv.z * vSign > 0.) mol[n].rv.z *= -1.;
		}
	}
}

void GridAverage (int opCode)
{
	VecR invWid, rs, va;
	VecI cc;
	double pSum; 
	int c, hSize, j, n;

	hSize = VProd (sizeHistGrid);
	if (opCode == 0) {
		for (j = 0; j < NHIST; j ++) { 
			for (n = 0; n < hSize; n ++) histGrid[j][n] = 0.;
		}
	} else if (opCode == 1) {
		VDiv (invWid, sizeHistGrid, region);
		DO_MOL { 
			VSAdd (rs, mol[n].r, 0.5, region);
			VMul (cc, rs, invWid);
			c = VLinear (cc, sizeHistGrid);
			++ histGrid[0][c];
			histGrid[1][c] += mol[n].mass * VLenSq (mol[n].rv); 
			histGrid[2][c] += mol[n].rv.x;
			histGrid[3][c] += mol[n].rv.y;
			histGrid[4][c] += mol[n].rv.z;
		}
	} else if (opCode == 2) { 
		pSum = 0.;
		for (n = 0; n < hSize; n ++) {
			if (histGrid[0][n] > 0.) {
				for (j = 1; j < NHIST; j ++) histGrid[j][n] /= histGrid[0][n];
				VSet (va, histGrid[2][n], histGrid[3][n], histGrid[4][n]); 
				histGrid[1][n] = (histGrid[1][n] - VLenSq (va)) / NDIM;
				pSum += histGrid[0][n];
			}
		}
		pSum /= hSize; 
		for (n = 0; n < hSize; n ++) histGrid[0][n] /= pSum;
	}
}

//The functions that extract the profiles from the grid data and output the results follow.
void EvalProfile ()
{
	int k, n;

	for (n = 0; n < sizeHistGrid.z; n ++) {
		profileV[n] = 0.;
		profileT[n] = 0.;
	}
	for (n = 0; n < VProd (sizeHistGrid); n ++) {
		k = n / (sizeHistGrid.x * sizeHistGrid.y);
		profileV[k] += histGrid[2][n];
		profileT[k] += histGrid[1][n];
	}
	for (n = 0; n < sizeHistGrid.z; n ++) {
		profileV[n] /= sizeHistGrid.x * sizeHistGrid.y;
		profileT[n] /= sizeHistGrid.x * sizeHistGrid.y;
	}
}

void PrintProfile ()
{
	double zVal;
	int n;

	if ((velprofile = fopen ("out/md/V_profile.dat", "a+")) == NULL){
		printf ("\nopen V_profile file error");
		getchar ();
		exit (1);
	}
	fprintf (velprofile, "position flowSpeed(km/s)\n");
	for (n = 0; n < sizeHistGrid.z; n ++) {
		zVal = (n + 0.5) / sizeHistGrid.z;
		fprintf (velprofile, "%.2f %.3f\n", zVal, profileV[n] * vUnit);
	}
	fclose (velprofile);

	if ((temprofile = fopen ("out/md/T_profile.dat", "a+")) == NULL){
		printf ("\nopen T_profile file error");
		getchar ();
		exit (1);
	}
	fprintf (temprofile, "position temperature(K)\n");
	for (n = 0; n < sizeHistGrid.z; n ++) {
		zVal = (n + 0.5) / sizeHistGrid.z;
		fprintf (temprofile, "%.2f %.3f\n", zVal, profileT[n] * TUnit);
	}
	fclose (temprofile);
}
