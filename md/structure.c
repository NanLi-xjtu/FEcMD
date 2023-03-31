#include "structure.h"

FILE *rdf;
FILE *ligancy;

SurfAtoms *surf;
double rRdf, nSurfCut, rSurfCut, rFlyCut, nFlyCut;
int countSurf, limitSurf, stepSurf, stepLigancy;

double *histRdf, rangeRdf;
int countRdf, limitRdf, sizeHistRdf, stepRdf;

int ligancy_flag = 0, flycount_flag = 0;

void ApplyBoundaryCondSurface()
{
	double vNew, vSignx, vSigny, vSignz, vvNew, vvOld;
	int n, i;
	enTransSum = 0.;
	DO_MOL {
		vSignx = 0.;
		vSigny = 0.;
		vSignz = 0.;
		if (mol[n].r.x >= 0.5 * region.x) vSignx = 1.;
		else if (mol[n].r.x < -0.5 * region.x) vSignx = -1.;
		if (mol[n].r.y >= 0.5 * region.y) vSigny = 1.;
		else if (mol[n].r.y < -0.5 * region.y) vSigny = -1.;
		if (mol[n].r.z >= 0.5 * region.z) vSignz = 1.;
		else if (mol[n].r.z < -0.5 * region.z) vSignz = -1.;
		if (vSignx != 0.) {
			mol[n].r.x = 0.49999 * vSignx * region.x;
//			vvOld = VLenSq (mol[n].rv); 
			i = (int)((mol[n].r.z / region.z + 0.5) * (double)sizeTempGrid);
			vNew = sqrt (NDIM * tempCell[i].temp);
			VRand3D (&mol[n].rv);
			VScale (mol[n].rv, vNew);
			vvNew = VLenSq (mol[n].rv);
//			enTransSum += 0.5 * vSignx * (vvNew - vvOld);
			if (mol[n].rv.x * vSignx > 0.) mol[n].rv.x *= -1.;
		}
		if (vSigny != 0.) {
			mol[n].r.y = 0.49999 * vSigny * region.y;
//			vvOld = VLenSq (mol[n].rv); 
			i = (int)((mol[n].r.z / region.z + 0.5) * (double)sizeTempGrid);
			vNew = sqrt (NDIM * tempCell[i].temp);
			VRand3D (&mol[n].rv);
			VScale (mol[n].rv, vNew);
			vvNew = VLenSq (mol[n].rv);
//			enTransSum += 0.5 * vSigny * (vvNew - vvOld);
			if (mol[n].rv.y * vSigny > 0.) mol[n].rv.y *= -1.;
		}
		if (vSignz != 0.) {
			mol[n].r.z = 0.49999 * vSignz * region.z;
			vvOld = VLenSq (mol[n].rv); 
			i = (int)((mol[n].r.z / region.z + 0.5) * (double)sizeTempGrid);
			vNew = sqrt (NDIM * tempCell[i].temp);
			VRand3D (&mol[n].rv);
			VScale (mol[n].rv, vNew);
			vvNew = VLenSq (mol[n].rv);
			enTransSum += 0.5 * vSignz * (vvNew - vvOld);
			if (mol[n].rv.z * vSignz > 0.) mol[n].rv.z *= -1.;
		}
	}
}

void PrintLigancy (FILE *fp)
{
	int n;

	fprintf (fp, "n     nebrnu\n");
	DO_MOL{
		fprintf(fp, "%5d %8.4f\n", n, surf[n].nebrnu);
	}
}

void SurfaceAtoms ()
{
	int n, m, N, i, id;
	double zmin, r;
	VecR dr;
	FILE *output;
	char filename[128];

	zmin = 1.e100;
	DO_MOL zmin = Min (zmin, mol[n].r.z);
	DO_MOL {
		if (surf[n].nebrnu < nSurfCut) {
			surf[n].surface = 1;
			if (fabs(mol[n].r.z - minzInit) < pedestal_thick) mol[n].flag = -1; //pedestal
			else mol[n].flag = 2; //surface
		} else {
			surf[n].surface = 0;
			mol[n].flag = 1; //bulk
		}
	}
	if (rFlyCut > rCut + rNebrShell) Error ("error(boundary.c): rFlyCut is too large\n");
	DO_MOL {
		N = 0;
		for (m = 0; m < mol[n].Nnebr; m ++) {
			id = mol[n].id_nebr[m];
			VSub (dr, mol[id].r, mol[n].r);
			if (strcmp(boundaryCond, "p") == 0) VWrapAll (dr);
			r = VLen (dr);
			if (r < rFlyCut) N ++;
		}
		if (N < nFlyCut + 0.000001) mol[n].flag = 3; //fly out of the surface
	}
	//count the number of atoms flying out of the surface
	if (boundary_fly == 1) {
		N = 0;
		DO_MOL {
			if (mol[n].flag == 3) N ++;
		}
		sprintf (filename, "out/md/atom_fly.dat");
		output = WriteFile (filename);
		if (flycount_flag == 0) {
			fprintf (output, "time(ps) Nflyatom\n");
			flycount_flag = 1;
		}
		if (N != 0) fprintf (output, "%f %d\n", TIMENOW_ps, N);
		fclose (output);
	}
}

void EvalLigancy ()
{
	VecR dr;
	double rr;
	int n, i, j1, j2;

	//calculate ligancy
	if (countSurf == 0) DO_MOL surf[n].nebrnu = 0;
	if (rSurfCut > rCut + rNebrShell) Error ("error(boundary.c): rSurfCut is too large\n");
	for (n = 0; n < nebrTabLen; n ++) {
		j1 = nebrTab[2 * n];
		j2 = nebrTab[2 * n + 1];
		VSub (dr, mol[j1].r, mol[j2].r);
		if (strcmp(boundaryCond, "p") == 0) VWrapAll (dr);
		rr = VLenSq (dr);
		if (rr < Sqr (rSurfCut)){
			++ surf[j1].nebrnu;
			++ surf[j2].nebrnu;
		}
	}
	DO_MOL mol[n].ligancy = surf[n].nebrnu;

	ligancy_flag = 1;
}

void MeasurementRdf (FILE *fp)
{
	double rb;
	int n;

	fprintf (fp, "d        rdf\n");
	for (n = 0; n < sizeHistRdf; n ++){
		rb = (n + 0.5) * rangeRdf / sizeHistRdf;
		fprintf (fp, "%8.4f %8.4f\n", rb, histRdf[n]);
	}

	//integral
	double intergral = 0, dx, x1 = 0, x2 = 1.59;
	int n1, n2;

	dx = rangeRdf / sizeHistRdf;
	n1 = x1 / dx - 0.5;
	n2 = x2 / dx - 0.5;
	for (n = n1; n <= n2; n ++){
		intergral += histRdf[n] * dx;
	}
//	printf ("\nRDF: r=%.4f to %.4f: intergral rdf= %.4f\n", x1, x2, intergral);

	//extremum
	int i, testnu = sizeHistRdf / 10 + 2, ifmax = 1, ifmin = 1;
	for (n = testnu; n < sizeHistRdf - testnu; n++){
		for (i = 1; i <= testnu; i ++){
			if (histRdf[n] <= histRdf[n-i] || histRdf[n] <= histRdf[n+i]) ifmax = 0;
			if (histRdf[n] >= histRdf[n-i] || histRdf[n] >= histRdf[n+i]) ifmin = 0;
		}
		if (ifmax == 1){ 
//			printf("maxr = %.4f, maxrdf = %.4f, ", (n + 0.5) * dx, histRdf[n]);
		}
		if (ifmin == 1){
//			printf("minr = %.4f, minrdf = %.4f, ", (n + 0.5) * dx, histRdf[n]);
			rRdf = (n + 0.5) * dx;
			break;
		}
		ifmax = 1;
		ifmin = 1;
	}
//	printf("\n\n");
}

void EvalRdf ()
{
	VecR dr;
	double deltaR, normFac, rr;
	int j1, j2, n;

	if (countRdf == 0) {
		for (n = 0; n < sizeHistRdf; n ++) histRdf[n] = 0.;
	}
	deltaR = rangeRdf / sizeHistRdf;
	for (j1 = 0; j1 < nMol - 1; j1 ++) {
		for (j2 = j1 + 1; j2 < nMol; j2 ++) {
			VSub (dr, mol[j1].r, mol[j2].r);
//			VWrapAll (dr);
			rr = VLenSq (dr);
			if (rr < Sqr (rangeRdf)) {
				n = sqrt (rr) / deltaR;
				++ histRdf[n];
			}
		}
	}
	++ countRdf;
	if (countRdf == limitRdf) {
		normFac = VProd (region) / (2. * M_PI * Cube (deltaR) * Sqr (nMol) * countRdf);
		for (n = 0; n < sizeHistRdf; n ++)
			histRdf[n] *= normFac / Sqr (n - 0.5);
		MeasurementRdf (rdf);
		countRdf = 0;
	}
}
