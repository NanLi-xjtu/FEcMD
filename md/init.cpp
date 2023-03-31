#include "init.h"

void SetParams (int argc, char **argv)
{
	FILE *input;
	char line[1024], filename[128], eam_elem[20][10];
	double temp;
	int n;

	//Unit
	lUnit = 1.e-10; //--m
	mUnit = M_C12; //--kg
	eUnit = 300. * kB; //--J
	tUnit = sqrt((mUnit * Sqr(lUnit)) / eUnit); //--s
	TUnit = eUnit / kB; //--K
	vUnit = sqrt(eUnit / mUnit); //--m/s
	fUnit = eUnit / lUnit; //--N
	PUnit = fUnit / Sqr (lUnit); //--Pa
	aUnit = fUnit / mUnit; //--m/s2

	//read in atom number and region from "./in/**"
	sprintf (mdXYZFilename, "in/mdlat.in.xyz");
	input = ReadFile (mdXYZFilename);
	fgets(line, 1024, input);
	sscanf (line, "%d", &nMol);
	fgets(line, 1024, input);
	sscanf (line, "Lattice=\"%lg %lg %lg %lg %lg %lg %lg %lg %lg", \
		&region.x, &temp, &temp, &temp, &region.y, &temp, &temp, &temp, &region.z); //--A
	VScale(region, 1. / (lUnit * 1.e10)); //--dimensionless
	if (nMol == 0 || fabs(VProd(region)) < 1.e-100) {
		printf ("error(init.c): in/md.xyz input format\n");
		exit (1);
	}
	fclose (input);

	deltaT = DELTAT_fs / (tUnit / 1.e-15); //--dimensionless
	latticePara = LATTIPARA_A / (lUnit * 1.e10); //--dimensionless
	if (strcmp(structure_type, "FCC") == 0) density = 4. / pow(latticePara, 3.); //--dimensionless
	else if (strcmp(structure_type, "BCC") == 0) density = 2. / pow(latticePara, 3.); //--dimensionless
	else if (strcmp(structure_type, "cubic") == 0) density = 1. / pow(latticePara, 3.); //--dimensionless
	else {
//		printf ("error(init.c): structure_type format\n");
//		exit (1);
	}
	temperature = TEMPERATURE_K / TUnit;
	vMax_set /= TUnit;
	vMin_set /= TUnit;

	if (strcmp(boundaryCond, "p") == 0) {
//		if (fabs (pow (nMol / density, 1. / 3.) - pow (VProd (region), 1. / 3.)) > 0.1 / (lUnit * 1.e10)) {
//			printf ("error(init.c): LATTIPARA_A or structure_type\n");
//			exit (1);
//		}
//		density = nMol / VProd (region);
	}

	//read potential parameters
	if (strcmp (force_type, "metal") == 0) {
		Coeff ();//read metal EAM potential file
		init_style(); //metal-convert read-in file(s) to arrays and spline them
		printfeam (0);//print eam potential from files
		printfspline (phi_spline);//frho_spline; rhor_spline; z2r_spline; phi_spline;//print eam spline
		printf ("%d mass = %.3f, lattice constant = %.3f A, %s\n", \
			funcfl->id, funcfl->mass, funcfl->lattice, funcfl->structure);
	} else if (strcmp (force_type, "alloy") == 0) {
		coeff_alloy(eamFilename); //read alloy EAM potential file
		init_style_alloy (); //alloy-convert read-in file(s) to arrays and spline them
		printfeam_alloy ();
		printfspline_alloy (rhor_spline);//frho_spline; rhor_spline; z2r_spline; phi_spline;
		for (n = 0; n < setfl->nelements; n ++) {
			printf ("%s mass = %.3f, lattice constant = %.3f A, %s\n", \
				setfl->elements[n], setfl->mass[n], setfl->lattice[n], setfl->structure[n]);
		}
	} else if (strcmp (force_type, "eamfs") == 0) {
		coeff_eamfs(); //read alloy EAM potential file
		init_style_eamfs (); //alloy-convert read-in file(s) to arrays and spline them
		printfeam_alloy ();
		for (n = 0; n < fs->nelements; n ++) {
			printf ("%s mass = %.3f, lattice constant = %.3f A, %s\n", \
				fs->elements[n], fs->mass[n], fs->lattice[n], fs->structure[n]);
		}
	} else if (strcmp (force_type, "lj") == 0) {
			cutmax = 3.8; //--A
			printf ("Lennard-Jones (LJ) potential\n");
	} else if (strcmp (force_type, "snap") == 0) {
		PairSNAP_coeff();
		PairSNAP_init_style();
		ComputeSNAAtom ();
		
		cutmax = sna_rcutmax; //--A
	} else if (strcmp (force_type, "MTP") == 0) {
			mlp.settings ();
			mlp.coeff();
			mlp.init_style ();
			cutmax = mlp.init_one (1, 1); //--A
	}

	rCut = cutmax / (lUnit * 1.e10); //--dimensionless

	//set cell size and box size
	if (strcmp (ensemble, "NPT") == 0) sprintf (interact_method, "cell");
	rNebrShell /= lUnit * 1.e10; //dimensionless
	if (isnormal (VProd(regionScale)) == 1) VVMul (region, regionScale);
	if (strcmp (interact_method, "nebr") == 0) {
		VSCopy (cells, 1. / (rCut + rNebrShell), region);
	} else if (strcmp (interact_method, "cell") == 0)
		VSCopy (cells, 1. / rCut, region);
	nebrTabMax = nebrTabFac * nMol;

	//simple EAM params
//	EvalEamParams ();

	//the initial velocity params
	velMag = sqrt (NDIM * (1. - 1. / nMol) * temperature);

	//NPT
	maxEdgeCells = 1.3 * cells.x;
	extPressure /= (PUnit / 1.e9); //--dimensionless

	//ligancy cut radius
	rSurfCut /= (lUnit * 1.e10); //--dimensionless
	rFlyCut /= (lUnit * 1.e10); //--dimensionless

	//Pedestal
	pedestal_thick /= lUnit * 1.e10; //dimensionless

	//heat transport measurement
	wallTempHi /= TUnit; //dimensionless
	wallTempLo /= TUnit; //dimensionless

	//Maxwell strss
	Maxwell_stress = Maxwell_begin;
}

void SetupJob(int argc, char **argv)
{
	int n, i;
	char filename[128];

	printf("\n\n-- init\n");
	InitRand (randSeed);
	AllocArrays ();
	AllocArraysMD ();
	stepCount = 0;

	/*******************************************************
	** InitCoords (); -- simple cubic lattice
	** InitCoordsFCC (); -- FCC lattice
	** InitCoordsReadin (); --read in lattice
	*******************************************************/
	InitCoordsReadin ();
	//DO_MOL printf ("%s %f %f %f\n", mol[n].elem, mol[n].r.x, mol[n].r.y, mol[n].r.z);
	//set atom mass according to the element and eam potential file
	if (strcmp (force_type, "metal") == 0) {
		DO_MOL mol[n].mass = funcfl->mass * 1.6606 * 1.e-27 / mUnit; //--dimensionless
	} else if (strcmp (force_type, "alloy") == 0) {
		DO_MOL {
			for (i = 0; i < setfl->nelements; i ++)
				if (strcmp(mol[n].elem, setfl->elements[i]) == 0) 
					mol[n].mass = setfl->mass[i] * 1.6606 * 1.e-27 / mUnit; //--dimensionless
		}
	} else if (strcmp (force_type, "eamfs") == 0) {
		DO_MOL {
			for (i = 0; i < fs->nelements; i ++)
				if (strcmp(mol[n].elem, fs->elements[i]) == 0) 
					mol[n].mass = fs->mass[i] * 1.6606 * 1.e-27 / mUnit; //--dimensionless
		}
	} else if (strcmp (force_type, "lj") == 0) {
		DO_MOL mol[n].mass = 39.95 * M_C12 / mUnit; //--dimensionless, Ar atom
	} else if (strcmp (force_type, "snap") == 0 || strcmp (force_type, "MTP") == 0) {
		DO_MOL {
			for (i = 0; i < Nelems; i ++)
				if (strcmp(mol[n].elem, elems[i]) == 0) 
					mol[n].mass = mass[i] * 1.6606 * 1.e-27 / mUnit; //--dimensionless
		}
	}

	DO_MOL {
		if (fabs(mol[n].mass) < 1.e-100) {
			if (argc > 1) {
				if (strcmp (argv[1], "ED") != 0) {
					printf ("error(init.c): set atom mass %s: %e\n", mol[n].elem, mol[n].mass);
					exit (1);
				}
			}
		}
	}

	AllocArraysNebr ();

	//find pedestal atoms' id -- pedest[]
	FindPedestal ();

	InitVels ();
	InitAccels ();
	DO_MOL {
		VZero(mol[n].F);
		VZero(mol[n].f);
		mol[n].flag = 0;
		mol[n].ligancy = 0;
	}

	AccumProps (0);
	
	//neigbor list params
	nebrNow = 1;

	//heat transport measurement initialization
	GridAverage (0);
	countGrid = 0;

	//tempearture adjust
	kinEnInitSum = 0;

	//scaling
	stepSca = 0;

	//surface
	countSurf = 0;

	//RDF
	rRdf = 0;
	countRdf = 0;

	//velocity distribution
	rangeVel /= vUnit*1.e-2;//dimesionless
	countVel = 0;

	//atom fly
	atomFlyCount = 0;
	DO_MOL mol[n].deleteFlag = 0;

	/*************************Print***************************/
	if ((summary = fopen ("out/md/summary.dat", "a+")) == NULL) {
		printf ("\nopen summary file error");
		getchar ();
		exit (1);
	}
	fprintf(summary, "SetParams -- data\n");
	fprintf(summary, "lUnit = %.4f A, mUnit = %.4e kg, eUnit = %.4f eV, tUnit = %.4f ps, TUnit = %.4f K, PUnit = %.4e Pa, \
		\nvUnit = %.4e m/s, fUnit = %.4e N, aUnit = %.4e m/s2\n", \
		lUnit * 1.e10, mUnit, eUnit / eleChar, tUnit * 1.e12, TUnit, PUnit, vUnit, fUnit, aUnit);
	fprintf(summary, "deltaT = %.4f ps, temperature = %.4f K, rCut = %.4f A, \
		\nregion = [%.4f %.4f %.4f] A, nMol = %d, cells = [%d %d %d]\n", \
		deltaT * tUnit * 1.e12, temperature * TUnit, rCut * lUnit * 1.e10, \
		region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10, \
		nMol, cells.x, cells.y, cells.z);
	fclose (summary);

	//print the initial state: stepcount = 0
	PrintMovie ();

	if (strcmp (ensemble, "NPT") == 0) {
		InitFeedbackVars ();
		ScaleCoords ();
		ScaleVels ();
	}
}

void AllocArrays ()
{
	
}

void AllocArraysMD()
{
	AllocMem (mol, nMol, Mol);
	AllocMem (histRdf, sizeHistRdf, double);
	AllocMem (histVel, sizeHistVel, double);
}

void AllocArraysNebr ()
{
	int n, k;

	if (strcmp (ensemble, "NPT") == 0) AllocMem (cellList, Cube (maxEdgeCells) + nMol, int);
	else AllocMem (cellList, VProd (cells) + nMol, int);

	AllocMem (nebrTab, 2 * nebrTabMax, int);
	AllocMem2 (histGrid, NHIST, VProd (sizeHistGrid), double);
	AllocMem (profileT, sizeHistGrid.z, double);
	AllocMem (profileV, sizeHistGrid.z, double);
	AllocMem (tempCell, sizeTempGrid, TempCell);
	AllocMem (surf, nMol, SurfAtoms);
	DO_MOL AllocMem (mol[n].id_nebr, 2 * nebrTabFac, int);
}

void InitCoords()
{
	VecR c, gap;
	int n, nx, ny, nz;

	VDiv(gap, region, initUcell);
	n = 0;
	for (nz = 0; nz < initUcell.z; nz ++){
		for (ny = 0; ny < initUcell.y; ny ++){
			for (nx = 0; nx < initUcell.x; nx ++){
				VSet (c, nx + 0.5, ny + 0.5, nz + 0.5);
				VMul (c, c, gap);
				VVSAdd (c, -0.5, region);
				mol[n].r = c;
				++ n;
			}
		}
	}
}

void InitCoordsFCC ()
{
	VecR c, gap;
	int j, n, nx, ny, nz;

	VDiv (gap, region, initUcell);
	n = 0;
	for (nz = 0; nz < initUcell.z; nz ++) {
		for (ny = 0; ny < initUcell.y; ny ++) {
			for (nx = 0; nx < initUcell.x; nx ++) {
				VSet (c, nx + 0.25, ny + 0.25, nz + 0.25);
				VMul (c, c, gap);
				VVSAdd (c, -0.5, region);
				for (j = 0; j < 4; j ++) {
					mol[n].r = c;
					if (j != 3) {
						if (j != 0) mol[n].r.x += 0.5 * gap.x;
						if (j != 1) mol[n].r.y += 0.5 * gap.y;
						if (j != 2) mol[n].r.z += 0.5 * gap.z;
					}
					++ n;
				}
			}
		}
	}
}

void InitCoordsReadin ()
{
	FILE *input;
	char line[1024];
	int n, i;
	VecR rMin, rMax, rCenter;
	Mol *molNew;

	input = ReadFile (mdXYZFilename);
	fgets(line, 1024, input);
	fgets(line, 1024, input);
	n = 0;
	while (1) {
		if (fgets(line, 1024, input) == NULL) break;
		sscanf (line, "%s %lg %lg %lg", mol[n].elem, &mol[n].r.x, &mol[n].r.y, &mol[n].r.z); //--A
//		printf ("%s %f %f %f\n", mol[n].elem, mol[n].r.x, mol[n].r.y, mol[n].r.z);
		VScale(mol[n].r, 1. / (lUnit * 1.e10)); //--dimensionless
		n ++;
	}
	//check
	if (n != nMol) {
		printf ("error: read in atom coords(init.c)\n");
		exit (1);
	}
	fclose (input);

	VSet(rMin, 1.e100, 1.e100, 1.e100);
	VSet(rMax, -1.e100, -1.e100, -1.e100);
	DO_MOL {
		if (rMin.x > mol[n].r.x) rMin.x = mol[n].r.x;
		if (rMin.y > mol[n].r.y) rMin.y = mol[n].r.y;
		if (rMin.z > mol[n].r.z) rMin.z = mol[n].r.z;
		if (rMax.x < mol[n].r.x) rMax.x = mol[n].r.x;
		if (rMax.y < mol[n].r.y) rMax.y = mol[n].r.y;
		if (rMax.z < mol[n].r.z) rMax.z = mol[n].r.z;
	}
	VSSAdd (rCenter, 0.5, rMin, 0.5, rMax);
	DO_MOL VVSub (mol[n].r, rCenter);
	
	//anode
/*
	DO_MOL {
		if (mol[n].r.z + 0.5 * region.z  > 800) mol[n].anode = 1;
		else mol[n].anode = 0;
	}
	//cathode priority
	AllocMem (molNew, nMol, Mol);
	i = 0;
	DO_MOL {
		if (mol[n].anode == 0) {
			molNew[i] = mol[n];
			i ++;
		}
	}
	DO_MOL {
		if (mol[n].anode == 1) {
			molNew[i] = mol[n];
			i ++;
		}
	}
	DO_MOL mol[n] = molNew[n];

	free (molNew);
*/
}

void InitVels()
{
	int n;
	double vFac;

	VZero (vSum);
	DO_MOL {
		VRand3D (&mol[n].rv);
		vFac = velMag / sqrt (mol[n].mass);
		VScale (mol[n].rv, vFac);
		VVAdd (vSum, mol[n].rv);
	}
	DO_MOL {
		VVSAdd (mol[n].rv, -1. / nMol, vSum);
	}
	VZero (vSum);
	DO_MOL VVAdd (vSum, mol[n].rv);
}

void InitAccels()
{
	int n;

	DO_MOL VZero (mol[n].ra);
}
