#include "namelist.h"

Mol *mol;
VecR region, regionScale, vSum, pSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy, tipHeight;
Prop KINENERGY_eV, PRESSURE_GPa, TOTENERGY_eV, TEMPAVG_K, TIPHEIGHT_A;
double deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum, latticePara, atomRadius, pedestal_thick;
double lUnit, mUnit, eUnit, tUnit, TUnit, PUnit, vUnit, fUnit, aUnit, DELTAT_fs, MASS_relative, DENSITY_g_cm3, TEMPERATURE_K, TIMENOW_ps, LATTIPARA_A, ATOMRADIUS_A, gap;
int moreCycles, nMol, stepAvg, stepMovie, stepCount, stepEquil, stepLimit;
int Nelems; //alloy
double *mass;
char **elems;

TempCell *tempCell;
int sizeTempGrid, stepTempGrid, stepSca;
double heatConduct, depositedHeat;

char mdXYZFilename[128];

char boundaryCond[16], structure_type[16], interact_method[16], ensemble[16];

int Nthreads;

FILE *namelist;
FILE *summary;
FILE *movie;
FILE *velprofile;
FILE *temprofile;
FILE *tempcell;
FILE *femocsin_file;

char mdPamsFilename[128], forceFilename[128], model[128];

int summary_flag = 0, stepNebr = 0;

int GetParaValue ()
{
	FILE *input, *fp;
	char line[1024], filename[128], buff[80];
	int N, k;

	sprintf (buff, "%s", mdPamsFilename);
	if ((fp = fopen (buff, "r")) == 0) return (0);
	printf ("\nMD NameList -- variables\n");
	while (1) {
		if (fgets (line, 1024, fp) == NULL) break;
		GetIntVariable (line, "Nelems", &Nelems);
	}
	AllocMem2 (elems, Nelems, 10, char);
	AllocMem (mass, Nelems, double);
	for (k = 0; k < Nelems; k ++) mass[k] = 0.;
	rewind (fp);
	while (1) {
		if (fgets (line, 1024, fp) == NULL) break;

		GetCharVariable (line, "force_file", forceFilename);
		GetCharVariable (line, "force_type", force_type);
		GetCharVariable (line, "tip_file", tipFilename);
		GetCharVariable (line, "boundary", boundaryCond);
		GetCharVariable (line, "cell_file", cellFilename);
		GetCharVariable (line, "cell_order", cellOrder);
		GetCharVariable (line, "structure_type", structure_type);
		GetCharVariable (line, "interact_method", interact_method);
		GetCharVariable (line, "ensemble", ensemble);
		N = GetCharVariables (line, "elems", elems);
		N = GetDoubleVariables (line, "mass", mass);
		if (N != 0) {
			sprintf (buff, "error(namelist.c): the number of elements %d %d\n", N, Nelems);
			if (N != Nelems) Error (buff);
		}
		GetDoubleVariable (line, "deltaT", &DELTAT_fs);
		GetDoubleVariable (line, "temperature", &TEMPERATURE_K);
		GetDoubleVariable (line, "Trate", &Trate);
		GetDoubleVariable (line, "vMax_set", &vMax_set);
		GetDoubleVariable (line, "vMin_set", &vMin_set);
		GetIntVariable (line, "stepAvg", &stepAvg);
		GetIntVariable (line, "stepMovie", &stepMovie);
		GetIntVariable (line, "stepEquil", &stepEquil);
		GetIntVariable (line, "stepLimit", &stepLimit);
		GetIntVariable (line, "randSeed", &randSeed);
		GetDoubleVector (line, "regionScale", &regionScale);
		GetDoubleVariable (line, "lattice", &LATTIPARA_A);
		GetIntVariable (line, "Nthreads", &Nthreads);
		GetDoubleVariable (line, "pedestal_thick", &pedestal_thick);
		GetIntVector (line, "initUcell", &initUcell);
		GetDoubleVariable (line, "tip_R", &tip_r);
		GetDoubleVariable (line, "tip_r", &tip_r1);
		GetDoubleVariable (line, "tip_h", &tip_h);
		GetDoubleVariable (line, "tip_theta", &tip_theta);
		GetIntVariable (line, "nebrTabFac", &nebrTabFac);
		GetDoubleVariable (line, "rNebrShell", &rNebrShell);
		GetIntVariable (line, "stepAdjustTemp", &stepAdjustTemp);
		GetIntVariable (line, "stepInitlzTemp", &stepInitlzTemp);
		GetDoubleVariable (line, "rSurfCut", &rSurfCut);
		GetDoubleVariable (line, "nSurfCut", &nSurfCut);
		GetDoubleVariable (line, "rFlyCut", &rFlyCut);
		GetDoubleVariable (line, "nFlyCut", &nFlyCut);
		GetIntVariable (line, "stepSurf", &stepSurf);
		GetIntVariable (line, "stepLigancy", &stepLigancy);
		GetDoubleVariable (line, "rangeVel", &rangeVel);
		GetIntVariable (line, "limitVel", &limitVel);
		GetIntVariable (line, "stepVel", &stepVel);
		GetIntVariable (line, "sizeHistVel", &sizeHistVel);
		GetDoubleVariable (line, "extPressure", &extPressure);
		GetDoubleVariable (line, "massS", &massS);
		GetDoubleVariable (line, "massV", &massV);
		GetDoubleVariable (line, "Maxwell_rate", &Maxwell_rate);
		GetDoubleVariable (line, "Maxwell_max", &Maxwell_max);
		GetDoubleVariable (line, "Maxwell_begin", &Maxwell_begin);
		GetDoubleVariable (line, "gravField", &gravField);
		GetDoubleVariable (line, "wallTempHi", &wallTempHi);
		GetDoubleVariable (line, "wallTempLo", &wallTempLo);
		GetIntVariable (line, "limitGrid", &limitGrid);
		GetIntVariable (line, "stepGrid", &stepGrid);
		GetIntVector (line, "sizeHistGrid", &sizeHistGrid);
		GetIntVariable (line, "sizeTempGrid", &sizeTempGrid);
		GetIntVariable (line, "stepTempGrid", &stepTempGrid);
		GetDoubleVariable (line, "heatConduct", &heatConduct);
		GetDoubleVariable (line, "depositedHeat", &depositedHeat);
		GetDoubleVariable (line, "lj_strength", &lj_e);
		GetDoubleVariable (line, "lj_length", &lj_r);
	}
	fclose (fp);

	if (strcmp(force_type, "metal") == 0 || strcmp(force_type, "alloy") == 0 || \
	    strcmp(force_type, "eamfs") == 0) strcpy (eamFilename, forceFilename);
	else if (strcmp(force_type, "snap") == 0) strcpy (snapFilename, forceFilename);

	sprintf (filename, "out/md/namelist.dat");
	namelist = WriteFile (filename);
	fclose (namelist);

	return (0);
}

NameList nameList[] = {
	//MD
	NameR (DELTAT_fs),//--fs
	NameR (TEMPERATURE_K),//--K
	NameR (Trate),//--K/step
	NameI (Nelems),
	NameI (initUcell),
	NameI (stepAvg),
	NameI (stepMovie),
	NameI (stepEquil),
	NameI (stepLimit),
	NameI (randSeed),
	NameI (Nthreads),
	NameR (LATTIPARA_A),
	NameR (regionScale),
	NameI (boundary_fly),
	NameR (vMax_set),
	NameR (vMin_set),

	//neigbor list
	NameI (nebrTabFac), //determines how much storage should be provided for the neighbor list (per atom)
	NameR (rNebrShell), //corresponding to delta r

	//adjust temperature
	NameI (stepAdjustTemp),
	NameI (stepInitlzTemp),

	//g : the external field driving the flow
	NameR (gravField),	

	//constant-temperature walls
	NameR (wallTempHi),
	NameR (wallTempLo),

	//gird computation
	NameI (limitGrid),
	NameI (sizeHistGrid),
	NameI (stepGrid),

	//heat diffusion grids
	NameI (sizeTempGrid),
	NameI (stepTempGrid),
	NameR (heatConduct),
	NameR (depositedHeat),

	//RDF
	NameI (limitRdf),
	NameR (rangeRdf),
	NameI (sizeHistRdf),
	NameI (stepRdf),

	//surface atoms' neighbors number and neighbor cut radius
	NameR (rSurfCut),
	NameR (nSurfCut),
	NameR (rFlyCut),
	NameR (nFlyCut),
	NameI (stepSurf),
	NameI (stepLigancy),

	//velocity distriburion
	NameI (limitVel),
	NameR (rangeVel),
	NameI (sizeHistVel),
	NameI (stepVel),

	//tip maker
	NameR (tip_r),
	NameR (tip_r1),
	NameR (tip_h),
	NameR (tip_theta),

	//fix pedestal
	NameR (pedestal_thick),

	//NPT
	NameR (extPressure),
	NameR (massS),
	NameR (massV),

	//Maxwell stress
	NameR (Maxwell_rate),
	NameR (Maxwell_max),
	NameR (Maxwell_begin),
};

int GetNameList (int argc, char **argv)
{
	int id, j, k, match, ok, length, N;
	char buff[80], *token, line[1024], filename[128];
	FILE *fp;

	sprintf (buff, "%s", mdPamsFilename);
	if ((fp = fopen (buff, "r")) == 0) return (0);
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++)
		nameList[k].vStatus = 0;
	ok = 1;
	while (1){
		fgets (buff, 80, fp);
		if (feof (fp)) break;
		token = strtok (buff, "\t\n");
		if (! token) break;
		match = 0;
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k++){
			if (strcmp (token, nameList[k].vName) == 0){
				match = 1;
				if (nameList[k].vStatus == 0){
					nameList[k].vStatus = 1;
					for (j = 0; j < nameList[k].vLen; j++){
						token = strtok (NULL, ", \t\n");
						if (token){
							switch (nameList[k].vType){
								case N_I:
									*NP_I = atol (token);
									break;
								case N_R:
									*NP_R = atof (token);
									break;
							}
						} else{
							nameList[k].vStatus = 2;
							ok = 0;
						}
					}
					token = strtok (NULL, ", \t\n");
					if (token){
						nameList[k].vStatus = 3;
						ok = 0;
					}
					break;
				}else{
					nameList[k].vStatus = 4;
					ok = 0;
				}
			}
		}
		if (! match) ok = 0;
	}
	fclose (fp);
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++){
		if (nameList[k].vStatus != 1) ok = 0;
	}

	sprintf (buff, "%s", mdPamsFilename);
	if ((fp = fopen (buff, "r")) == 0) return (0);
	printf ("\nMD NameList -- variables\n");
	AllocMem2 (elems, Nelems, 10, char);
	AllocMem (mass, Nelems, double);
	for (k = 0; k < Nelems; k ++) mass[k] = 0.;
	while (1) {
		if (fgets (line, 1024, fp) == NULL) break;

		GetCharVariable (line, "force_file", forceFilename);
		GetCharVariable (line, "force_type", force_type);
		GetCharVariable (line, "tip_file", tipFilename);
		GetCharVariable (line, "boundary", boundaryCond);
		GetCharVariable (line, "cell_file", cellFilename);
		GetCharVariable (line, "cell_order", cellOrder);
		GetCharVariable (line, "structure_type", structure_type);
		GetCharVariable (line, "interact_method", interact_method);
		GetCharVariable (line, "ensemble", ensemble);
		N = GetCharVariables (line, "elems", elems);
		N = GetDoubleVariables (line, "mass", mass);
		if (N != 0) {
			sprintf (buff, "error(namelist.c): the number of elements %d %d\n", N, Nelems);
			if (N != Nelems) Error (buff);
		}
	}
	fclose (fp);

	if (strcmp(force_type, "metal") == 0 || strcmp(force_type, "alloy") == 0 || \
	    strcmp(force_type, "eamfs") == 0) strcpy (eamFilename, forceFilename);
	else if (strcmp(force_type, "snap") == 0) strcpy (snapFilename, forceFilename);
	PrintNameList (stdout);
	sprintf (filename, "out/md/namelist.dat");
	namelist = WriteFile (filename);
	PrintNameList (namelist);
	fclose (namelist);

	return (ok);
}

void PrintNameList (FILE *fp)
{
	int j, k;

	fprintf (fp, "\nMD NameList -- data\n");
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k ++){
		fprintf (fp, "%s\t", nameList[k].vName);
		if (strlen (nameList[k].vName) < 8) fprintf (fp, "\t");
		if (nameList[k].vStatus > 0){
			for (j = 0; j < nameList[k].vLen; j ++){
				switch (nameList[k].vType){
					case N_I:
						fprintf (fp, "%d ", *NP_I);
						break;
					case N_R:
						fprintf (fp, "%#g ", *NP_R);
						break;
				}
			}
		}
		switch (nameList[k].vStatus){
			case 0:
				fprintf (fp, "** no data");
				break;
			case 1:
				break;
			case 2:
				fprintf (fp, "** missing data");
				break;
			case 3:
				fprintf (fp, "** extra data");
				break;
			case 4:
				fprintf (fp, "** multiply defined");
				break;
		}
		fprintf (fp, "\n");
	}
	fprintf (fp, "----\n");
}

void PRINTSUMMARYTITLE ()
{
	printf ("\nsummary -- data    \n");
	printf ("step    time(ps)  <p>(M m/s)  <E>(eV)   <Ek>(eV)  <T>(K)    <P>(GPa)    h(A)    stepNebr NatomFly ");
	if (strcmp(nonq_mode, "heattransport") == 0) printf ("conductivity(W/m k) ");
	printf ("\n");

	if ((summary = fopen ("out/md/summary.dat", "a+")) == NULL){
		printf ("\nopen summary file error");
		getchar ();
		exit (1);
	}
	if (summary_flag == 0) {
		fprintf(summary, "step    time(ps)   <v>(m/s)  <E>(eV)   stddva    <Ek>(eV)  stddva    <T>(K)    stddva    <P>(GPa)  stddva    h(A)      stddva  NatomFly ");
		if (strcmp(nonq_mode, "heattransport") == 0) fprintf(summary, "conductivity(W/m k) ");
		fprintf(summary, "\n");
	}

	fclose (summary);
}

void PRINTSUMMARY ()
{
	printf("%6d %9.4f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6d %10d ",
		stepCount, TIMENOW_ps, VCSum (pSum) / nMol * vUnit * mUnit / M_C12, TOTENERGY_eV.sum,
		KINENERGY_eV.sum, TEMPAVG_K.sum, PRESSURE_GPa.sum, TIPHEIGHT_A.sum, stepNebr, atomFlyCount);
	if (strcmp(nonq_mode, "heattransport") == 0) printf ("%16f ", thermalCond.sum * (eUnit / (tUnit * lUnit * TUnit)));
	printf ("\n");

	if ((summary = fopen ("out/md/summary.dat", "a+")) == NULL){
		printf ("\nopen summary file error");
		getchar ();
		exit (1);
	}
	fprintf(summary,
		"%6d %9.4f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %d ",
		stepCount, TIMENOW_ps, VCSum (vSum) / nMol * vUnit, PropEst (TOTENERGY_eV),
		PropEst (KINENERGY_eV), PropEst (TEMPAVG_K), PropEst (PRESSURE_GPa), PropEst (TIPHEIGHT_A), atomFlyCount);
	if (strcmp(nonq_mode, "heattransport") == 0) fprintf (summary, "%e ", thermalCond.sum * (eUnit / (tUnit * lUnit * TUnit)));
	fprintf (summary, "\n");

	summary_flag = 1;

	fclose (summary);
}

void PrintMovie ()
{
	if ((movie = fopen ("out/md/md.movie", "a+")) == NULL){
		printf ("\nopen movie file error");
		getchar ();
		exit (1);
	}

	fprintf (movie, "%d\nLattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" SolutionReader properties=id:I:1:species:S:1:pos:R:3:ligancy:I:1:flag:I:1:surf_F(V/A):R:1:f(V/A):R:3:|f|(V/A):R:1:a(m/s2):R:1:Ek(eV):R:1:T(K):R:1\n", nMol, region.x * (lUnit * 1.e10), region.y * (lUnit * 1.e10), region.z * (lUnit * 1.e10));

	int n;

	DO_MOL fprintf (movie, "%d %s %e %e %e %d %d %e %e %e %e %e %e %e %e\n", \
			n, mol[n].elem, (mol[n].r.x + 0.5 * region.x) * (lUnit * 1.e10), (mol[n].r.y + 0.5 * region.y) * (lUnit * 1.e10), \
			(mol[n].r.z + 0.5 * region.z) * (lUnit * 1.e10), mol[n].ligancy, mol[n].flag, \
			VLen(mol[n].F) * (eUnit / eleChar) / (lUnit * 1.e10), \
			mol[n].f.x * (eUnit / eleChar) / (lUnit * 1.e10), \
			mol[n].f.y * (eUnit / eleChar) / (lUnit * 1.e10), \
			mol[n].f.z * (eUnit / eleChar) / (lUnit * 1.e10), \
			VLen(mol[n].f) * (eUnit / eleChar) / (lUnit * 1.e10), \
			VLen(mol[n].ra) * aUnit, \
			0.5 * mol[n].mass * (VLenSq(mol[n].rv)) * (eUnit / eleChar), \
			1. / NDIM * mol[n].mass * (VLenSq(mol[n].rv)) * (TUnit));

	fclose (movie);
}

void PrintFemocsin_xyz ()
{
	int n, N;
	char filename[128];
	double xmid, ymid, xmax, xmin, ymax, ymin;

	xmin = ymin = 1.e100;
	xmax = ymax = -1.e100;
	DO_MOL {
		if (fabs(mol[n].r.z - minzInit) < 0.1) {
			if (xmin > mol[n].r.x) xmin = mol[n].r.x;
			if (xmax < mol[n].r.x) xmax = mol[n].r.x;
			if (ymin > mol[n].r.y) ymin = mol[n].r.y;
			if (ymax < mol[n].r.y) ymax = mol[n].r.y;
		}
	}
	xmid = (xmin + xmax) / 2.;
	ymid = (ymin + ymax) / 2.;

	sprintf (filename, "out/md/femocs.in.xyz");
	if ((femocsin_file = fopen ("out/md/femocs.in.xyz", "w")) == NULL){
		printf ("\nopen out/md/femocs.in.xyz file error");
		getchar ();
		exit (1);
	}
	N = 0;
	DO_MOL {
		if (mol[n].anode == 0) N ++;
	}
	fprintf (femocsin_file, "%d\n", N);
	fprintf (femocsin_file, "%f 0. 0. 0. %f 0. 0. 0. %f ", region.x  * lUnit * 1.e10, \
		 region.y * lUnit * 1.e10, region.z * lUnit * 1.e10); //--A
	fprintf (femocsin_file, "Properties=species:S:1:pos:R:3:structure_type:I:1\n");
	DO_MOL {
		if (mol[n].anode == 0) {
			fprintf (femocsin_file, "%s %f %f %f %d\n", mol[n].elem, (mol[n].r.x - xmid) * lUnit * 1.e10, \
				 (mol[n].r.y - ymid) * lUnit * 1.e10, (mol[n].r.z - minzInit) * lUnit * 1.e10, mol[n].flag); //--A
		}
	}
	fclose (femocsin_file);
}

void PrintOpen()
{
	char filename[128];

//	printf ("=== print open...\n");
	system ("rm -rf out/md");
	system ("mkdir out/md");
	
	if ((tempcell = fopen ("out/md/tempcell.dat", "w")) == NULL){
		printf ("\nopen tempcell file error");
		getchar ();
		exit (1);
	}
	if ((rdf = fopen ("out/md/rdf.dat", "w")) == NULL){
		printf ("\nopen rdf file error");
		getchar ();
		exit (1);
	}
	if ((ligancy = fopen ("out/md/ligancy.dat", "w")) == NULL){
		printf ("\nopen ligancy file error");
		getchar ();
		exit (1);
	}

//	printf("PrintOpen -- ok\n");
}

void PrintClose()
{
	fclose (tempcell);
	fclose (rdf);
	fclose (ligancy);
}

void GetIntVariable (char *line, char *name, int *variable)
{
	int size;
	char *words;

	if (*variable != 0) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);

		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		*variable = (int)(atof(words));
	}

	if (*variable != 0) printf ("%s = %d\n", name, *variable);
}

void GetCharVariable (char *line, char *name, char *variable)
{
	int size;
	char *words;

	if (variable[0] >= 65 && variable[0] <= 122) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);

		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		sscanf (words, "%s", variable);
	}

	if (variable[0] >= 65 && variable[0] <= 122) printf ("%s = %s\n", name, variable);
}

void GetDoubleVariable (char *line, char *name, double *variable) //the initial value of variable must be zero
{
	int size;
	char *words;

	if (fabs(*variable) > 1.e-60) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);
		
		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		*variable = atof(words);
	}

	if (fabs(*variable) > 1.e-60) printf ("%s = %f\n", name, *variable);
}

int GetCharVariables (char *line, char *name, char **variable)
{
	int size, i, N;
	char *words;

	if (variable[0][0] >= 65 && variable[0][0] <= 122) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL || words[0] == 35) break;
			strcpy (variable[N], words);
			N ++;
		}
	}
	if (variable[0][0] >= 65 && variable[0][0] <= 122) {
		printf ("%s =", name);
		for (i = 0; i < N; i ++) printf (" %s", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}

int GetDoubleVariables (char *line, char *name, double *variable) //the initial value of variable must be zero
{
	int size, i, N;
	char *words;

	if (fabs(variable[0]) > 1.e-60) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL || words[0] == 35) break;
			variable[N] = atof (words);
			N ++;
		}
	}
	if (fabs(variable[0]) > 1.e-60) {
		printf ("%s =", name);
		for (i = 0; i < N; i ++) printf (" %.3f", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}

int GetIntVariables (char *line, char *name, int *variable) //the initial value of variable must be zero
{
	int size, i, N;
	char *words;

	if (fabs(variable[0]) != 0) return (0);
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		N = 0;
		while (1) {
			words = strtok(NULL,"' \t\n\r\f");
			if (words == NULL || words[0] == 35) break;
			variable[N] = atoi (words);
			N ++;
		}
	}
	if (fabs(variable[0]) != 0) {
		printf ("%s =", name);
		for (i = 0; i < N; i ++) printf (" %d", variable[i]);
		printf ("\n");
		return (N);
	}

	return (0);
}

void GetIntVector (char *line, char *name, VecI *vector)
{
	int size;
	char *words;

	if (VProd(*vector) != 0) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);

		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		vector->x = (int)(atof(words));
		words = strtok(NULL,"' \t\n\r\f");
		vector->y = (int)(atof(words));
		words = strtok(NULL,"' \t\n\r\f");
		vector->z = (int)(atof(words));
	}

	if (VProd(*vector) != 0) printf ("%s = %d %d %d\n", name, vector->x, vector->y, vector->z);
}

void GetDoubleVector (char *line, char *name, VecR *vector)
{
	int size;
	char *words;

	if (fabs(VProd(*vector)) > 1.e-60) return;
	size = strlen (name);
	if (strncasecmp (line, name, size) == 0) {
		strcpy (line, line + size);

		words = strtok(line,"' \t\n\r\f");
		words = strtok(NULL,"' \t\n\r\f");
		vector->x = atof(words);
		words = strtok(NULL,"' \t\n\r\f");
		vector->y = atof(words);
		words = strtok(NULL,"' \t\n\r\f");
		vector->z = atof(words);
	}

	if (fabs(VProd(*vector)) > 1.e-60) printf ("%s = %f %f %f\n", name, vector->x, vector->y, vector->z);
}
