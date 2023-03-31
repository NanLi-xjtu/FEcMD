#ifndef NAMELIST_H
#define NAMELIST_H

#include "main.h"

typedef enum {N_I, N_R} VType;

typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;

//how to use NameList
/****************************
NameList nameList[] = {
	NameI (intVariable),
	NameR (doubleVector),
};
*****************************/

#define NameI(x) {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x) {#x, &x, N_R, sizeof (x) / sizeof (double)}
#define NP_I ((int *) (nameList[k].vPtr) + j)
#define NP_R ((double *) (nameList[k].vPtr) + j)

//global variables
extern Mol *mol;
extern VecR region, regionScale, vSum, pSum;
extern VecI initUcell;
extern Prop kinEnergy, pressure, totEnergy, tipHeight;
extern Prop KINENERGY_eV, PRESSURE_GPa, TOTENERGY_eV, TEMPAVG_K, TIPHEIGHT_A;
extern double deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum, latticePara, atomRadius, pedestal_thick;
extern double lUnit, mUnit, eUnit, tUnit, TUnit, PUnit, vUnit, fUnit, aUnit, DELTAT_fs, MASS_relative, DENSITY_g_cm3, TEMPERATURE_K, TIMENOW_ps, LATTIPARA_A, ATOMRADIUS_A, gap;
extern int moreCycles, nMol, stepAvg, stepMovie, stepCount, stepEquil, stepLimit;
extern int Nelems; //alloy
extern double *mass;
extern char **elems;

extern TempCell *tempCell;
extern int sizeTempGrid, stepTempGrid, stepSca;
extern double heatConduct, depositedHeat;

extern double minzInit;
extern char mdXYZFilename[128];

extern char boundaryCond[16], structure_type[16], interact_method[16], ensemble[16];

extern int Nthreads;

extern FILE *namelist;
extern FILE *summary;
extern FILE *movie;
extern FILE *velprofile;
extern FILE *temprofile;
extern FILE *tempcell;
extern FILE *veldist;
extern char mdPamsFilename[128], forceFilename[128], model[128];
extern int stepNebr;

int GetParaValue ();
int GetNameList (int argc, char **argv);
void PrintNameList (FILE *fp);
void PRINTSUMMARYTITLE ();
void PRINTSUMMARY ();
void PrintMovie();
void PrintFemocsin_xyz ();
void PrintOpen();
void PrintClose();
void GetIntVariable (char *line, char *name, int *variable);
void GetCharVariable (char *line, char *name, char *variable);
void GetDoubleVariable (char *line, char *name, double *variable);
int GetIntVariables (char *line, char *name, int *variable);
int GetCharVariables (char *line, char *name, char **variable);
int GetDoubleVariables (char *line, char *name, double *variable);
void GetIntVector (char *line, char *name, VecI *vector);
void GetDoubleVector (char *line, char *name, VecR *vector);

#endif
