#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "main.h"

typedef struct{
	double nebrnu;
	int surface;
} SurfAtoms;

#define SurfZero(v)		\
	{v.nebrnu = 0;		\
	 v.surface = 0;}

extern SurfAtoms *surf;
extern double rRdf, nSurfCut, rSurfCut, rFlyCut, nFlyCut;
extern int countSurf, limitSurf, stepSurf, stepLigancy;

extern double *histRdf, rangeRdf;
extern int countRdf, limitRdf, sizeHistRdf, stepRdf;
extern FILE *rdf, *ligancy;

extern int ligancy_flag;

void ApplyBoundaryCondSurface ();
void EvalRdf ();
void EvalLigancy ();
void SurfaceAtoms ();

#endif
