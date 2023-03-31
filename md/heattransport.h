#ifndef HEATTRANSPORT_H
#define HEATTRANSPORT_H

#include "main.h"

#define NHIST 5

extern char nonq_mode[128];
extern double gravField;
extern VecI sizeHistGrid;
extern double **histGrid;
extern int countGrid, limitGrid, stepGrid;
extern double *profileT, *profileV;
extern double enTransSum, wallTempHi, wallTempLo;
extern Prop thermalCond;

void BuildNebrListHeat ();
void ComputeExternalForce ();
void ApplyBoundaryCondFlow ();
void ApplyBoundaryCondHeat ();
void GridAverage (int opCode);
void EvalProfile ();
void PrintProfile ();

#endif
