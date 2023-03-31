#ifndef INTERFACE_H
#define INTERFACE_H

#include "main.h"

typedef struct{
	double T, Eksum, vFac;
	int N, Natoms;
} TCELL;

#define DO_TCELL for (i = 0; i < Ncell.x; i ++) { \
		for (j = 0; j < Ncell.y; j ++) { \
			for (k = 0; k < Ncell.z; k ++)

#define END_TCELL }}

extern int force_flag, T_flag;
extern FILE *force_file, *T_file;
extern TCELL ***Tcell;
extern VecI Ncell;
extern char force_sort[128];
extern double Maxwell_rate, Maxwell_max, Maxwell_stress, Maxwell_begin;

void ForceInterface ();
void TemperatureInterface (int step);
void PrintTcellMovie ();
void MaxwellStress (double theta); //--GPa;

#endif
