#ifndef FORCES_H
#define FORCES_H

#include "main.h"

#define N_OFFSET 14  
#define OFFSET_VALS						\
	{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1},	\
	{1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, 		\
	{-1,-1,1}, {0,-1,1}, {1,-1,1}}
#define VCellWrap(t)						\
	if (m2v.t >= cells.t){					\
		m2v.t = 0;					\
		shift.t = region.t;				\
	}else if (m2v.t < 0){					\
		m2v.t = cells.t - 1;				\
		shift.t = - region.t;				\
	}
#define VCellWrapAll()						\
	{VCellWrap (x);						\
	 VCellWrap (y);						\
	 VCellWrap (z);}
#define DO_CELL(j, m)						\
	for (j = cellList[m]; j >= 0; j = cellList[j])

extern VecI cells;
extern double dispHi, rNebrShell, lj_r, lj_e;
extern int *cellList, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax;

extern double embedWt, rSwitch, splineA2, splineA3;

void ComputeForces ();
void ComputeForcesCell ();
void ComputeForcesNebr ();
void ComputeForcesEam ();
void EvalEamParams ();

#endif
