#ifndef BOUNDARY
#define BOUNDARY

#include "main.h"

#define VWrap(v, t)					\
	if(v.t >= 0.5 * region.t) v.t -= region.t;	\
	else if (v.t < -0.5 * region.t) v.t += region.t
#define VWrapAll(v)					\
	{VWrap (v, x);					\
	 VWrap (v, y);					\
	 VWrap (v, z);}

#define VWrap_NPT(v, t)					\
	if(v.t >= 0.5) v.t -= 1.;			\
	else if (v.t < -0.5) v.t += 1.
#define VWrapAll_NPT(v)					\
	{VWrap_NPT (v, x);					\
	 VWrap_NPT (v, y);					\
	 VWrap_NPT (v, z);}

extern int atomFlyCount;
extern int atomFlyFlag, boundary_fly;
extern double minzInit, maxzInit;

void ApplyBoundaryCond ();
void AtomFlyBox ();
void AtomFlySurface ();
void DeleteAtoms ();
void FindPedestal ();

#endif
