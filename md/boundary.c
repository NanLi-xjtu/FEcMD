#include "boundary.h"

int atomFlyCount = 0;
int atomFlyFlag = 0, boundary_fly = 0;
double minzInit, maxzInit;

void ApplyBoundaryCond ()
{
	int n;

	if (strcmp(ensemble, "NPT") == 0) {
		#pragma omp parallel for num_threads(Nthreads)
		DO_MOL VWrapAll_NPT(mol[n].r);
	} else {
		#pragma omp parallel for num_threads(Nthreads)
		DO_MOL VWrapAll(mol[n].r)
	}
}

void AtomFlyBox ()
{
	int n, flag = 0;

	//delete the atom flying out of the box
	#pragma omp parallel for num_threads(Nthreads)
	DO_MOL {
		if (mol[n].r.x >= 0.5 * region.x || mol[n].r.x < - 0.5 * region.x || \
		    mol[n].r.y >= 0.5 * region.y || mol[n].r.y < - 0.5 * region.y || \
		    mol[n].r.z >= 0.5 * region.z || mol[n].r.z < - 0.5 * region.z) {
			mol[n].deleteFlag = 1;
			flag = 1;
		}
	}
	if (flag == 1) DeleteAtoms ();
}

void AtomFlySurface ()
{
	int n, i;

	//delete the atom flying out of the surface
	i = 0;
	DO_MOL {
		if (mol[n].flag == 3) {
			atomFlyCount ++;
			atomFlyFlag = 1;
			continue;
		}
		mol[i] = mol[n];
		i ++;
	}
	nMol = i;
}

void DeleteAtoms ()
{
	int n, i;

	i = 0;
	DO_MOL {
		if (mol[n].deleteFlag == 0) {
			mol[i] = mol[n];
			i ++;
		}
	}
	atomFlyCount += nMol - i;
	nMol = i;

	
	nebrListFlag = 0;
	atomFlyFlag = 1;
}

void FindPedestal ()
{
	int n, i;

	minzInit = 1.e100;
	maxzInit = -1.e100;
	DO_MOL{
		minzInit = Min(minzInit, mol[n].r.z);
		maxzInit = Max(maxzInit, mol[n].r.z);
	}
	
	i = 0;
	DO_MOL {
		if (fabs (mol[n].r.z - minzInit) < pedestal_thick) mol[n].pedest = 1;
		else if (fabs (mol[n].r.z - maxzInit) < pedestal_thick && \
			 mol[n].anode == 1) {
			mol[n].pedest = 1;
		} else mol[n].pedest = 0;
	}

}
