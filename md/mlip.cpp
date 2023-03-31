#include "mlip.h"

int mlip_init_flag = 0;
AtomLMP atom[1];
NeighList list[1];

PairMLIP mlp;

PairMLIP::PairMLIP()
{
}

PairMLIP::~PairMLIP()
{
	MLIP_finalize();
}

void PairMLIP::compute(int eflag, int vflag)
{
	int n, i;
	double lattice[9];
	int j, k;

	list->inum = nMol;

	AllocMem (list->ilist, nMol, int);
	AllocMem (list->numneigh, nMol, int);
	AllocMem2 (list->firstneigh, nMol, nebrTabFac, int);
	AllocMem (atom->type, nMol, int);
	AllocMem2 (atom->x, nMol, 3, double);
	AllocMem2 (atom->v, nMol, 3, double);
	AllocMem2 (atom->f, nMol, 3, double);

	lattice[0] = region.x * lUnit * 1.e10; //--A
	lattice[1] = 0.0;
	lattice[2] = 0.0;
	lattice[3] = 0.0;
	lattice[4] = region.y * lUnit * 1.e10;
	lattice[5] = 0.0;
	lattice[6] = 0.0;
	lattice[7] = 0.0;
	lattice[8] = region.z * lUnit * 1.e10;

	DO_MOL {
		atom->type[n] = -1;
		for (i = 0; i < Nelems; i ++) {
			if (strcmp(mol[n].elem, elems[i]) == 0) atom->type[n] = i + 1;
		}
//		printf ("%s, %d\n", mol[n].elem, atom->type[n]);
//		atom->type[n] = 2;
		atom->x[n][0] = mol[n].r.x * lUnit * 1.e10; //--A
		atom->x[n][1] = mol[n].r.y * lUnit * 1.e10;
		atom->x[n][2] = mol[n].r.z * lUnit * 1.e10;
		list->ilist[n] = n;
		list->numneigh[n] = mol[n].Nnebr;
		for (i = 0; i < mol[n].Nnebr; i ++) list->firstneigh[n][i] = mol[n].id_nebr[i];
	}

	if (mode == 0) { // nbh version
		double energy = 0;
		double *p_site_en = NULL;
		double **p_site_virial = NULL;
		
		DO_MOL atom->f[n][0] = atom->f[n][1] = atom->f[n][2] = 0.;
		energy = 0.;
		MLIP_calc_nbh(boundaryCond,
			lattice,
			list->inum, 
			list->ilist, 
			list->numneigh, 
			list->firstneigh,
			atom->nlocal,
			atom->nghost,
			atom->x,  //--A
			atom->type,
			atom->f,  //--V/A
			energy,
			p_site_en,      // if NULL no site energy is calculated
			p_site_virial); // if NULL no virial stress per atom is calculated

//		DO_MOL printf ("%d %f %f %f %f %f %f %f %d\n", n, atom->x[n][0], atom->x[n][1], atom->x[n][2], \
				atom->f[n][0], atom->f[n][1], atom->f[n][2], sqrt(Sqr(atom->f[n][0]) + Sqr(atom->f[n][1]) + Sqr(atom->f[n][2])), \
				list->numneigh[n]);
//		for (i = 0; i < list->numneigh[1]; i ++) printf ("%d\n", list->firstneigh[1][i]);
//		printf ("energy = %f eV\n", energy);
//		if (stepCount == 10) exit (1);
		
		virSum = 0.;
		DO_MOL {
			VSet (mol[n].f, atom->f[n][0] / ((eUnit / eleChar) / (lUnit * 1.e10)), \
					atom->f[n][1] / ((eUnit / eleChar) / (lUnit * 1.e10)), \
					atom->f[n][2] / ((eUnit / eleChar) / (lUnit * 1.e10)));  //dimensionless
		}
		uSum = energy; //--eV
		uSum = uSum / (eUnit / eleChar);//dimensionless
	} else { //cfg version
/*
	    double en = 0.0;
	    double virstr[9];

	    MLIP_calc_cfg(list->inum, lattice, atom->x, atom->type, atom->tag, en, atom->f, virstr);
*/
	}

	free (list->ilist);
	free (list->numneigh);
	free2 (list->firstneigh);
	free (atom->type);
	free2 (atom->x);
	free2 (atom->v);
	free2 (atom->f);
}

void PairMLIP::allocate()
{

}

void PairMLIP::settings()
{
	FILE *input;
	char line[1024], *words;

	sprintf (MLIPsettings_filename, "in/mlip.ini");
	sprintf (MLIPlog_filename, "out/mlip.log");

	input = ReadFile (forceFilename);
	cutoff = 0.;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;
		words = strtok(line,"' \t\n\r\f");
		if (strcmp (words, "max_dist") == 0) {
			words = strtok(NULL,"' \t\n\r\f");
			words = strtok(NULL,"' \t\n\r\f");
			cutoff = atof (words);
		}
	}
	fclose (input);

	atom->nlocal = atom->nghost = 0;
	atom->ntypes = Nelems;
	mode = 0;
}

void PairMLIP::coeff()
{
	allocate();
}

void PairMLIP::init_style()
{
	int i, n;

	if (mlip_init_flag == 0) {
		printf ("\ninit ...\n");
		printf ("ntypes = %d, cutoff = %.3f A, mode = %d\n", atom->ntypes, cutoff, mode);
		mlip_init_flag = 1;
	}

	MLIP_init(MLIPsettings_filename, MLIPlog_filename, atom->ntypes, cutoff, mode);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIP::init_one(int i, int j)
{
  return cutoff;
}

void ComputeForcesMTP ()
{
	mlp.compute (0, 0);
}
