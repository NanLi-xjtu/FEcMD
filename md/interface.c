#include "interface.h"

int force_flag = 0, T_flag = 0, step_heat = 0, femocsin_flag = 0;
FILE *force_file;
TCELL ***Tcell;
VecI Ncell;
double heat_dt, md_timestep;
char force_sort[128], Tmode[128];
double Maxwell_rate, Maxwell_max, Maxwell_stress = 0., Maxwell_begin;

void ForceInterface ()
{
	char filename[128], line[1024];
	int N = 0; //the number of surface atoms
	int n, id;
	double temp, Fnorm;
	VecR F;

	printf ("\n=== force interface...\n");
	sprintf (filename, "out/forces.xyz");
	force_file = ReadFile (filename);
	fgets (line, 1024, force_file);
	sscanf (line, "%d", &N);
	printf ("the number of surface atom is %d\n", N);
	if (N == 0) {
		fclose (force_file);
		return;
	}
	fgets (line, 1024, force_file);
	DO_MOL VZero (mol[n].F);
	for (n = 0; n < N; n ++) {
		if (fgets (line, 1024, force_file) == 0) break;
		sscanf (line, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &id, &temp, \
			&temp, &temp, &temp, &F.x, &F.y, &F.z, &Fnorm, &temp, &temp); //--V/A
		VScale(F, lUnit * 1.e10 / (eUnit / eleChar)); //--dimensionless
		if (isnormal(Fnorm) != 0 && fabs(Fnorm) > 1.e-100) VCopy (mol[id].F, F); //--dimensionless
	}

	fclose (force_file);
}

void TemperatureInterface (int step)
{
	TCELL *cell;
	char filename[128], line[1024];
	int N = 0; //the number of T points from femocs
	int n, id, i, j, k;
	double temp, *T, *T_electron, *heat_tot, *heat_nottingham, *heat_joule, coff;
	double heat_tot_max, heat_nottingham_max, heat_joule_max, Tmax, T1_3avg, Tele_max, Tele1_3avg;
	VecR *r, rmax, rmin, threshold, box;
	FILE *input;
	FILE *output;

	printf ("\n=== temperature interface...\n");
	//read the size of cell
	if (femocsin_flag == 0) {
		sprintf (filename, "in/femocs.in");
		input = ReadFile (filename);
		while (1) {
			if (fgets (line, 1024, input) == 0) break;
			if (strncasecmp (line, "Ncell", 5) == 0) sscanf (line, "Ncell = %d %d %d", &Ncell.x, &Ncell.y, &Ncell.z);
			GetDoubleVariable (line, "heat_dt =", &heat_dt); //--fs
			GetDoubleVariable (line, "md_timestep =", &md_timestep); //--fs
                        GetCharVariable (line, "temperature_mode =", Tmode);
		}
		fclose (input);
		femocsin_flag = 1;
	}
	if ((double)step * md_timestep < (double)step_heat * heat_dt) return;
	else step_heat ++;
	//read teamperature points from femocs
	if (strcmp (Tmode, "single") == 0) sprintf (filename, "out/ch_solver.xyz");
	else if (strcmp (Tmode, "double") == 0) sprintf (filename, "out/temperature_phonon.xyz");
	if ((input = fopen (filename, "r")) == NULL) {
		printf ("%s does not exsit\n", filename);
		return;
	}
	fgets (line, 1024, input);
	sscanf (line, "%d", &N);
	printf ("the number of temperature points is %d\n", N);
	fgets (line, 1024, input);

	//allocate memory
	AllocMem3 (Tcell, Ncell.x, Ncell.y, Ncell.z, TCELL);
	AllocMem (T, N, double);
        AllocMem (T_electron, N, double);
	AllocMem (heat_tot, N, double);
	AllocMem (heat_nottingham, N, double);
	AllocMem (heat_joule, N, double);
	AllocMem (r, N, VecR);
	DO_TCELL {
		Tcell[i][j][k].N = 0;
		Tcell[i][j][k].Natoms = 0;
		Tcell[i][j][k].T = 0.;
		Tcell[i][j][k].Eksum = 0.;
	} END_TCELL

	//read heat diffusion temperature from femocs
	for (n = 0; n < N; n ++) {
		if (fgets (line, 1024, input) == 0) break;
		if (strcmp (Tmode, "single") == 0) {
			sscanf (line, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &id, &r[n].x, &r[n].y, &r[n].z, \
				&temp, &temp, &temp, &heat_tot[n], &heat_nottingham[n], &heat_joule[n], &temp, &T[n], &temp, &temp);
			T_electron[n] = 0.;
		} else if (strcmp (Tmode, "double") == 0) {
			sscanf (line, "%d %lg %lg %lg %lg %lg %lg %lg %lg", &id, &r[n].x, &r[n].y, &r[n].z, \
				&heat_tot[n], &heat_nottingham[n], &heat_joule[n], &T_electron[n], &T[n]);
		}
	} //--A, K
	fclose (input);

	//find rmax and rmin for temperature cell
	VSet (rmax, -1.e100, -1.e100, -1.e100);
	VSet (rmin, 1.e100, 1.e100, 1.e100);
	for (n = 0; n < N; n ++) {
		if (r[n].z < 0) continue;
		rmax.x = Max (rmax.x, r[n].x); //--A
		rmax.y = Max (rmax.y, r[n].y);
		rmax.z = Max (rmax.z, r[n].z);
		rmin.x = Min (rmin.x, r[n].x);
		rmin.y = Min (rmin.y, r[n].y);
	}
	rmin.z = 0.;
	VSet (threshold, 2., 2., 2.); //--A
	VVAdd (rmax, threshold); //--A
	VVSub (rmin, threshold); //--A
	VSub (box, rmax, rmin); //--A

	//output temperature.dat, calculate Tmax, heat_max, Tavg of 1/3 tip top.
	sprintf (filename, "out/temperature.dat");
	output = WriteFile (filename);
	if (T_flag == 0) {
		if (strcmp (Tmode, "single") == 0) {
			fprintf (output, "time(fs) Tmax(K) T1_3avg(K) heat_tot_max heat_nottingham_max heat_joule_max\n");
		} else if (strcmp (Tmode, "double") == 0) {
			fprintf (output, "time(fs) Tmax_phonon(K) T1_3avg_phonon(K) Tmax_electron(K) T1_3avg_electron(K) heat_tot_max heat_nottingham_max heat_joule_max\n");
		}
		T_flag = 1;
	}
	Tmax = Tele_max = heat_tot_max = heat_nottingham_max = heat_joule_max = -1.e100;
	T1_3avg = Tele1_3avg = 0.;
	i = 0;
	for (n = 0; n < N; n ++) {
		if (Tmax < T[n]) Tmax = T[n]; //--K
		if (Tele_max < T_electron[n]) Tele_max = T_electron[n]; //--K
		if (heat_tot_max < heat_tot[n]) heat_tot_max = heat_tot[n];
		if (heat_nottingham_max < heat_nottingham[n]) heat_nottingham_max = heat_nottingham[n];
		if (heat_joule_max < heat_joule[n]) heat_joule_max = heat_joule[n];	
		if (r[n].z > 2. / 3. * rmax.z) {
			T1_3avg += T[n]; //--K
			Tele1_3avg += T_electron[n]; //--K
			i ++;
		}
	}
	T1_3avg /= (double)i;
	Tele1_3avg /= (double)i;
	if (strcmp (Tmode, "single") == 0) {
		fprintf (output, "%f %f %f %e %e %e\n", (double)step * md_timestep, Tmax, T1_3avg, \
			 heat_tot_max, heat_nottingham_max, heat_joule_max);
	} else if (strcmp (Tmode, "double") == 0) {
		fprintf (output, "%f %f %f %f %f %e %e %e\n", (double)step * md_timestep, Tmax, T1_3avg, \
			 Tele_max, Tele1_3avg, heat_tot_max, heat_nottingham_max, heat_joule_max);
	}
	fclose (output);

	//compute the average temperature of every cell from femocs
	for (n = 0; n < N; n ++) {
		if (r[n].z > 0) {
			i = (int)((r[n].x - rmin.x) * Ncell.x / box.x);
			j = (int)((r[n].y - rmin.y) * Ncell.y / box.y);
			k = (int)((r[n].z - rmin.z) * Ncell.z / box.z);
			if (i < 0 || i >= Ncell.x || j < 0 || j >= Ncell.y || k < 0 || k >= Ncell.z) {
				printf ("error1(Main.cpp): temperature cell set [%d %d %d]  box: [%d %d %d]\n", \
					i, j, k, Vec(Ncell));
				exit (1);
			}
			Tcell[i][j][k].T += T[n]; //--K
			Tcell[i][j][k].N ++;
		}
	}
	DO_TCELL {
		cell = &Tcell[i][j][k];
		if (cell->N > 0) cell->T /= cell->N; //average
		else cell->T = 0.;
		cell->T /= TUnit; //dimensionless
	} END_TCELL

	//md kinetic energy
	VScale (rmax, 1. / (lUnit * 1.e10)); //dimensionless
	VScale (rmin, 1. / (lUnit * 1.e10)); //dimensionless
	VScale (box, 1. / (lUnit * 1.e10)); //dimensionless
	DO_MOL {
		if (mol[n].pedest != 1 && \
		    mol[n].r.z > minzInit && mol[n].r.z - minzInit < box.z-0.000001 && \
		    mol[n].r.x > rmin.x && mol[n].r.x < rmax.x-0.000001 && \
		    mol[n].r.y > rmin.y && mol[n].r.y < rmax.y-0.000001) {
			i = (int)((mol[n].r.x - rmin.x) * Ncell.x / box.x);
			j = (int)((mol[n].r.y - rmin.y) * Ncell.y / box.y);
			k = (int)((mol[n].r.z - minzInit) * Ncell.z / box.z);
			if (i < 0 || i >= Ncell.x || j < 0 || j >= Ncell.y || k < 0 || k >= Ncell.z) {
				printf ("error2(Main.cpp): temperature cell set [%d %d %d]  box: [%d %d %d]\n", \
					i, j, k, Vec(Ncell));
				exit (1);
			}
			Tcell[i][j][k].Eksum += 0.5 * mol[n].mass * VLenSq (mol[n].rv); //--dimensionless
			Tcell[i][j][k].Natoms ++;
		}
	}

	//temperature scaling
	DO_TCELL {
		Tcell[i][j][k].vFac = 1.;
	} END_TCELL
	DO_MOL {
		if (mol[n].pedest != 1 && \
		    mol[n].r.z > minzInit && mol[n].r.z - minzInit < box.z-0.000001 && \
		    mol[n].r.x > rmin.x && mol[n].r.x < rmax.x-0.000001 && \
		    mol[n].r.y > rmin.y && mol[n].r.y < rmax.y-0.000001) {
			i = (int)((mol[n].r.x - rmin.x) * Ncell.x / box.x);
			j = (int)((mol[n].r.y - rmin.y) * Ncell.y / box.y);
			k = (int)((mol[n].r.z - minzInit) * Ncell.z / box.z);
			cell = &Tcell[i][j][k];
			coff = sqrt (NDIM * (1. - 1. / cell->Natoms) * cell->T);
			cell->vFac = coff / sqrt (2. * cell->Eksum / cell->Natoms);
			if (isnormal(cell->vFac) == 0 || fabs(cell->vFac) < 1.e-100 || \
			    cell->Natoms == 0 || cell->N == 0) cell->vFac = 1.; 
			VScale (mol[n].rv, cell->vFac); 
		}
	}

	//calculate md kinetic energy after temperature scaling
	DO_TCELL {
		Tcell[i][j][k].Eksum = 0.;
	} END_TCELL
	DO_MOL {
		if (mol[n].pedest != 1 && \
		    mol[n].r.z > minzInit && mol[n].r.z - minzInit < box.z-0.000001 && \
		    mol[n].r.x > rmin.x && mol[n].r.x < rmax.x-0.000001 && \
		    mol[n].r.y > rmin.y && mol[n].r.y < rmax.y-0.000001) {
			i = (int)((mol[n].r.x - rmin.x) * Ncell.x / box.x);
			j = (int)((mol[n].r.y - rmin.y) * Ncell.y / box.y);
			k = (int)((mol[n].r.z - minzInit) * Ncell.z / box.z);
			Tcell[i][j][k].Eksum += 0.5 * mol[n].mass * VLenSq (mol[n].rv); //--dimensionless
		}
	}
	PrintTcellMovie ();

	free (Tcell);
	free (T);
	free (T_electron);
	free (heat_tot);
	free (heat_nottingham);
	free (heat_joule);
	free (r);
}

void PrintTcellMovie ()
{
	FILE *output;
	char filename[128], line[1024];
	int i, j, k, n;
	TCELL *cell;
	double T;

	sprintf (filename, "out/Tcell.movie");
	output = WriteFile (filename);
	fprintf (output, "%d\nproperties=pos:I:3:T(K):R:1:N:I:1:Ek(K):R:1:Natoms:I:1:scale:R:1\n", VProd (Ncell));
	DO_TCELL {
		cell = &Tcell[i][j][k];
		if (cell->Natoms > 0) T = 2. / NDIM * cell->Eksum / cell->Natoms * TUnit;
		else T = 0.;
		fprintf (output, "%d %d %d %f %d %f %d %f\n", i, j, k, cell->T * TUnit, cell->N, \
			 T, cell->Natoms, cell->vFac);
	} END_TCELL
	fclose (output);
}

void MaxwellStress (double theta) //--GPa
{
	int n;
	double zmax, zmin;

	theta /= (PUnit / 1.e9); //--dimensionless
	zmax = -1.e100;
	zmin = 1.e100;
	DO_MOL {
		VZero (mol[n].F);
		zmax = Max (zmax, mol[n].r.z);
		zmin = Min (zmin, mol[n].r.z);
	}
	DO_MOL {
		if (mol[n].flag == 2 && mol[n].r.z > maxzInit - (maxzInit - minzInit) / 5.) {
			VSet (mol[n].F, 0., 0., theta * M_PI * Sqr(rSurfCut) / mol[n].ligancy);
		} else VZero (mol[n].F);
	}
}
