#include "tipMaker.h"

Elem *unce, *suce;
VecI cellSize;
VecR unitRegion;
char CONTCAR_file[128], tipFilename[128], cellFilename[128], cellOrder[128], elemName[20][8];
double tip_r, tip_r1, tip_h, tip_theta;
int Nelems;

float Mould(double x, double y, double z)
{
	return sqrt(x*x+y*y+z*z);
}   

void CreatSuperCell ()
{
	FILE *input, *output;
	char line[1024], cellName[128], filename[128], command[128];
	int i, j, k, l, m, n, row, sum, N;
	int NatomsUnCe, NatomsSuCe;
	double scaling, xdata, ydata, zdata;
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, a, b, c; //lattice parameters
	float X_1, Y_1, Z_1; // Atomic positions in Cartesian coordinates
	float G_1x, G_1y, G_1z; // Atomic positions in supercell
	VecR region;

	//read input file
	input = ReadFile (CONTCAR_file);
	row = 1;
	while (1) {
		if (fgets (line, 1024, input) == NULL) break;
		switch (row) {
			case 1: 
				sscanf (line, "%s", cellName); 
				puts (cellName);
				break;
			case 2: 
				sscanf (line, "%lg", &scaling); 
				printf ("%f\n", scaling);
				break;
			case 3: 
				sscanf (line, "%lg %lg %lg", &a1, &a2, &a3);
				printf ("%f\t%f\t%f\n", a1, a2, a3);
				break;
			case 4: 
				sscanf (line, "%lg %lg %lg", &b1, &b2, &b3);
				printf ("%f\t%f\t%f\n", b1, b2, b3);
				break;
			case 5: 
				sscanf (line, "%lg %lg %lg", &c1, &c2, &c3);
				printf ("%f\t%f\t%f\n", c1, c2, c3);
				break;
			case 6:
				n = 0;
				while (line[0] == ' ') strcpy (line, line+1);
				while (1) {
					sscanf (line, "%s", elemName[n]);
					while (line[0] != ' ' && line[0] != '\0') strcpy (line, line+1);
					n ++;
					while (line[0] == ' ') strcpy (line, line+1);
					if (line[0] < 65 || line[0] > 122) break;
				}
				Nelems = n;
				for (n = 0; n < Nelems; n ++) printf ("%s ", elemName[n]);
				printf ("\n");
				AllocMem (unce, Nelems, Elem);
				AllocMem (suce, Nelems, Elem);
				for (n = 0; n < Nelems; n ++) sprintf (suce[n].name, "%s", elemName[n]);
				for (n = 0; n < Nelems; n ++) sprintf (unce[n].name, "%s", elemName[n]);
				break;
			case 7:
				for (n = 0; n < Nelems; n ++) {
					while (line[0] == ' ') strcpy (line, line+1);
					sscanf (line, "%d", &unce[n].Natoms);
					while (line[0] != ' ' && line[0] != '\0') strcpy (line, line+1);
					if (line[0] == '\0') break;
				}
				N = 0;
				for (n = 0; n < Nelems; n ++) {
					AllocMem (unce[n].r, unce[n].Natoms, VecR);
					N +=  unce[n].Natoms;
					printf("%d ", unce[n].Natoms);
				}
				printf ("\n");
				break;
			case 8: 
				puts (line);
				i = j = 0;
				break;
			default: 
				if (i == unce[j].Natoms) {
					j ++;
					i = 0;
				}
				sscanf (line, "%lg %lg %lg", &unce[j].r[i].x, &unce[j].r[i].y, &unce[j].r[i].z);
				printf ("%f %f %f\n", Vec(unce[j].r[i]));
				i ++;
				break;
		}
		row ++;
	}
	fclose (input);

	//output unit cell
	sprintf(filename,"out/Unitcell.xyz");
	output = WriteFile (filename);
	fprintf (output, "%d\nProperties=pos:R:3\n", N);
	VSet (unitRegion, a1, b2, c3);
	for (i = 0; i < Nelems; i ++) {
		for (j = 0; j < unce[i].Natoms; j ++) {
			xdata = unce[i].r[j].x;
			ydata = unce[i].r[j].y;
			zdata = unce[i].r[j].z;
			X_1 = xdata * a1 + ydata * b1 + zdata * c1;
			Y_1 = xdata * a2 + ydata * b2 + zdata * c2;
			Z_1 = xdata * a3 + ydata * b3 + zdata * c3;
//			printf ("%f\t%f\t%f\n", xdata, ydata, zdata);
			fprintf (output, "%f\t%f\t%f\n", X_1, Y_1, Z_1);
		}
	}
	fclose (output);

	//super cell initialize
	printf("The size of supercell are defined by three integers:\n");
	printf("Enter the size of supercell:%d %d %d\n", cellSize.x, cellSize.y, cellSize.z);
	sum = 0;
	for (i = 0; i < Nelems; i ++) sum += unce[i].Natoms;
	NatomsUnCe = sum; //the number of unit cell atoms
	NatomsSuCe = NatomsUnCe * VProd (cellSize);
	for (i = 0; i < Nelems; i ++) {
		suce[i].Natoms = unce[i].Natoms * VProd (cellSize);
		AllocMem (suce[i].r, suce[i].Natoms, VecR);
	}
	//creat super cell
	for (m = 0; m < Nelems; m ++) {
		i = 0;
		for (n = 0; n < unce[m].Natoms; n ++) {
			for(j=0;j<cellSize.x;j++){
				for(k=0;k<cellSize.y;k++){
					for(l=0;l<cellSize.z;l++){
						xdata = unce[m].r[n].x + j;
						ydata = unce[m].r[n].y + k;
						zdata = unce[m].r[n].z + l;
						X_1 = xdata * a1 + ydata * b1 + zdata * c1;
						Y_1 = xdata * a2 + ydata * b2 + zdata * c2;
						Z_1 = xdata * a3 + ydata * b3 + zdata * c3;
						VSet (suce[m].r[i], X_1, Y_1, Z_1); //--A
						VScale (suce[m].r[i], scaling);
						i ++;
					}
				}
			}
		}
	}

	//output super cell
	sprintf(filename,"out/Supercell.xyz");
	output = WriteFile (filename);
	fprintf(output,"%d\nSolutionReader properties=id:I:1:pos:R:3\n", NatomsSuCe);
	k = 0;
	for (i = 0; i < Nelems; i ++) {
		for (j = 0; j < suce[i].Natoms; j ++) {
			fprintf (output, "%d %f %f %f %s\n", k, suce[i].r[j].x, suce[i].r[j].y, suce[i].r[j].z, elemName[i]);
			k ++;
		}
	}
	fclose (output);

	//output POSCAR
	sprintf(filename,"out/POSCAR");
	output = WriteFile (filename);
	fprintf (output, "%s\n", cellName);
	fprintf (output, "\t1.0\n"); 
	fprintf (output, "\t%.10f\t%.10f\t%.10f\n", scaling*a1*cellSize.x, scaling*a2*cellSize.x, scaling*a3*cellSize.x);
	fprintf (output, "\t%.10f\t%.10f\t%.10f\n", scaling*b1*cellSize.y, scaling*b2*cellSize.y, scaling*b3*cellSize.y);
	fprintf (output, "\t%.10f\t%.10f\t%.10f\n", scaling*c1*cellSize.z, scaling*c2*cellSize.z, scaling*c3*cellSize.z);
	for (i = 0; i < Nelems; i ++) fprintf (output, "%s ", elemName[i]);
	fprintf (output, "\n");
	for (i = 0; i < Nelems; i ++) fprintf (output, "%d ", suce[i].Natoms);
	fprintf(output,"\nDirect\n");
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < unce[m].Natoms; n ++) {
			for(j=0;j<cellSize.x;j++){
				for(k=0;k<cellSize.y;k++){
					for(l=0;l<cellSize.z;l++){
						xdata = (unce[m].r[n].x + j) / cellSize.x;
						ydata = (unce[m].r[n].y + k) / cellSize.y;
						zdata = (unce[m].r[n].z + l) / cellSize.z;
						fprintf (output, "     %.10f     %.10f     %.10f\n", xdata, ydata, zdata);
					}
				}
			}
		}
	}
	fclose (output);
}

void TipMaker ()
{
	Atom *atom, atom_exchange;
	VecR vec;
	double r, h, theta, dz, x1, x2, h1, L, x, y, z, offset, temp, xmin, ymin, zmin;
	FILE *input, *output;
	char line[1024], filename[128];
	int m, n, i, j, k, range, N, *Natom, *count, Nold;

	printf ("\n=== tip maker...\n");

	r = tip_r; //--A
	h = tip_h; //--A
	theta = 90. - tip_theta; //degree
	theta *= M_PI / 180.; //rad

	//unit: A
	dz = r * cos(theta);
	x1 = r * sin(theta);
	h1 = h + dz - r;
	x2 = x1 + h1 * cot(theta);
	L = 2. * x2;
	VSet (vec, L, L, h);
	printf ("x2 = %f\n", x2);

	input = ReadFile (CONTCAR_file);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line,"%lg %lg %lg", &unitRegion.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &unitRegion.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &unitRegion.z);
	fclose (input);

	VDiv (cellSize, vec, unitRegion);
	VSet (vec, 1, 1, 1);
	VVAdd (cellSize, vec);
	VMul (region, unitRegion, cellSize);

	//creat rectangular super cell
	CreatSuperCell ();
	
	N = 0;
	for (m = 0; m < Nelems; m ++) N += suce[m].Natoms;
	AllocMem (atom, N, Atom);
	i = 0;
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < suce[m].Natoms; n ++) {
			VCopy (atom[i].r, suce[m].r[n]);
			sprintf (atom[i].element, "%s", suce[m].name);
			i ++;
		}
	}
	offset = 0.5; //--A
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= x2 + offset;
		atom[n].r.y -= x2 + offset;
		atom[n].r.z -= h1 - dz + offset;
	}

	//creat tip top
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (z > sqrt (r * r - x * x - y * y) && z > dz) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//creat tip body
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (z > (x2 - sqrt(x * x + y * y)) * tan(theta) - (h1 - dz)) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//move
	xmin = ymin = zmin = 1.e100;
	for (n = 0; n < N; n ++) {
		xmin = Min(xmin, atom[n].r.x);
		ymin = Min(ymin, atom[n].r.y);
		zmin = Min(zmin, atom[n].r.z);
	}
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= xmin;
		atom[n].r.y -= ymin;
		atom[n].r.z -= zmin;
	}

	//the initial elements are randomly distrubuted
	for (n = 0; n < N; n ++) {
		i = rand () % N;
		atom_exchange = atom[n];
		atom[n] = atom[i];
		atom[i] = atom_exchange;
	}
	if (Nelems > 1 && strcmp (cellOrder, "disorder") == 0) {
		Nold = 0;
		for (i = 0; i < Nelems; i ++) Nold += suce[i].Natoms;
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		srand (randSeed);
		for (n = 0; n < N; n ++) {
			j = rand () % Nold;
			range = 0;
			for (k = 0; k < Nelems; k ++) {
				if (j >= range && j < range + suce[k].Natoms) break;
				range += suce[k].Natoms;
			}
			if (count[k+1] >= (double)N*(double)suce[k].Natoms/(double)Nold) {
				n --;
				continue;
			}
			strcpy (atom[n].element, suce[k].name);
			count[k+1] ++;
			count[0] ++;
		}
		free (count);
	}

	//atom number
	AllocMem (Natom, Nelems, int);
	for (n = 0; n < Nelems; n ++) Natom[n] = 0;
	for (n = 0; n < N; n ++) {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp (atom[n].element, elemName[i]) == 0) Natom[i] ++;
		}
	}
	printf ("total atoms is %d ", N);
	for (i = 0; i < Nelems; i ++) printf ("%s %d %.2f%% ", elemName[i], Natom[i], (double)Natom[i]/(double)N*100.);
	printf ("\n");

	if ((output = fopen (tipFilename, "w")) == NULL){
		printf ("\ncreat %s error\n", tipFilename);
		getchar ();
		exit (1);
	}
	fprintf (output, "%d\n", N);
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");
	for (n = 0; n < N; n ++) {
		fprintf (output, "%s %.10f %.10f %.10f\n", atom[n].element, atom[n].r.x, atom[n].r.y, atom[n].r.z);
	}
	fclose (output);

	free (unce);
	free (suce);
	free (atom);
	free (Natom);
}

void MushroomTipMaker ()
{
	Atom *atom, atom_exchange;
	VecR vec;
	double r, r1, h, theta, dz, x1, x2, h1, L, x, y, z, offset, temp, xmin, ymin, zmin;
	FILE *input, *output;
	char line[1024], filename[128];
	int m, n, i, j, k, range, N, *Natom, *count, Nold;

	printf ("\n=== mushroom tip maker...\n");

	r = tip_r; //--A
	r1 = tip_r1; //--A
	if (r <= r1) {
		printf ("r should > r1\n");
		exit (1);
	}
	h = tip_h; //--A
	theta = 90. - tip_theta; //degree
	theta *= M_PI / 180.; //rad

	//unit: A
	dz = sqrt (r * r - r1 * r1);
	h1 = h - dz - r;
	x2 = r1 + h1 * cot(theta);
	L = 2. * Max (r, x2);
	VSet (vec, L, L, h);
	printf ("x2 = %f\n", x2);

	input = ReadFile (CONTCAR_file);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line,"%lg %lg %lg", &unitRegion.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &unitRegion.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &unitRegion.z);
	fclose (input);

	VDiv (cellSize, vec, unitRegion);
	VSet (vec, 1, 1, 1);
	VVAdd (cellSize, vec);
	VMul (region, unitRegion, cellSize);

	//creat rectangular super cell
	CreatSuperCell ();
	
	N = 0;
	for (m = 0; m < Nelems; m ++) N += suce[m].Natoms;
	AllocMem (atom, N, Atom);
	i = 0;
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < suce[m].Natoms; n ++) {
			VCopy (atom[i].r, suce[m].r[n]);
			sprintf (atom[i].element, "%s", suce[m].name);
			i ++;
		}
	}
	offset = 0.5; //--A
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= x2 + offset;
		atom[n].r.y -= x2 + offset;
		atom[n].r.z -= h1 - dz + offset;
	}

	//creat tip top
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (sqrt (x * x + y * y + Sqr (z - dz)) > r && z > 0) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//creat tip body
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (z > (x2 - sqrt(x * x + y * y)) * tan(theta) - h1 && z <= 0) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//move
	xmin = ymin = zmin = 1.e100;
	for (n = 0; n < N; n ++) {
		xmin = Min(xmin, atom[n].r.x);
		ymin = Min(ymin, atom[n].r.y);
		zmin = Min(zmin, atom[n].r.z);
	}
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= xmin;
		atom[n].r.y -= ymin;
		atom[n].r.z -= zmin;
	}

	//the initial elements are randomly distrubuted
	for (n = 0; n < N; n ++) {
		i = rand () % N;
		atom_exchange = atom[n];
		atom[n] = atom[i];
		atom[i] = atom_exchange;
	}
	if (Nelems > 1 && strcmp (cellOrder, "disorder") == 0) {
		Nold = 0;
		for (i = 0; i < Nelems; i ++) Nold += suce[i].Natoms;
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		srand (randSeed);
		for (n = 0; n < N; n ++) {
			j = rand () % Nold;
			range = 0;
			for (k = 0; k < Nelems; k ++) {
				if (j >= range && j < range + suce[k].Natoms) break;
				range += suce[k].Natoms;
			}
			if (count[k+1] >= (double)N*(double)suce[k].Natoms/(double)Nold) {
				n --;
				continue;
			}
			strcpy (atom[n].element, suce[k].name);
			count[k+1] ++;
			count[0] ++;
		}
		free (count);
	}

	//atom number
	AllocMem (Natom, Nelems, int);
	for (n = 0; n < Nelems; n ++) Natom[n] = 0;
	for (n = 0; n < N; n ++) {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp (atom[n].element, elemName[i]) == 0) Natom[i] ++;
		}
	}
	printf ("total atoms is %d ", N);
	for (i = 0; i < Nelems; i ++) printf ("%s %d %.2f%% ", elemName[i], Natom[i], (double)Natom[i]/(double)N*100.);
	printf ("\n");

	if ((output = fopen (tipFilename, "w")) == NULL){
		printf ("\ncreat %s error\n", tipFilename);
		getchar ();
		exit (1);
	}
	fprintf (output, "%d\n", N);
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");
	for (n = 0; n < N; n ++) {
		fprintf (output, "%s %.10f %.10f %.10f\n", atom[n].element, atom[n].r.x, atom[n].r.y, atom[n].r.z);
	}
	fclose (output);

	free (unce);
	free (suce);
	free (atom);
	free (Natom);
}

void SuperCellMaker ()
{
	Atom *atom, atom_exchange;
	FILE *input, *output;
	int m, n, i, j, k, range, N, *Natom, *count;
	char line[1024];
	double temp, xmin, ymin, zmin;

	printf ("\n=== super cell maker...\n");

	input = ReadFile (CONTCAR_file);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line,"%lg %lg %lg", &unitRegion.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &unitRegion.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &unitRegion.z);
	fclose (input);

	VCopy (cellSize, initUcell);
	VMul (region, unitRegion, cellSize);
	//creat rectangular super cell
	CreatSuperCell ();

	N = 0;
	for (m = 0; m < Nelems; m ++) N += suce[m].Natoms;
	AllocMem (atom, N, Atom);
	i = 0;
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < suce[m].Natoms; n ++) {
			VCopy (atom[i].r, suce[m].r[n]);
			sprintf (atom[i].element, "%s", suce[m].name);
			i ++;
		}
	}

	//move
	xmin = ymin = zmin = 1.e100;
	for (n = 0; n < N; n ++) {
		xmin = Min(xmin, atom[n].r.x);
		ymin = Min(ymin, atom[n].r.y);
		zmin = Min(zmin, atom[n].r.z);
	}
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= xmin;
		atom[n].r.y -= ymin;
		atom[n].r.z -= zmin;
	}

	//the initial elements are randomly distrubuted
	for (n = 0; n < N; n ++) {
		i = rand () % N;
		atom_exchange = atom[n];
		atom[n] = atom[i];
		atom[i] = atom_exchange;
	}
	if (Nelems > 1 && strcmp (cellOrder, "disorder") == 0) {
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		srand (randSeed);
		for (n = 0; n < N; n ++) {
			j = rand () % N;
			range = 0;
			for (k = 0; k < Nelems; k ++) {
				if (j >= range && j < range + suce[k].Natoms) break;
				range += suce[k].Natoms;
			}
			if (count[k+1] >= (double)suce[k].Natoms) {
				n --;
				continue;
			}
			strcpy (atom[n].element, suce[k].name);
			count[k+1] ++;
			count[0] ++;
		}
		free (count);
	}

	//atom number
	AllocMem (Natom, Nelems, int);
	for (n = 0; n < Nelems; n ++) Natom[n] = 0;
	for (n = 0; n < N; n ++) {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp (atom[n].element, elemName[i]) == 0) Natom[i] ++;
		}
	}
	printf ("total atoms is %d ", N);
	for (i = 0; i < Nelems; i ++) printf ("%s %d %.2f%% ", elemName[i], Natom[i], (double)Natom[i]/(double)N*100.);
	printf ("\n");

	if ((output = fopen (cellFilename, "w")) == NULL){
		printf ("\ncreat %s error\n", cellFilename);
		getchar ();
		exit (1);
	}
	fprintf (output, "%d\n", N);
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");
	for (n = 0; n < N; n ++) {
		fprintf (output, "%s %.10f %.10f %.10f\n", atom[n].element, atom[n].r.x, atom[n].r.y, atom[n].r.z);
	}
	fclose (output);

	free (Natom);
	free (atom);
}

void ExtensionMaker ()
{
	double angle;
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	double rad,x,y,z;
	double zmax=0.00000000;
	double zmin;
	printf("Enter the zmin:\n");
	scanf("%lf",&zmin);
//	zmin = -20000;
	//the radius of the top
	double R0= 0.;  
	double lat=0.; 
	printf ("Enter the radius of extension top:\n");
	scanf ("%lg", &R0);
//	R0 = 10537.801253;
	printf ("Enter the lattice parameter:\n");
	scanf ("%lg", &lat);
//	lat = 300;
//	double theta = 3.;
	double theta;
	printf ("Enter the angle of extension top [deg]:\n");
	scanf ("%lg", &theta);
	double width = lat * 100.;
// 54.52762682;
	double zd2=zmin+50.00000000*lat/3.615;
	double zd1=zmin+120.00000000*lat/3.615;
	double zd0=zmin+150.00000000*lat/3.615;	
	printf ("%.3f %.3f %.3f\n", zd2, zd1, zd0);
	double R,Rd,a,r,d;
	double a1=3.41116400*lat/3.615;
	double a2=20.80240900*lat/3.615;
	double a3=91.73219800*lat/3.615;	
	char file[128];
	FILE *output;
	sprintf(file,"in/extension.xyz");
	output=fopen(file,"w");
//print natoms	
	z=zmax;
	R=R0;
	for (k=0;k<4;k++){
		if (k>0){
			R=R+tan(theta*M_PI/180.)*lat*k/4.00000000;
			z=z-lat*k/4.00000000;
			}
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			j=j+1;
                               		    	    	 }
			}
//cone part
	z=z-lat;
	for (z;z>=zd0;z=z-lat){
		R=R+tan(theta*M_PI/180.)*lat;
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			j=j+1;
                               		    	     }
			 	 }
	Rd=R;
	a=(zd0-zd1)/pow(a1,2);
	for (z;z>=zd1;z=z-lat){
		R=Rd+a1-sqrt((z-zd1)/a);
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			j=j+1;
                               		    	     }
			 	 	}
	Rd=R;
	a=(zd1-zd2)/pow(a2,2);
	for (z;z>=zd2;z=z-lat){
		R=Rd+a2-sqrt((z-zd2)/a);
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			j=j+1;
                               		    	     }
	
			 	 }
	Rd=R;
	a=(zd2-zmin)/pow(a3,2);
	for (z;z>=zmin;z=z-lat){
		r=R;
		R=Rd+a3-sqrt((z-zmin)/a);
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			j=j+1;
                               		    	     }
			 	 }
	d=R-r;
	for(R;R<=width*sqrt(2.00000000);){
		R=R+d;
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			if (x<=width&&x>=-width&&y<=width&&y>=-width){
				j=j+1;
				   }
                               			    	     }
			 	 		}
	printf("%d\n",j);

	/*-------------------------------------------------------------------*/
	/*-------------------------production--------------------------------*/
	/*-------------------------------------------------------------------*/

	fprintf(output,"%d\n",j);
	fprintf(output,"Properties=id:I:1:pos:R:3:Type:R:1\n");
	
//print atoms
	z=zmax;
	R=R0;

//Connect part
	for (k=0;k<4;k++){
		if (k>0){
			R=R+tan(theta*M_PI/180.)*lat*k/4.00000000;
			z=z-lat*k/4.00000000;
			}
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			i=i+1;
			fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
                               		    	    	 }
			}

//cone part
	z=z-lat;
	for (z;z>=zd0;z=z-lat){
		R=R+tan(theta*M_PI/180.)*lat;
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			i=i+1;
			fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
                               		    	     }
			 	 }
//base connect part1
	Rd=R;
	a=(zd0-zd1)/pow(a1,2);
	for (z;z>=zd1;z=z-lat){
		R=Rd+a1-sqrt((z-zd1)/a);
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			i=i+1;
			fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
                               		    	     }
			 	 	}
//base connect part2
	Rd=R;
	a=(zd1-zd2)/pow(a2,2);
	for (z;z>=zd2;z=z-lat){
		R=Rd+a2-sqrt((z-zd2)/a);
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			i=i+1;
			fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
                               		    	     }
	
			 	 }

//base connect part3
	Rd=R;
	a=(zd2-zmin)/pow(a3,2);
	for (z;z>=zmin;z=z-lat) {
		r=R;
		R=Rd+a3-sqrt((z-zmin)/a);
		for(angle=0;angle<360;angle=angle+4.80000000) {
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			i=i+1;
			fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
        	}
	}

//base part1
	d=R-r;
	for(R;R<=width*sqrt(2.00000000);){
		R=R+d;
		for(angle=0;angle<360;angle=angle+4.80000000){
			rad = angle*M_PI/180;
			x=R*cos(rad);
			y=R*sin(rad);
			z=z;
			if (x<=width&&x>=-width&&y<=width&&y>=-width){
				i=i+1;
				fprintf(output,"%d\t%lf\t%lf\t%lf\t%d\n",i,x,y,z,1);
				   }
                               			    	     }
			 	 		}
	fclose(output);
}

void ProlateSpheroidalMaker ()
{
	FILE *output;
	char filename[128];
	double a, u, v, phi, x, y, z, theta, H, r;
	int i, j, k, N, N1, N2;

	sprintf (filename, "%s", tipFilename);
	output = fopen (filename, "w");

//	printf ("input the radius of tip r [A]:\n");
//	scanf ("%lg", &r);
	if (fabs(tip_r1) > 1.e-60) r = tip_r1;
	if (fabs(tip_r) > 1.e-60) r = tip_r;
//	printf ("input tip half-angle theta [deg]:\n");
//	scanf ("%lg", &theta);
	theta = tip_theta;
	theta *= M_PI / 180.; //--rads
	a = r / (sin(theta) * tan(theta));
	printf ("the half of the foci distance a [A]: %.3f\n", a);
	printf ("the tip to anode distance d[A]: %.3f\n", a * cos(theta));
//	printf ("input tip hight H [A]:\n");
//	scanf ("%lg", &H);
	H = tip_h;
//	printf ("the number of points [Nz Nxy]:(eg. 1000 100)\n");
//	scanf ("%d %d", &N1, &N2);
	N1 = 200;
	N2 = 100;

	v = M_PI - theta;
	N = N1 * N2;
	fprintf (output, "%d\n", N);
	u = (a * cos(v) - H) / (a * cos(v));
	u = log (u + sqrt(u * u - 1.));
	region.x = region.y = a * sinh(u) * sin(v) * 2. + 10.;
	region.z = H + 10.;
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");

	for (i = 1; i <= N1; i ++) {
		for (k = 1; k <= N2; k ++) {
			u = (a * cos(v) - H) / (a * cos(v));
			u = log (u + sqrt(u * u - 1.)) * (double)i / (double)N1;
			phi = k * 2. * M_PI / (double)N2;
			x = a * sinh(u) * sin(v) * cos(phi) + region.x / 2.;
			y = a * sinh(u) * sin(v) * sin(phi) + region.y / 2.;
			z = a * cosh(u) * cos(v) - (a * cos(v) - H);
			fprintf (output, "point %e %e %e\n", x, y, z);
		}
	}
	printf ("x2 = %f A\n", a * sinh(u) * sin(v));

	fclose (output);
}

void EllipticMaker ()
{
	FILE *output;
	char filename[128];
	double a, b, c, x, y, z, H, Rb, theta, phi;
	int i, j, k, N, N1, N2;

	sprintf (filename, "%s", tipFilename);
	output = fopen (filename, "w");

//	printf ("input tip hight H [A]:\n");
//	scanf ("%lg", &H);
	H = tip_h;
//	printf ("the radius of bottom Rb [A]:\n");
//	scanf ("%lg", &Rb);
	if (fabs(tip_r1) > 1.e-60) Rb = tip_r1;
	if (fabs(tip_r) > 1.e-60) Rb = tip_r;
//	printf ("the number of points [Nz Nxy]:(eg. 200 200)\n");
//	scanf ("%d %d", &N1, &N2);
	N1 = 200;
	N2 = 100;
	a = b = Rb; //--A
	c = H; //--A

	N = N1 * N2;
	fprintf (output, "%d\n", N);
	region.x = region.y = Rb * 2. + 10.;
	region.z = H + 10.;
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");

	for (i = 1; i <= N1; i ++) {
		for (k = 1; k <= N2; k ++) {
			theta = 2. * M_PI * i / N1;
			phi = M_PI / 2. * k / N2;
			x = a * sin(theta) * sin(phi) + a;
			y = b * cos(theta) * sin(phi) + b;
			z = c * cos(phi);
			fprintf (output, "point %e %e %e\n", x, y, z);
		}
	}
	printf ("x2 = %f A\n", Rb);

	fclose (output);
}

void ConeTipMaker ()
{
	Atom *atom, atom_exchange;
	VecR vec;
	double r, r1, h, theta, dz, x1, x2, h1, L, x, y, z, offset, temp, xmin, ymin, zmin;
	FILE *input, *output;
	char line[1024], filename[128];
	int m, n, i, j, k, range, N, *Natom, *count, Nold;

	printf ("\n=== Cone tip maker...\n");

	h = tip_h; //--A
	theta = 90. - tip_theta; //degree
	theta *= M_PI / 180.; //rad

	//unit: A
	h1 = tip_h;
	x2 = r = h * cot (theta);
	L = 2. * r;
	VSet (vec, L, L, h);
	printf ("x2 = %f\n", x2);

	input = ReadFile (CONTCAR_file);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line,"%lg %lg %lg", &unitRegion.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &unitRegion.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &unitRegion.z);
	fclose (input);

	VDiv (cellSize, vec, unitRegion);
	VSet (vec, 1, 1, 1);
	VVAdd (cellSize, vec);
	VMul (region, unitRegion, cellSize);

	//creat rectangular super cell
	CreatSuperCell ();
	
	N = 0;
	for (m = 0; m < Nelems; m ++) N += suce[m].Natoms;
	AllocMem (atom, N, Atom);
	i = 0;
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < suce[m].Natoms; n ++) {
			VCopy (atom[i].r, suce[m].r[n]);
			sprintf (atom[i].element, "%s", suce[m].name);
			i ++;
		}
	}

	offset = - 0.5; //--A
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= x2 + offset;
		atom[n].r.y -= x2 + offset;
		atom[n].r.z -= h1 + offset;
	}

	//creat tip body
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (z > (x2 - sqrt(x * x + y * y)) * tan(theta) - h1) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//move
	xmin = ymin = zmin = 1.e100;
	for (n = 0; n < N; n ++) {
		xmin = Min(xmin, atom[n].r.x);
		ymin = Min(ymin, atom[n].r.y);
		zmin = Min(zmin, atom[n].r.z);
	}
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= xmin;
		atom[n].r.y -= ymin;
		atom[n].r.z -= zmin;
	}

	//the initial elements are randomly distrubuted
	for (n = 0; n < N; n ++) {
		i = rand () % N;
		atom_exchange = atom[n];
		atom[n] = atom[i];
		atom[i] = atom_exchange;
	}
	if (Nelems > 1 && strcmp (cellOrder, "disorder") == 0) {
		Nold = 0;
		for (i = 0; i < Nelems; i ++) Nold += suce[i].Natoms;
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		srand (randSeed);
		for (n = 0; n < N; n ++) {
			j = rand () % Nold;
			range = 0;
			for (k = 0; k < Nelems; k ++) {
				if (j >= range && j < range + suce[k].Natoms) break;
				range += suce[k].Natoms;
			}
			if (count[k+1] >= (double)N*(double)suce[k].Natoms/(double)Nold) {
				n --;
				continue;
			}
			strcpy (atom[n].element, suce[k].name);
			count[k+1] ++;
			count[0] ++;
		}
		free (count);
	}

	//atom number
	AllocMem (Natom, Nelems, int);
	for (n = 0; n < Nelems; n ++) Natom[n] = 0;
	for (n = 0; n < N; n ++) {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp (atom[n].element, elemName[i]) == 0) Natom[i] ++;
		}
	}
	printf ("total atoms is %d ", N);
	for (i = 0; i < Nelems; i ++) printf ("%s %d %.2f%% ", elemName[i], Natom[i], (double)Natom[i]/(double)N*100.);
	printf ("\n");

	if ((output = fopen (tipFilename, "w")) == NULL){
		printf ("\ncreat %s error\n", tipFilename);
		getchar ();
		exit (1);
	}
	fprintf (output, "%d\n", N);
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");
	for (n = 0; n < N; n ++) {
		fprintf (output, "%s %.10f %.10f %.10f\n", atom[n].element, atom[n].r.x, atom[n].r.y, atom[n].r.z);
	}
	fclose (output);

	free (unce);
	free (suce);
	free (atom);
	free (Natom);
}

void CylinderTipMaker ()
{
	Atom *atom, atom_exchange;
	VecR vec;
	double r, r1, h, theta, dz, x1, x2, h1, L, x, y, z, offset, temp, xmin, ymin, zmin;
	FILE *input, *output;
	char line[1024], filename[128];
	int m, n, i, j, k, range, N, *Natom, *count, Nold;

	printf ("\n=== Cylinder tip maker...\n");

	h = tip_h; //--A
	if (fabs(tip_r1) > 1.e-60) r = tip_r1;
	if (fabs(tip_r) > 1.e-60) r = tip_r;

	//unit: A
	h1 = tip_h;
	x2 = r;
	L = 2. * r;
	VSet (vec, L, L, h);
	printf ("x2 = %f\n", x2);

	input = ReadFile (CONTCAR_file);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	fgets (line, 1024, input);
	sscanf (line,"%lg %lg %lg", &unitRegion.x, &temp, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &unitRegion.y, &temp);
	fgets (line, 1024, input);
	sscanf (line, "%lg %lg %lg", &temp, &temp, &unitRegion.z);
	fclose (input);

	VDiv (cellSize, vec, unitRegion);
	VSet (vec, 1, 1, 1);
	VVAdd (cellSize, vec);
	VMul (region, unitRegion, cellSize);

	//creat rectangular super cell
	CreatSuperCell ();
	
	N = 0;
	for (m = 0; m < Nelems; m ++) N += suce[m].Natoms;
	AllocMem (atom, N, Atom);
	i = 0;
	for (m = 0; m < Nelems; m ++) {
		for (n = 0; n < suce[m].Natoms; n ++) {
			VCopy (atom[i].r, suce[m].r[n]);
			sprintf (atom[i].element, "%s", suce[m].name);
			i ++;
		}
	}

	offset = 0.5; //--A
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= x2 + offset;
		atom[n].r.y -= x2 + offset;
		atom[n].r.z -= h1 + offset;
	}

	//creat tip body
	for (n = 0; n < N; n ++) {
		x = atom[n].r.x;
		y = atom[n].r.y;
		z = atom[n].r.z;
		if (x * x + y * y > r * r) atom[n].delete_flag = 1;
		else atom[n].delete_flag = 0;
	}
	i = 0;
	for (n = 0; n < N; n ++) {
		if (atom[n].delete_flag == 0) {
			atom[i] = atom[n];
			i ++;
		}
	}
	N = i;

	//move
	xmin = ymin = zmin = 1.e100;
	for (n = 0; n < N; n ++) {
		xmin = Min(xmin, atom[n].r.x);
		ymin = Min(ymin, atom[n].r.y);
		zmin = Min(zmin, atom[n].r.z);
	}
	for (n = 0; n < N; n ++) {
		atom[n].r.x -= xmin;
		atom[n].r.y -= ymin;
		atom[n].r.z -= zmin;
	}

	//the initial elements are randomly distrubuted
	for (n = 0; n < N; n ++) {
		i = rand () % N;
		atom_exchange = atom[n];
		atom[n] = atom[i];
		atom[i] = atom_exchange;
	}
	if (Nelems > 1 && strcmp (cellOrder, "disorder") == 0) {
		Nold = 0;
		for (i = 0; i < Nelems; i ++) Nold += suce[i].Natoms;
		AllocMem (count, Nelems+1, int);
		for (i = 0; i <= Nelems; i ++) count[i] = 0;
		srand (randSeed);
		for (n = 0; n < N; n ++) {
			j = rand () % Nold;
			range = 0;
			for (k = 0; k < Nelems; k ++) {
				if (j >= range && j < range + suce[k].Natoms) break;
				range += suce[k].Natoms;
			}
			if (count[k+1] >= (double)N*(double)suce[k].Natoms/(double)Nold) {
				n --;
				continue;
			}
			strcpy (atom[n].element, suce[k].name);
			count[k+1] ++;
			count[0] ++;
		}
		free (count);
	}

	//atom number
	AllocMem (Natom, Nelems, int);
	for (n = 0; n < Nelems; n ++) Natom[n] = 0;
	for (n = 0; n < N; n ++) {
		for (i = 0; i < Nelems; i ++) {
			if (strcmp (atom[n].element, elemName[i]) == 0) Natom[i] ++;
		}
	}
	printf ("total atoms is %d ", N);
	for (i = 0; i < Nelems; i ++) printf ("%s %d %.2f%% ", elemName[i], Natom[i], (double)Natom[i]/(double)N*100.);
	printf ("\n");

	if ((output = fopen (tipFilename, "w")) == NULL){
		printf ("\ncreat %s error\n", tipFilename);
		getchar ();
		exit (1);
	}
	fprintf (output, "%d\n", N);
	fprintf (output, "Lattice=\"%.10f 0. 0. 0. %.10f 0. 0. 0. %.10f\" ", region.x, region.y, region.z);
	fprintf (output, "SolutionReader properties=species:S:1:pos:R:3\n");
	for (n = 0; n < N; n ++) {
		fprintf (output, "%s %.10f %.10f %.10f\n", atom[n].element, atom[n].r.x, atom[n].r.y, atom[n].r.z);
	}
	fclose (output);

	free (unce);
	free (suce);
	free (atom);
	free (Natom);
}
