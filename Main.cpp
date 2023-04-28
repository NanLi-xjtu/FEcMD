/*
 *  Main function for FEcMD
 *
 *  Field emission coupled with molecular dynamics simulation (FEcMD) software package is a computational tool 
 *  for studying the electron emission characteristics and the atomic structure evolution of micro- and nano-protrusions made of 
 *  pure metals or multi-component alloys by means of multi-physics and multi-scale methodology. 
 *
 *  Created on: 2023/4/1
 *  Author: NanLi
 */
#include "main.h"

using namespace std;
using namespace femocs;

int time_flag = 0;

int main (int argc, char **argv) 
{
	double start, start1, end, movie_timestep, md_timestep, time_limit;
	FILE *femocs_in, *output;
	char line[1024], name[128], Tmode[128];

	start = omp_get_wtime();

      printf("\n");
      printf("                                \\\\--//          \n");
      printf("                              (|-@ @--|)\n");
      printf("                               \\ (_)  /      \n");
      printf("                                  -         \n");
      printf("+------------o00o-------------------------------+\n");
      printf("|                                               |\n");
      printf("|                    FEcMD                      |\n");
      printf("|        Written by Nan Li && Bing Xiao         |\n");
      printf("|                  2023-04-1                    |\n");
      printf("| ln906061119@stu.xjtu.edu.cn, Nan Li from XJTU |\n");
      printf("|  bingxiao84@xjtu.edu.cn, Bing Xiao from XJTU  |\n");
      printf("|                                               |\n");
      printf("+----------------------------------oOOo---------+\n");
      printf("\n");

/*************************************************************************************/
/******************************************MD*****************************************/
/***************************************configure*************************************/
	int time, coreNum = omp_get_num_procs ();
	int n, frame;
	double startMD, endMD;

	system ("rm -rf out");
	system ("mkdir out");
	sprintf (mdXYZFilename, "in/mdlat.in.xyz");
	if (strcmp (model, "ED") != 0) {
		sprintf (CONTCAR_file, "in/CONTCAR");
		sprintf (mdPamsFilename, "in/md.in");
		PrintOpen ();
		GetParaValue ();
		switch (argc)
		{
			case 1:
				sprintf (model, "ED-MD");
				break;
			case 2:	
				if (strcmp (argv[1], "PH") == 0) {
					TipMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "MH") == 0) {
					MushroomTipMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "extension") == 0) { 
					ExtensionMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "PS") == 0) {
					ProlateSpheroidalMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "HE") == 0) {
					EllipticMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "Cone") == 0) {
					ConeTipMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "Cylinder") == 0) {
					CylinderTipMaker ();
					system ("rm -rf out");
				} else if(strcmp (argv[1], "supercell") == 0) { 
					SuperCellMaker ();
					system ("rm -rf out");
				} else if (strcmp (argv[1], "ED") == 0) {
					sprintf (model, "ED");
					break;
				} else if (strcmp (argv[1], "MD") == 0) {
					sprintf (model, "MD");
					break;
				} else if (strcmp (argv[1], "Maxwell") == 0) {
					sprintf (model, "MD");
					sprintf (force_sort, "Maxwell");
					break;
				} else if (strcmp (argv[1], "heattransport") == 0) {
					sprintf (model, "MD");
					sprintf (nonq_mode, "heattransport");
					break;
				} else {
					printf ("error: input format\n");
					printf ("./run\n      PH\n      MH\n      extension\n      PS\n      HE\n      Cone\n      Cylinder\n");
					printf ("      ED\n      MD\n      Maxwell\n      heattransport\n");
					exit (1);
				}
				return 0;
			default:
				printf ("error: input format\n");
				exit (1);
		}
		moreCycles = 1;
	}
	SetParams (argc, argv);
	SetupJob (argc, argv);
	if (strcmp (model, "MD") == 0) {
		printf ("=== MD simulation...\n");
		printf ("\ntotal steps: %d\n", stepLimit);
		PRINTSUMMARYTITLE ();
		for (int i = 1; i <= stepLimit; i ++) SingleStep ();
		return 0;
	}

/*************************************************************************************/
/*****************************************FEMOCS**************************************/
/***************************************configure*************************************/
	string mode = "default";
	char arg[128];
	int success = 0;
	bool add_rnd_noise = false;
	int n_iterations = 1;// determine number of iterations
	int n_atoms = 0;

	PrintFemocsin_xyz ();
	system ("cp out/md/femocs.in.xyz in/");
	string filename = "in/femocs.in";
	femocs::Femocs femocs(filename);
	success = system("rm -rf md.in.tmp");
	string cmd1 = "infile"; string infile = "";
	
	success = femocs.parse_command(cmd1, infile);
	print_progress("\n> reading " + cmd1, infile != "");

	sprintf (name, "in/femocs.in");
	femocs_in = ReadFile (name);
	while (1) {
		if (fgets (line, 1024, femocs_in) == NULL) break;
		GetDoubleVariable (line, "movie_timestep", &movie_timestep);
		GetDoubleVariable (line, "timestep", &md_timestep);
		GetDoubleVariable (line, "timelimit", &time_limit);
		GetCharVariable (line, "temperature_mode", Tmode);
	}
	fclose (femocs_in);
	system ("mkdir out/mesh && mkdir out/md");
	system ("cp in/femocs.in.xyz out/md/");
		success = 0;
	if (infile != "") read_n_atoms(infile, n_atoms);

/**************************************************************************************/
/**************************************Main Loop***************************************/
/**************************************************************************************/
	n_iterations = (int)(time_limit * 1000. / md_timestep);
	printf ("\ntotal steps: %d\n", n_iterations);	
	if (strcmp (model, "ED") != 0) {
		printf ("=== MD simulation...\n");
		PRINTSUMMARYTITLE ();
		for (n = 0; n < stepEquil; n ++) SingleStep ();
	}
	frame = 0;
	for (int i = 1; i <= n_iterations; ++i) {
		start1 = omp_get_wtime ();
		if (n_iterations > 1) cout << "\n> iteration " << i << endl;
/*------------------------------------------------------------------------------------*/
/*-----------------------------------------MD-----------------------------------------*/
		if (strcmp (model, "ED-MD") == 0) {
			printf ("=== MD simulation...\n");
			startMD = omp_get_wtime ();
			PRINTSUMMARYTITLE ();
			SingleStep ();
			PrintFemocsin_xyz ();
			endMD = omp_get_wtime ();
			printf ("time: %.3f\n", endMD - startMD);
		}
/*------------------------------------------------------------------------------------*/
/*---------------------------------------FEMOCS---------------------------------------*/
		success = femocs.import_atoms(infile, add_rnd_noise);
		success += femocs.run();
/*------------------------------------------------------------------------------------*/
/*-------------------------------------INTERFACE--------------------------------------*/
		if (strcmp (model, "ED-MD") == 0) {
			ForceInterface ();
			TemperatureInterface (i);
		}
		if ((double)i * md_timestep > (double)frame * movie_timestep) {
			system ("cat out/emission.xyz >> out/emission.movie");
			system ("cat out/electrons.xyz >> out/electrons.movie");
			system ("cat out/surface_fields.xyz >> out/surface_fields.movie");
			system ("cat out/surface_temperatures.xyz >> out/surface_temperatures.movie");
			system ("cat out/forces.xyz >> out/forces.movie");
			system ("cat out/fields.xyz >> out/fields.movie");
			system ("cat out/ch_solver.xyz >> out/ch_solver.movie");
			system ("cat out/temperature_phonon.xyz >> out/temperature_phonon.movie");
			system ("cat out/exchange-correlation.xyz >> out/exchange-correlation.movie");
			system ("cat out/getelec.dat >> out/getelec.movie.dat");
			system ("mv out/*.vtk out/mesh/");
			frame ++;
		}
		system ("rm out/*.txt");
		

		end = omp_get_wtime();
		time = (int)(end-start);
		printf("\none step (iteration) time:%.3f,  time taken = %d h %d min %d s\n", end - start1,  time/ 3600, \
			(time % 3600) / 60, time % 60);
	}

	//md
	end = omp_get_wtime();
	time = (int)(end-start);
	printf("\ntime taken = %d h %d min %d s\n",  time/ 3600, \
		(time % 3600) / 60, time % 60);
	PrintClose ();

	//femocs
	print_progress("\n> full run of Femocs", success == 0);

	return 0;
}

void print_progress(const string& message, const bool contition) {
	cout << message << ":  ";
	if (contition) cout << "passed" << endl;
	else cout << "failed" << endl;
}

void read_xyz(const string &file_name, double* x, double* y, double* z) {
	ifstream in_file(file_name, ios::in);
	require(in_file.is_open(), "Did not find a file " + file_name);

	int n_atoms, type, id;
	string elem, line;
	istringstream iss;

	getline(in_file, line); // Read number of atoms
	iss.clear();
	iss.str(line);
	iss >> n_atoms;

	getline(in_file, line); // Skip comments line

	id = -1;
	// keep storing values from the text file as long as data exists:
	while (++id < n_atoms && getline(in_file, line)) {
		iss.clear();
		iss.str(line);
		iss >> elem >> x[id] >> y[id] >> z[id] >> type;
	}
}

void read_ckx(const string &file_name, double* x, double* y, double* z) {
	ifstream in_file(file_name, ios::in);
	require(in_file.is_open(), "Did not find a file " + file_name);

	int n_atoms, type, id;
	string line;
	istringstream iss;

	getline(in_file, line); // Read number of atoms
	iss.clear();
	iss.str(line);
	iss >> n_atoms;

	getline(in_file, line);// Skip comments line

	id = -1;
	// keep storing values from the text file as long as data exists:
	while (++id < n_atoms && getline(in_file, line)) {
		iss.clear();
		iss.str(line);
		iss >> type >> x[id] >> y[id] >> z[id];
	}
}

void read_atoms(const string& file_name, double* x, double* y, double* z) {
	string file_type = get_file_type(file_name);

	if (file_type == "xyz")
	read_xyz(file_name, x, y, z);
	else if (file_type == "ckx")
	read_ckx(file_name, x, y, z);
	else
	require(false, "Unsupported file type: " + file_type);
}

void read_n_atoms(const string& file_name, int& n_atoms) {
	ifstream in_file(file_name, ios::in);
	require(in_file.is_open(), "Did not find a file " + file_name);

	string line;
	istringstream iss;

	getline(in_file, line); // Read number of atoms
	iss.clear();
	iss.str(line);
	iss >> n_atoms;
}

