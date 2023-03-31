extern "C"{
#include "main.h"
}

#include "init.h"
#include "mlip.h"
#include "singlestep.h"
#include "heatdiffusion.h"

using namespace std;

void CheckParallel ()
{
	int coreNum = omp_get_num_procs ();

	printf("Core Num is %d \n", coreNum);
}

int main (int argc, char **argv)
{
	double start, end;
	double T;

	system ("rm -rf out");
	system ("mkdir out");

	sprintf (CONTCAR_file, "in/CONTCAR");
	sprintf (mdPamsFilename, "in/md.in");
	sprintf (mdXYZFilename, "in/mdlat.in.xyz");

	CheckParallel ();
	PrintOpen ();
//	GetNameList (argc, argv);
	GetParaValue ();

	switch (argc)
	{
		case 1:
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
			} else if (strcmp (argv[1], "Maxwell") == 0) {
				sprintf (force_sort, "Maxwell");
				break;
			} else if (strcmp (argv[1], "heattransport") == 0) {
				sprintf (nonq_mode, "heattransport");
				break;
			} else {
				printf ("error: input format\n");
				printf ("./md\n      tip\n      supercell\n      extension\n      prolate\n      elliptic\n");
				printf ("      Maxwell\n      heattransport\n");
				exit (1);
			}
			return 0;
		default:
			printf ("error: input format\n");
			exit (1);
	}

	SetParams (argc, argv);
	SetupJob (argc, argv);
//	TempDistriInit();
	moreCycles = 1;

//	PrintSummaryTitle ();
	PRINTSUMMARYTITLE ();
	PrintFemocsin_xyz ();
	start = omp_get_wtime();
	T = TEMPERATURE_K;
	while (moreCycles){

//		TempDistri ();
		SingleStep ();

		if (stepCount >= stepLimit) moreCycles = 0;
	}
	end = omp_get_wtime();
	printf("\ntime taken = %d h %d min %d s\n", (int)((end-start)/3600.), (int)(((int)(end-start)%3600)/60.), (int)(end-start)%60);

	printf("\nnumbers of atoms removed is %d\n", atomFlyCount);

//	TempPrintClose ();
	PrintClose ();
}
