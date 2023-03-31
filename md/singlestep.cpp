#include "singlestep.h"
#include "mlip.h"

int stepCountOld = 0, stepMaxwell = 0;
double start_md, end_md;

void SingleStep ()
{
	int n;

	if (stepCount == 0) start_md = omp_get_wtime ();

	++ stepCount;
	timeNow = stepCount * deltaT;
	TIMENOW_ps = stepCount * DELTAT_fs * 1.e-3;

	//move the atoms
	if (strcmp(ensemble, "NPT") == 0) {
		PredictorStep ();
		PredictorStepPT (); //NPT
	} else LeapfrogStep (1);

	/*******************************************
	** ApplyBoundaryCond (); --periodic boundary
	** ApplyBoundaryCondFlow (); --viscous flow
	** ApplyBoundaryCondHeat (); --heat transport
	** ApplyBoundaryCondSurface (); --surface
	** AtomFlySurface (); --delete the atom flying out of the surface
	*******************************************/
	if (stepCount >= stepEquil && strcmp(nonq_mode, "heattransport") == 0) ApplyBoundaryCondHeat ();
	else if (strcmp(boundaryCond, "n") == 0) {
		if (ligancy_flag == 1 && boundary_fly == 1) AtomFlySurface ();
		AtomFlyBox ();
	} else if (strcmp(boundaryCond, "p") == 0) ApplyBoundaryCond ();

	if (strcmp(ensemble, "NPT") == 0) {
		UpdateCellSize (); //NPT
		UnScaleCoords (); //NPT
	}

	//neighbor list
	if (strcmp(interact_method, "nebr") == 0 && (nebrNow || atomFlyFlag)) {
		nebrNow = 0;
		atomFlyFlag = 0;
		dispHi = 0.;

		//elete the atom flying out of the Box
//		AtomFlyBox ();
		/************************************************************
		** BuildNebrList ();
		** BuildNebrListHeat (); --heat transport or viscous flow
		** BuildNebrListCommon (); --use common method to search all pairs
		*************************************************************/
		if (stepCount >= stepEquil && strcmp(nonq_mode, "heattransport") == 0) BuildNebrListHeat ();
		else if (nMol < 200) BuildNebrListCommon ();
		else BuildNebrList ();

		stepNebr = stepCount - stepCountOld;
		stepCountOld = stepCount;
	} else if (strcmp(interact_method, "cell") == 0) BuildCellList ();

	/******************************************************
	** the different approaches to computing interactions.
	*******************************************************
	** ComputeForces (); --all pairs, lj
	** ComputeForcesCell (); --cell subdivision, lj
	** ComputeForcesNebr (); --neighbor/cell list, lj
	** ComputeForcesEam (); --EAM potential, simple
	** ComputeForcesEamPoten (); -- EAM potential files
	** ComputeForcesEamPoten_alloy (); -- alloy EAM potential files
	** ComputeForcesMTP (); --MTP machine learning potential
	******************************************************/
	if (strcmp (force_type, "metal") == 0) ComputeForcesEamPoten ();
	else if (strcmp (force_type, "alloy") == 0 || strcmp (force_type, "eamfs") == 0) {
		ComputeForcesEamPoten_alloy ();
	} else if (strcmp (force_type, "lj") == 0) {
		if (strcmp(interact_method, "allpairs") == 0) ComputeForces ();
		else ComputeForcesNebr ();
	} else if (strcmp (force_type, "snap") == 0) {
		PairSNAP_compute_regular(1, 1);
	} else if (strcmp (force_type, "MTP") == 0) {
		ComputeForcesMTP ();
	} else {
		printf ("error(singlestep.c): force_type format\n");
		exit (1);
	}

	if (strcmp (force_sort, "Maxwell") == 0) {
		Maxwell_stress += Maxwell_rate / 1000. * (DELTAT_fs); //--GPa/ps
		if (Maxwell_stress > Maxwell_max && fabs(Maxwell_max) > 1.e-100) Maxwell_stress = Maxwell_max;
		if (fabs(Maxwell_rate) < 1.e-100) Maxwell_stress = Maxwell_begin;
		if (stepCount % stepAvg == 0) printf ("\nMaxwell stress = %f Gpa\n", Maxwell_stress);
		MaxwellStress (Maxwell_stress); //--GPa
	}

	//surface atoms are disturbed by the Coulomb force and Lolorz force, or Maxwell stress
	#pragma omp parallel for num_threads(Nthreads)
	DO_MOL VVAdd (mol[n].f, mol[n].F);

	//compute accelerated velocity
	#pragma omp parallel for num_threads(Nthreads)
	DO_MOL VSCopy (mol[n].ra, 1. / mol[n].mass, mol[n].f); //--dimensionless

	//viscous flow
//	ComputeExternalForce ();

	if (strcmp(ensemble, "NPT") == 0) {
		ComputeDerivsPT (); //NPT
		CorrectorStep ();
		CorrectorStepPT (); //NPT
	} else LeapfrogStep (2);

	if (strcmp(nonq_mode, "heattransport") != 0) {
		if (strcmp(boundaryCond, "n") == 0 && ligancy_flag == 1) {
			if (ligancy_flag == 1 && boundary_fly == 1) AtomFlySurface ();
			AtomFlyBox ();
		} else if (strcmp(boundaryCond, "p") == 0) ApplyBoundaryCond ();
	}

	if (strcmp(ensemble, "NPT") == 0) {
		UnScaleCoords ();
		UnScaleVels ();
	}

	//calculate and print summary
	EvalProps ();
	AccumProps (1);
	if (stepCount % stepAvg == 0){
		AccumProps (2);
		PRINTSUMMARY ();
		AccumProps (0);
	}

	//ligancy analysis
	if (stepCount % stepLigancy == 0) EvalLigancy ();

	//mark surface atoms
	if (stepCount % stepSurf == 0) SurfaceAtoms ();

	if (stepCount % stepMovie == 0) {
		PrintMovie ();//"md.movie"
	}

	//initial temperature adjustment
	if ((strcmp(ensemble, "NVT") == 0 || strcmp(ensemble, "nonEquil") == 0) \
	    && stepCount < stepEquil && stepCount % stepInitlzTemp == 0) {
		AdjustInitTemp ();
	}

	//adjust temperature
	TEMPERATURE_K += Trate; //--K
	if (strcmp(ensemble, "NVT") == 0 && stepCount > stepEquil \
	    && (stepCount - stepEquil) % stepAdjustTemp == 0) {
		AdjustTemp (TEMPERATURE_K);
	}

	//velocity distribution
	if (stepCount % stepVel == 0) EvalVelDist ();


	//heat transport measurement and print profile
	if ((stepCount - stepEquil) % stepGrid == 0){
		++ countGrid;
		GridAverage (1);
		if (countGrid % limitGrid == 0){
			GridAverage (2);
			EvalProfile ();
			PrintProfile ();
			GridAverage (0);
		}
	} 

/*
	//RDF
	if (stepCount >= stepEquil && (stepCount - stepEquil) % stepRdf == 0){
		EvalRdf();
	}

*/

/*
	//scaling
	CalcuTempCell ();
	if (stepCount > stepEquil && (stepCount - stepEquil) % stepTempGrid == 0){
		Scaling ();
	}
	if (stepCount > stepEquil && (stepCount - stepEquil) % (stepGrid * limitGrid) == 0){
		PrintTempCell();
	}
*/

	if (strcmp(ensemble, "NPT") == 0) {
		ScaleCoords ();
		ScaleVels ();
	}
}

