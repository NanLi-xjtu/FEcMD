#include "measurements.h"

double *histVel, rangeVel;
int countVel, limitVel, sizeHistVel, stepVel;
double hFunction;
int flag_vDis = 0;

//the following variables are dimensionless
void EvalProps ()
{
	double vv;
	int n;
	double vvMax;
	double zmin, zmax;
	double EkSum;

	VZero (vSum);
	VZero (pSum);
	vvSum = 0.;
	vvMax = 0.;
	EkSum = 0.;
	tipHeight.val = 0.;
	zmin = 1.e100;
	zmax = -1.e100;
	DO_MOL{
		VVAdd (vSum, mol[n].rv);
		VVSAdd (pSum, mol[n].mass, mol[n].rv);
		vv = VLenSq (mol[n].rv);
		vvSum += vv;
		EkSum += 0.5 * mol[n].mass * vv;
		vvMax = Max (vvMax, vv);
		if (mol[n].anode == 0) {
			zmin = Min (zmin, mol[n].r.z);
			zmax = Max (zmax, mol[n].r.z);
		}
	}

	tipHeight.val = zmax - zmin;
	kinEnergy.val = EkSum / nMol;
	
	pressure.val = density * (2. * EkSum + virSum) / (nMol * NDIM);
	totEnergy.val = kinEnergy.val + uSum / nMol;
	if (strcmp (ensemble, "NPT") == 0) {
		totEnergy.val += (0.5 * (massS * Sqr (varSv) + massV * Sqr (varVv)) / 
//				 Sqr (varS) + extPressure * varV) /
				 Sqr (varS) + extPressure * nMol / density) /
				 nMol + 3. * temperature * log (varS);
//		pressure.val = (2. * EkSum + virSum) / (3. * varV);
	}

	KINENERGY_eV.val = kinEnergy.val * (eUnit / eleChar);
	TEMPAVG_K.val = ((2. / NDIM) * kinEnergy.val) * TUnit;
	TOTENERGY_eV.val = totEnergy.val * (eUnit / eleChar);
	PRESSURE_GPa.val = pressure.val * (PUnit / 1.e9);
	TIPHEIGHT_A.val = tipHeight.val * (lUnit * 1.e10);
	
	//heat transport
	thermalCond.val = 0.5 * enTransSum / (deltaT * region.x * region.y * ((wallTempHi - wallTempLo) / region.z));

	//neighbor list
	dispHi += sqrt (vvMax) * deltaT;//the maximum movement of the atoms
	if (dispHi > 0.5 * rNebrShell) nebrNow = 1;
}

void AccumProps(int icode)
{
	if (icode == 0){
		PropZero (totEnergy);
		PropZero (kinEnergy);
		PropZero (pressure);
		PropZero (tipHeight);
		PropZero (TOTENERGY_eV);
		PropZero (KINENERGY_eV);
		PropZero (TEMPAVG_K);
		PropZero (PRESSURE_GPa);
		PropZero (TIPHEIGHT_A);
		
		//heat transport
		PropZero (thermalCond);
	
	}else if (icode == 1){
		PropAccum (totEnergy);
		PropAccum (kinEnergy);
		PropAccum (pressure);
		PropAccum (tipHeight);
		PropAccum (TOTENERGY_eV);
		PropAccum (KINENERGY_eV);
		PropAccum (TEMPAVG_K);
		PropAccum (PRESSURE_GPa);
		PropAccum (TIPHEIGHT_A);

		//heat transport
		PropAccum (thermalCond);

	}else if (icode == 2){
		PropAvg (totEnergy, stepAvg);
		PropAvg (kinEnergy, stepAvg);
		PropAvg (pressure, stepAvg);
		PropAvg (tipHeight, stepAvg);
		PropAvg (TOTENERGY_eV, stepAvg);
		PropAvg (KINENERGY_eV, stepAvg);
		PropAvg (TEMPAVG_K, stepAvg);
		PropAvg (PRESSURE_GPa, stepAvg);
		PropAvg (TIPHEIGHT_A, stepAvg);

		//heat transport
		PropAvg (thermalCond, stepAvg);
	}
}

void PrintVelDist ()
{
	double vBin;
	int n;
	FILE *fp, *fp1;
	char filename[128];

	sprintf (filename, "out/md/v_distribution.dat");
	fp = WriteFile (filename);
	sprintf (filename, "out/md/H_function.dat");
	fp1 = WriteFile (filename);
	fprintf (fp, "vdist(A/ps) (%fps)\n", timeNow*tUnit*1.e12);
	for (n = 0; n < sizeHistVel; n ++) {
		vBin = (n + 0.5) * rangeVel / sizeHistVel;
		fprintf (fp, "%f %f\n", vBin*vUnit*1.e-2, histVel[n]);
	}
	if (flag_vDis == 0) fprintf (fp1, "t(ps) h-function\n");
	fprintf (fp1, "%f %f\n", timeNow*tUnit*1.e12, hFunction);

	flag_vDis = 1;

	fclose (fp);
	fclose (fp1);
}

void EvalVelDist ()
{
	double deltaV, histSum;
	int j, n;

	if (countVel == 0) {
		for (j = 0; j < sizeHistVel; j ++){
			histVel[j] = 0.;
		}
	}
	deltaV = rangeVel / sizeHistVel;
	DO_MOL{
		j = VLen (mol[n].rv) / deltaV;
		++ histVel[Min (j, sizeHistVel - 1)];
	}
	++ countVel;
	if (countVel == limitVel) {
		histSum = 0.;
		for (j = 0; j < sizeHistVel; j ++){
			histSum += histVel[j];
		}
		for (j = 0; j < sizeHistVel; j ++){
			histVel[j] /= histSum;
		}
		//Boltzmann H-function
		hFunction = 0.;
		for (j = 0; j < sizeHistVel; j ++) {
			if (histVel[j] > 0.) hFunction += histVel[j] * log (histVel[j] / 
							  pow(((j + 0.5) * deltaV), (NDIM-1.)));
		}
		PrintVelDist ();
		countVel = 0;
	}
}

