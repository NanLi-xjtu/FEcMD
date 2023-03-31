#include "vscale.h"

//tempearture adjust
int stepAdjustTemp, stepInitlzTemp;
double kinEnInitSum, Trate;

void AdjustTemp (double tempAdjust_K)
{
	double vFac, tempAdjust, k, EkSum;
	int n;

	tempAdjust = tempAdjust_K / TUnit; //dimensionless
	k = sqrt (NDIM * (1. - 1. / nMol) * tempAdjust);

	EkSum = 0.;
	DO_MOL EkSum += 0.5 * mol[n].mass * VLenSq (mol[n].rv);
	vFac = k / sqrt (2. * EkSum / nMol);
	DO_MOL VScale (mol[n].rv, vFac);
}

void AdjustInitTemp ()
{
	double vFac;
	int n;

	kinEnInitSum += kinEnergy.val;
	if (stepCount % stepInitlzTemp == 0){
		kinEnInitSum /= stepInitlzTemp;
		vFac = velMag / sqrt (2. * kinEnInitSum);
		DO_MOL {
			if (mol[n].pedest != 1) VScale (mol[n].rv, vFac); 
		}
		kinEnInitSum = 0.;
	}
}

void CalcuTempCell ()
{
	int n, i;

	DO_TEMPCELL{
		tempCell[i].vvSum = 0;
		tempCell[i].nmol = 0;
	} 
	DO_MOL{
//		if (tipHeight.sum == 0){
			i = (int)(mol[n].r.z / region.z * (double)sizeTempGrid);
//		}else{
//			i = (int)((mol[n].r.z / tipHeight.sum + 0.5) * (double)sizeTempGrid);
//		}
		if (i >= 0 && i < sizeTempGrid){
			tempCell[i].vvSum += VLenSq(mol[n].rv);
			tempCell[i].nmol ++;
		}
	}
}

void PrintTempCell()
{
	fprintf (tempcell, "   i   T(K)   EkSum(eV) nmol\n");

	int i;

	DO_TEMPCELL{
		fprintf (tempcell, "%4d %.3f %f %d\n", i, tempCell[i].temp * TUnit, 0.5 * tempCell[i].vvSum * (eUnit / eleChar), tempCell[i].nmol);
	}
}

void Scaling ()
{
	double vFac, vNew;
	int n, i;

	stepSca ++;
	DO_MOL{
//		if (tipHeight.sum == 0){
			i = (int)((mol[n].r.z / region.z + 0.5) * (double)sizeTempGrid);
//		}else{
//			i = (int)((mol[n].r.z / tipHeight.sum + 0.5) * (double)sizeTempGrid);
//		}
		if (i >= 0 && i < sizeTempGrid){
			vNew = sqrt (NDIM * (1. - 1. / nMol) * tempCell[i].temp);
			vFac = vNew / sqrt (tempCell[i].vvSum / tempCell[i].nmol);
			VScale (mol[n].rv, vFac);
		}
	}
	
}
