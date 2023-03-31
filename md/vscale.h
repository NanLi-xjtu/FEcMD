#ifndef VSCALE_H
#define VSCALE_H

#include "main.h"

#define DO_TEMPCELL				\
	for (i = 0; i < sizeTempGrid; i ++)

extern int stepAdjustTemp, stepInitlzTemp;
extern double kinEnInitSum, Trate;

void AdjustTemp (double tempAdjust_K);
void AdjustInitTemp ();
void CalcuTempCell ();
void PrintTempCell ();
void Scaling ();

#endif
