#ifndef INIT_H
#define INIT_H

extern "C"{
#include "main.h"
}
#include "mlip.h"

void SetParams (int argc, char **argv);
void SetupJob(int argc, char **argv);
void AllocArrays ();
void AllocArraysMD();
void AllocArraysNebr ();
void InitCoords();
void InitCoordsFCC ();
void InitCoordsReadin ();
void InitVels();
void InitAccels();

#endif
