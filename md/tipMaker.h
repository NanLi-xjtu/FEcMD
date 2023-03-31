#ifndef TIPMAKER_H
#define TIPMAKER_H

#include "main.h"

typedef struct {
	int Natoms;
	VecR *r;
	char name[8];
} Elem;

typedef struct {
	VecR r;
	char element[8];
	int delete_flag;
} Atom;

extern VecI cellSize;
extern char CONTCAR_file[128], tipFilename[128], cellFilename[128], cellOrder[128];
extern double tip_r, tip_r1, tip_h, tip_theta;

void CreatSuperCell ();
void TipMaker ();
void MushroomTipMaker ();
void SuperCellMaker ();
void ExtensionMaker ();
void ProlateSpheroidalMaker ();
void EllipticMaker ();
void ConeTipMaker ();
void CylinderTipMaker ();

#endif
