#ifndef MAIN_H
#define MAIN_H

#include "Femocs.h"
#include "Macros.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

extern "C"
{
#include "define.h"
#include "utilityfunctions.h"
#include "namelist.h"
#include "measurements.h"
#include "forces.h"
#include "heattransport.h"
#include "vscale.h"
#include "structure.h"
#include "eamforce.h"
#include "eam_alloy.h"
#include "integral.h"
#include "nebrList.h"
#include "tipMaker.h"
#include "npt.h"
#include "boundary.h"
#include "interface.h"
#include "snap.h"
}
#include "heatdiffusion.h"
#include "singlestep.h"
#include "mlip.h"
#include "init.h"

void print_progress(const string& message, const bool contition);
void read_xyz(const string &file_name, double* x, double* y, double* z);
void read_ckx(const string &file_name, double* x, double* y, double* z);
void read_atoms(const string& file_name, double* x, double* y, double* z);
void read_n_atoms(const string& file_name, int& n_atoms);

#endif

