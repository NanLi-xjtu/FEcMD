#ifndef SNAP_H
#define SNAP_H

#include "main.h"

typedef struct {
  double* bvec, ** dbvec;
  double** rij;
  int* inside;
  double* wj;
  double* rcutij;
} SNA;

typedef struct {
  int j1, j2, j;
} SNA_LOOPINDICES;

#define NEIGHMASK 0x3FFFFFFF
#define MY_PI M_PI

extern double sna_rcutmax;               // max cutoff for all elements
extern char snapFilename[128];

void PairSNAP_read_files(char *coefffilename, char *paramfilename);
void PairSNAP_coeff();
void PairSNAP_init_style();
void PairSNAP_compute_regular(int eflag, int vflag);
void ComputeSNAAtom ();
void SNA_init ();

#endif
