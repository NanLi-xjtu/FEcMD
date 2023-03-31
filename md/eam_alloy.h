#ifndef EAM_ALLOY_H
#define EAM_ALLOY_H

#include "main.h"

typedef struct {
  char **elements, **structure;
  int nelements,nrho,nr;
  double drho,dr,cut;
  double *mass, *lattice;
  double **frho,**rhor,***z2r;
} Setfl;

typedef struct {
    char **elements, **structure;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass, *lattice;
    double **frho,***rhor,***z2r;
} Fs;

extern Setfl *setfl;
extern Fs *fs;

void coeff_alloy(char *filename);
void read_file_alloy(char *filename);
void file2array_alloy ();
void init_style_alloy();
void printfeam_alloy ();
void printfspline_alloy (double ***xxxspline);
void ComputeForcesEamPoten_alloy ();
void coeff_eamfs();
void read_file_eamfs(char *filename);
void init_style_eamfs ();

#endif
