#ifndef EAMFORCE_H
#define EAMFORCE_H

#include "main.h"

#define MAXLINE 1024

// potentials as file data

typedef struct {
  char *file, structure[10];
  int nrho,nr;
  int id; //--atomic number
  double drho,dr,cut,mass,lattice;
  double *frho,*rhor,*zr;
} Funcfl;

extern Funcfl *funcfl;

extern double cutmax;
extern double ***rhor_spline,***frho_spline,***z2r_spline, ***phi_spline;
extern char eamFilename[128];
extern int allocated;
extern int ntypes;
extern char force_type[10];

extern int *map;                   // which element each atom type maps to
extern int nmax;                   // allocated size of per-atom arrays
extern double cutforcesq;
extern double **scale;

extern int nrho,nr;
extern int nfrho,nrhor,nz2r;
extern double **frho,**rhor,**z2r;
extern int *type2frho,**type2rhor,**type2z2r;

// potentials in spline form used for force computation

extern double dr,rdr,drho,rdrho,rhomax;
extern double ***rhor_spline,***frho_spline,***z2r_spline, ***phi_spline;
extern double *spline;

// per-atom arrays

extern double *rho,*fp;

//total energy
extern double frho_tot, pairpoten_tot, uSum;

void allocate();
void potential_date(FILE *fp, const char *name);
FILE * open_potential(const char *name);
void grab(FILE *fptr, int n, double *list);
void read_file(char *filename);
void Coeff();
void file2array();
void interpolate(int n, double delta, double *f, double **spline);
void array2spline();
void init_style();
double init_one(int i, int j);
void printFuncfl (int i);
void printfeam(int i);
void printfspline (double ***xxxspline);
void ComputeForcesEamPoten ();

#endif
