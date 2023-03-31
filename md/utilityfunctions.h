#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include "main.h"

#define IADD	453806245
#define IMUL	314159269
#define MASK	2147483647
#define SCALE	0.4656612873e-9

enum {ERR_NONE, ERR_BOND_SNAPPED, ERR_CHECKPT_READ, ERR_CHECKPT_WRITE,
	ERR_COPY_BUFF_FULL, ERR_EMPTY_ECPOOL, ERR_MSG_BUFF_FULL,
	ERR_OUTSIDE_REGION, ERR_SNAP_READ, ERR_SNAP_WRITE,
	ERR_SUBDIV_UNFIN, ERR_TOO_MANY_CELLS, ERR_TOO_MANY_COPIES, 
	ERR_TOO_MANY_LAYERS, ERR_TOO_MANY_LEVELS, ERR_TOO_MANY_MOLS,
	ERR_TOO_MANY_MOVES, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_REPLICAS};

#define smalloc(nbytes, name) malloc(nbytes)
#define bigint long
#define sfree free
#define free2(a) \
	free(a[0]); \
	free(a)
#define free3(a) \
	free(a[0][0]); \
	free(a[0]); \
	free(a)
#define free4(a) \
	free(a[0][0][0]); \
	free(a[0][0]); \
	free(a[0]); \
	free(a);
#define free5(a) \
	free(a[0][0][0][0]);\
	free(a[0][0][0]); \
	free(a[0][0]); \
	free(a[0]); \
	free(a)

extern int randSeed;

double RandR ();
void InitRand (int randSeedI);
void VRand2D (VecR *p);
void VRand3D (VecR *p);
void ErrExit (int code);
double *CreateDouble1 (int n);
void DestroyDouble1(double *array);
double **CreateDouble2 (int n1, int n2);
void DestroyDouble2(double **array);
double ***CreateDouble3 (int n1, int n2, int n3);
void DestroyDouble3 (double ***array);
double ****CreateDouble4 (int n1, int n2, int n3, int n4);
void DestroyDouble4 (double ****array);
double *****CreateDouble5 (int n1, int n2, int n3, int n4, int n5);
void DestroyDouble5 (double *****array);
int *CreateInt1 (int n);
void DestroyInt1(int *array);
int **CreateInt2 (int n1, int n2);
void DestroyInt2(int **array);
int ***CreateInt3 (int n1, int n2, int n3);
void DestroyInt3 (int ***array);
int ****CreateInt4 (int n1, int n2, int n3, int n4);
void DestroyInt4 (int ****array);
int *****CreateInt5 (int n1, int n2, int n3, int n4, int n5);
void DestroyInt5 (int *****array);
void ZeroDouble1 (double *a, int n1);
void ZeroDouble2 (double **a, int n1, int n2);
void ZeroDouble3 (double ***a, int n1, int n2, int n3);
void ZeroDouble4 (double ****a, int n1, int n2, int n3, int n4);
void ZeroDouble5 (double *****a, int n1, int n2, int n3, int n4, int n5);
void ZeroInt1 (int *a, int n1);
void ZeroInt2 (int **a, int n1, int n2);
void ZeroInt3 (int ***a, int n1, int n2, int n3);
void ZeroInt4 (int ****a, int n1, int n2, int n3, int n4);
void ZeroInt5 (int *****a, int n1, int n2, int n3, int n4, int n5);

FILE *ReadFile (char filename[128]);
FILE *WriteFile (char filename[128]);
void Warning (char *line);
void Error (char *line);

#endif
