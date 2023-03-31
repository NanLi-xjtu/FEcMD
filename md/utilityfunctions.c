#include "utilityfunctions.h"

char *errorMsg[] = {"", "bond snapped", "read checkpoint data",
	"write checkpoint data", "copy buffer full", "empty event pool",
	"message buffer full", "outside region", "read snap data",
	"write snap data", "subdivision unfinished", "too many cells",
	"too many copied mols", "too many layers", "too many levels",
	"too many mols", "too many moved mols", "too many neighbors",
	"too many replicas"};

int randSeedP = 17;
int randSeed;

double RandR ()
{
	randSeedP = (randSeedP * IMUL + IADD) & MASK;
	return (randSeedP * SCALE);
}

void InitRand (int randSeedI)
{
	struct timeval tv;

	if(randSeedI != 0) randSeedP = randSeedI;
	else{
		gettimeofday (&tv, 0);
		randSeedP = tv.tv_usec;
	}
}


//2D
void VRand2D (VecR *p)
{
	double s;

	s = 2. * M_PI * RandR ();
	p->x = cos (s);
	p->y = sin (s);
}

//3D
void VRand3D (VecR *p)
{
	double s, x, y;

	s = 2.;
	while (s > 1.){
		x = 2. * RandR () - 1.;
		y = 2. * RandR () - 1.;
		s = Sqr (x) + Sqr (y);
	}
	p->z = 1. - 2. * s;
	s = 2. * sqrt (1. - s);
	p->x = s * x;
	p->y = s * y;
}


void ErrExit (int code)
{
	printf ("Error: %s\n", errorMsg[code]);
	exit (0);
}

FILE *ReadFile (char *filename)
{
	FILE *input;
	
	if ((input = fopen (filename, "r")) == NULL) {
		printf ("%s does not exist\n", filename);
		exit (1);
	}
	
	return (input);
}

FILE *WriteFile (char *filename)
{
	FILE *output;
	
	if ((output = fopen (filename, "a+")) == NULL) {
		printf ("openning %s failed\n", filename);
		exit (1);
	}
	
	return (output);
}

void Warning (char *line)
{
	FILE *output;
	char filename[128];

	sprintf (filename, "out/md/warning.out");
	output = WriteFile (filename);
	printf ("%s", line);
	fprintf (output, "%s", line);

	fclose (output);
}

void Error (char *line)
{
	FILE *output;
	char filename[128];

	sprintf (filename, "out/md/error.out");
	output = WriteFile (filename);
	printf ("%s", line);
	fprintf (output, "%s", line);

	fclose (output);
	exit (1);
}

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

  double *CreateDouble1 (int n)
  {
    double *array;

    bigint nbytes = ((bigint) sizeof(double)) * n;
    array = (double *) smalloc(nbytes,name);
    return array;
  }

  void DestroyDouble1 (double *array)
  {
    sfree(array);
    array = NULL;
  }

  int *CreateInt1 (int n)
  {
    int *array;

    bigint nbytes = ((bigint) sizeof(int)) * n;
    array = (int *) smalloc(nbytes,name);
    return array;
  }

  void DestroyInt1 (int *array)
  {
    sfree(array);
    array = NULL;
  }

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

  double **CreateDouble2 (int n1, int n2)
  {
    double **array;
    bigint nbytes = ((bigint) sizeof(double)) * n1*n2;
    double *data = (double *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double *)) * n1;
    array = (double **) smalloc(nbytes,name);

    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  void DestroyDouble2(double **array)
  {
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

  int **CreateInt2 (int n1, int n2)
  {
    int **array;
    bigint nbytes = ((bigint) sizeof(int)) * n1*n2;
    int *data = (int *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int *)) * n1;
    array = (int **) smalloc(nbytes,name);

    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  void DestroyInt2(int **array)
  {
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

/* ----------------------------------------------------------------------
   create a 3d array
------------------------------------------------------------------------- */

  double ***CreateDouble3 (int n1, int n2, int n3)
  {
    double ***array;
    bigint nbytes = ((bigint) sizeof(double)) * n1*n2*n3;
    double *data = (double *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double *)) * n1*n2;
    double **plane = (double **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double **)) * n1;
    array = (double ***) smalloc(nbytes,name);

    int i,j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }

  void DestroyDouble3 (double ***array)
  {
    if (array == NULL) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

  int ***CreateInt3 (int n1, int n2, int n3)
  {
    int ***array;
    bigint nbytes = ((bigint) sizeof(int)) * n1*n2*n3;
    int *data = (int *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int *)) * n1*n2;
    int **plane = (int **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int **)) * n1;
    array = (int ***) smalloc(nbytes,name);

    int i,j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }

  void DestroyInt3 (int ***array)
  {
    if (array == NULL) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

/* ----------------------------------------------------------------------
   create a 4d array
------------------------------------------------------------------------- */

  double ****CreateDouble4 (int n1, int n2, int n3, int n4)
  {
    double ****array;
    bigint nbytes = ((bigint) sizeof(double)) * n1*n2*n3*n4;
    double *data = (double *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double *)) * n1*n2*n3;
    double **cube = (double **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double **)) * n1*n2;
    double ***plane = (double ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double ***)) * n1;
    array = (double ****) smalloc(nbytes,name);

    int i,j,k;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 + j*n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 + j*n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }

  void DestroyDouble4 (double ****array)
  {
    if (array == NULL) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

  int ****CreateInt4 (int n1, int n2, int n3, int n4)
  {
    int ****array;
    bigint nbytes = ((bigint) sizeof(int)) * n1*n2*n3*n4;
    int *data = (int *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int *)) * n1*n2*n3;
    int **cube = (int **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int **)) * n1*n2;
    int ***plane = (int ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int ***)) * n1;
    array = (int ****) smalloc(nbytes,name);

    int i,j,k;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 + j*n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 + j*n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }

  void DestroyInt4 (int ****array)
  {
    if (array == NULL) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

/* ----------------------------------------------------------------------
   create a 5d array
------------------------------------------------------------------------- */

  double *****CreateDouble5 (int n1, int n2, int n3, int n4, int n5)
  {
    double *****array;
    bigint nbytes = ((bigint) sizeof(double)) * n1*n2*n3*n4*n5;
    double *data = (double *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double *)) * n1*n2*n3*n4;
    double **level4 = (double **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double **)) * n1*n2*n3;
    double ***level3 = (double ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double ***)) * n1*n2;
    double ****level2 = (double ****) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(double ****)) * n1;
    array = (double *****) smalloc(nbytes,name);

    int i,j,k,l;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &level2[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3;
        level2[m1] = &level3[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3 + k;
          m2 = ((bigint) i) * n2*n3*n4 +
            ((bigint) j) * n3*n4 + ((bigint) k) * n4;
          level3[m1] = &level4[m2];
          for (l = 0; l < n4; l++) {
            m1 = ((bigint) i) * n2*n3*n4 +
              ((bigint) j) * n3*n4 + ((bigint) k) * n4 + l;
            level4[m1] = &data[n];
            n += n5;
          }
        }
      }
    }
    return array;
  }

  void DestroyDouble5 (double *****array)
  {
    if (array == NULL) return;
    sfree(array[0][0][0][0]);
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

  int *****CreateInt5 (int n1, int n2, int n3, int n4, int n5)
  {
    int *****array;
    bigint nbytes = ((bigint) sizeof(int)) * n1*n2*n3*n4*n5;
    int *data = (int *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int *)) * n1*n2*n3*n4;
    int **level4 = (int **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int **)) * n1*n2*n3;
    int ***level3 = (int ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int ***)) * n1*n2;
    int ****level2 = (int ****) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(int ****)) * n1;
    array = (int *****) smalloc(nbytes,name);

    int i,j,k,l;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &level2[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3;
        level2[m1] = &level3[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3 + k;
          m2 = ((bigint) i) * n2*n3*n4 +
            ((bigint) j) * n3*n4 + ((bigint) k) * n4;
          level3[m1] = &level4[m2];
          for (l = 0; l < n4; l++) {
            m1 = ((bigint) i) * n2*n3*n4 +
              ((bigint) j) * n3*n4 + ((bigint) k) * n4 + l;
            level4[m1] = &data[n];
            n += n5;
          }
        }
      }
    }
    return array;
  }

  void DestroyInt5 (int *****array)
  {
    if (array == NULL) return;
    sfree(array[0][0][0][0]);
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = NULL;
  }

void ZeroDouble1 (double *a, int n1)
{
	int i;

	for (i = 0; i < n1; i ++)
			a[i] = 0.;
}
void ZeroDouble2 (double **a, int n1, int n2)
{
	int i, j;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			a[i][j] = 0.;
}
void ZeroDouble3 (double ***a, int n1, int n2, int n3)
{
	int i, j, k;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				a[i][j][k] = 0.;
}
void ZeroDouble4 (double ****a, int n1, int n2, int n3, int n4)
{
	int i, j, k, l;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				for (l = 0; l < n4; l ++)
					a[i][j][k][l] = 0.;
}
void ZeroDouble5 (double *****a, int n1, int n2, int n3, int n4, int n5)
{
	int i, j, k, l, m;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				for (l = 0; l < n4; l ++)
					for (m = 0; m < n5; m ++)
						a[i][j][k][l][m] = 0.;
}
void ZeroInt1 (int *a, int n1)
{
	int i;

	for (i = 0; i < n1; i ++)
			a[i] = 0;
}
void ZeroInt2 (int **a, int n1, int n2)
{
	int i, j;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			a[i][j] = 0;
}
void ZeroInt3 (int ***a, int n1, int n2, int n3)
{
	int i, j, k;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				a[i][j][k] = 0;
}
void ZeroInt4 (int ****a, int n1, int n2, int n3, int n4)
{
	int i, j, k, l;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				for (l = 0; l < n4; l ++)
					a[i][j][k][l] = 0;
}
void ZeroInt5 (int *****a, int n1, int n2, int n3, int n4, int n5)
{
	int i, j, k, l, m;

	for (i = 0; i < n1; i ++)
		for (j = 0; j < n2; j ++)
			for (k = 0; k < n3; k ++)
				for (l = 0; l < n4; l ++)
					for (m = 0; m < n5; m ++)
						a[i][j][k][l][m] = 0;
}
