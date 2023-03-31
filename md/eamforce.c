#include "eamforce.h"

FILE *spline_file;
FILE *frhoeam_file;
FILE *zream_file;
FILE *rhoream_file;

int allocated = 0;
int ntypes;
char force_type[10];

// public variables so USER-ATC package can access them

double cutmax;

// potentials as array data

int nfuncfl = 0;
Funcfl *funcfl;
Fs *fs;

int nrho,nr;
int nfrho,nrhor,nz2r;
double **frho,**rhor,**z2r;
int *type2frho,**type2rhor,**type2z2r;

// potentials in spline form used for force computation

double dr,rdr,drho,rdrho,rhomax;
double ***rhor_spline,***frho_spline,***z2r_spline, ***phi_spline;
double *spline;

int nmax;                   // allocated size of per-atom arrays
double cutforcesq;
double **scale;

// per-atom arrays

double *rho,*fp;

int *map;                   // which element each atom type maps to

//total energy
double frho_tot, pairpoten_tot, uSum;

char eamFilename[128];

void allocate()
{
  int i, j, k;
  int n;

  n = ntypes = 1;
  allocated = 1;
  AllocMem (map, n+1, int);
  for (int i = 1; i <= n; i++) map[i] = -1;
  
  AllocMem (type2frho, n+1, int);
  AllocMem2 (type2rhor, n+1, n+1, int);
  AllocMem2 (type2z2r, n+1, n+1, int);
  AllocMem2 (scale, n+1, n+1, double);
}

void potential_date(FILE *fp, const char *name)
{
  char line[MAXLINE];
  char *ptr = fgets(line,MAXLINE,fp);
  if (ptr == NULL) return;

  char *word;
  word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"DATE:") == 0) {
      word = strtok(NULL," \t\n\r\f");
      if (word == NULL) return;
      printf("Reading potential file %s with DATE: %s\n",name,word);
      return;
    }
    word = strtok(NULL," \t\n\r\f");
  }
}

FILE * open_potential(const char *name)
{
  FILE *fp;

  if (name == NULL) return NULL;

  // attempt to open file directly
  // if successful, return ptr

  fp = fopen(name,"r");
  if (fp) {
    potential_date(fp,name);
    rewind(fp);
    return fp;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl-1];

  FILE *fptr;
  char line[MAXLINE];

    fptr = open_potential(filename);
    if (fptr == NULL) {
      printf("Cannot open EAM potential file %s\n",filename);
    }

  int tmp,nwords;
  fgets(line,MAXLINE,fptr);
  printf("\n%s\n", line);
  fgets(line,MAXLINE,fptr);
//  printf("\n%s\n", line);
  sscanf(line,"%d %lg %lg %s",&file->id,&file->mass,&file->lattice,file->structure);
  fgets(line,MAXLINE,fptr);
  nwords = sscanf(line,"%d %lg %d %lg %lg",
         &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);

  if ((nwords != 5) || (file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
    printf("Invalid EAM potential file");

  AllocMem(file->frho, file->nrho + 1, double);
  AllocMem(file->rhor, file->nr + 1, double);
  AllocMem(file->zr, file->nr + 1, double);

  grab(fptr,file->nrho,&file->frho[1]);
  grab(fptr,file->nr,&file->zr[1]);
  grab(fptr,file->nr,&file->rhor[1]);

  fclose(fptr);
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void Coeff ()
{
  char arg[2][128];
  sprintf (arg[1], "%s", eamFilename);

  if (!allocated) allocate();

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[1],funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    ReAllocMem(funcfl, nfuncfl, Funcfl);
    read_file(arg[1]);
    int n = strlen(arg[1]) + 1;
    AllocMem (funcfl[ifuncfl].file, n, char);
    strcpy (funcfl[ifuncfl].file, arg[1]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i = j
  int i = 1;
  int j = 1;
  map[i] = ifuncfl;
  scale[i][j] = 1.0;

}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void file2array()
{
  int i,j,k,m,n;
  double sixth = 1.0/6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

    Funcfl *file = &funcfl[i];
  

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1) * file->dr);
    rhomax = MAX(rhomax,(file->nrho-1) * file->drho);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = (int) (rmax/dr + 0.5);
  nrho = (int) (rhomax/drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  AllocMem2(frho, nfrho, nrho+1, double);

  // interpolate each file's frho to a single grid and cutoff

  double r,p,cof1,cof2,cof3,cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/file->drho + 1.0;
      k = (int) (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      frho[n][m] = cof1*file->frho[k-1] + cof2*file->frho[k] +
        cof3*file->frho[k+1] + cof4*file->frho[k+2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;
  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  AllocMem2(rhor, nrhor, nr + 1, double);

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/file->dr + 1.0;
      k = (int) (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      rhor[n][m] = cof1*file->rhor[k-1] + cof2*file->rhor[k] +
        cof3*file->rhor[k+1] + cof4*file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;
  AllocMem2 (z2r, nz2r, nr + 1, double);

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  double zri,zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
        r = (m-1)*dr;

        p = r/ifile->dr + 1.0;
        k = (int) (p);
        k = MIN(k,ifile->nr-2);
        k = MAX(k,2);
        p -= k;
        p = MIN(p,2.0);
        cof1 = -sixth*p*(p-1.0)*(p-2.0);
        cof2 = 0.5*(p*p-1.0)*(p-2.0);
        cof3 = -0.5*p*(p+1.0)*(p-2.0);
        cof4 = sixth*p*(p*p-1.0);
        zri = cof1*ifile->zr[k-1] + cof2*ifile->zr[k] +
          cof3*ifile->zr[k+1] + cof4*ifile->zr[k+2];

        p = r/jfile->dr + 1.0;
        k = (int) (p);
        k = MIN(k,jfile->nr-2);
        k = MAX(k,2);
        p -= k;
        p = MIN(p,2.0);
        cof1 = -sixth*p*(p-1.0)*(p-2.0);
        cof2 = 0.5*(p*p-1.0)*(p-2.0);
        cof3 = -0.5*p*(p+1.0)*(p-2.0);
        cof4 = sixth*p*(p*p-1.0);
        zrj = cof1*jfile->zr[k-1] + cof2*jfile->zr[k] +
          cof3*jfile->zr[k+1] + cof4*jfile->zr[k+2];

        z2r[n][m] = 27.2*0.529 * zri*zrj;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

void interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

void array2spline ()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  int j, k;
  AllocMem3 (frho_spline, nfrho, (nrho + 1), 7, double);
  AllocMem3 (rhor_spline, nrhor, (nr + 1), 7, double);
  AllocMem3 (z2r_spline, nz2r, (nr + 1), 7, double);

  for (int i = 0; i < nfrho - 1; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void init_style ()
{
  // convert read-in file(s) to arrays and spline them

  file2array ();
  array2spline ();

  cutmax = init_one (1, 1);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax,funcfl[m].cut);
  } 

  cutforcesq = cutmax*cutmax;

  nmax = nMol;

  return cutmax;
}

void printFuncfl (int i)
{
  int n;
  printf("%s %d %lg %d %lg %lg %lg\n", (funcfl[i].file), funcfl[i].nrho, funcfl[i].drho, funcfl[i].nr, funcfl[i].dr, funcfl[i].cut, funcfl[i].mass);
  for (n = 1; n < funcfl[i].nrho + 1; n ++){
    printf("%f ", funcfl[i].frho[n]);
  }
    printf("\n\n");
  for (n = 1; n < funcfl[i].nr + 1; n ++){
    printf("%f ", funcfl[i].zr[n]);
  }
    printf("\n\n");
  for (n = 1; n < funcfl[i].nr + 1; n ++){
    printf("%f ", funcfl[i].rhor[n]);
  }
    printf("\n\n");
}

void printfeam(int i)
{
  int n;

  if ((frhoeam_file = fopen ("out/md/eamfrho.dat", "a+")) == NULL){
    printf ("\nopen frho file error");
    getchar ();
    exit (1);
  }
  if ((zream_file = fopen ("out/md/eamzr.dat", "a+")) == NULL){
    printf ("\nopen zream file error");
    getchar ();
    exit (1);
  }
  if ((rhoream_file = fopen ("out/md/eamrhor.dat", "a+")) == NULL){
    printf ("\nopen rhoream file error");
    getchar ();
    exit (1);
  }

  fprintf(frhoeam_file, "rho f(eV)\n");
  for (n = 1; n < funcfl[i].nrho + 1;n ++){
    fprintf(frhoeam_file, "%f %f\n", funcfl[i].drho * (n - 1), funcfl[i].frho[n]);
  } 
  fprintf(zream_file, "r(A) z\n");
  for (n = 1; n < funcfl[i].nr + 1;n ++){
    fprintf(zream_file, "%f %f\n", funcfl[i].dr * (n - 1), funcfl[i].zr[n]);
  } 
  fprintf(rhoream_file, "r(A) rho\n");
  for (n = 1; n < funcfl[i].nr + 1;n ++){
    fprintf(rhoream_file, "%f %f\n", funcfl[i].dr * (n - 1), funcfl[i].rhor[n]);
  } 

  fclose (frhoeam_file);
  fclose (zream_file);
  fclose (rhoream_file);
}

void printfspline (double ***xxxspline)
{
  double *coeff;
  double p, rsqspline, rhospline, frhop, r, z2p, z2, recip, phi, phip;
  int m, n, i, itype, jtype, num;

  if ((spline_file = fopen ("out/md/eamspline.dat", "a+")) == NULL){
    printf ("\nopen spline file error");
    getchar ();
    exit (1);
  }

  num = 5000;
  AllocMem (spline, num, double);
  if (xxxspline == rhor_spline){
    for (n = 0; n < num; n++){
        if (n == 0) fprintf(spline_file, "r(A) rho\n");
        rsqspline = Sqr (n * funcfl[0].cut / num);
        itype = 1;
        jtype = 1;
        p = sqrt(rsqspline)*rdr + 1.0;
        m = (int) (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = xxxspline[type2rhor[jtype][itype]][m];
        spline[n] = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        fprintf (spline_file,"%e %e\n", sqrt(rsqspline), spline[n]);
    }
  }else if (xxxspline == frho_spline){
    for (i = 0; i < num; i ++){
      if (i == 0) fprintf(spline_file, "rho F(eV) dF_drho\n");
      rhospline = i * funcfl[0].nrho * funcfl[0].drho / num; 
      p = rhospline*rdrho + 1.0;
      m = (int) (p);
      m = MAX(1,MIN(m,nrho-1));
      p -= m;
      p = MIN(p,1.0);
      coeff = xxxspline[type2frho[1]][m];
      frhop = (coeff[0]*p + coeff[1])*p + coeff[2];
      spline[i] = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (rhospline > rhomax) spline[i] += frhop * (rhospline-rhomax);
      spline[i] *= scale[1][1];
      fprintf (spline_file,"%e %e %e\n", rhospline, spline[i], frhop);
    }
  }else if (xxxspline == z2r_spline){
    for (n = 0; n < num; n++){
        if (n == 0)  fprintf (spline_file, "r(A) z2\n");
        rsqspline = Sqr (n * funcfl[0].cut / num);
        itype = 1;
        jtype = 1;
        r = sqrt(rsqspline);
        p = r*rdr + 1.0;
        m = (int) (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = xxxspline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        fprintf (spline_file,"%e %e\n", sqrt(rsqspline), z2);
    }
  }else if (xxxspline == phi_spline){
    for (n = 1; n < num; n++){
        if (n == 1)  fprintf (spline_file, "r(A) phi(eV) dphi_dr\n");
        rsqspline = Sqr (n * funcfl[0].cut / num);
        itype = 1;
        jtype = 1;
        r = sqrt(rsqspline);
        p = r*rdr + 1.0;
        m = (int) (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        recip = 1.0/r;
        phi = z2*recip;
	phip = z2p*recip - phi*recip; //--eV/A
        fprintf (spline_file,"%e %e %e\n", sqrt(rsqspline), phi, phip);
    }
  }

  fclose (spline_file);
}

void ComputeForcesEamPoten ()
{
  VecR drVec, DR_A;
  double fcVal;
  double *coeff;
  double rsq,rr,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi,fpair;
  int m, n, itype, jtype, i, j, j1, j2;
  double mass;
  double *fx, *fy, *fz;

  nmax = nMol;
  AllocMem (rho, nmax, double);
  AllocMem (fp, nmax, double);

  // zero out density

  DO_MOL rho[n] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  #pragma omp parallel for private(j1, j2, i, j, drVec, DR_A, rsq, itype, jtype, p, m, coeff) num_threads(Nthreads)
  for (n = 0; n < nebrTabLen; n ++) {
    j1 = nebrTab[2 * n];
    j2 = nebrTab[2 * n + 1];
    i = j1;
    j = j2;
    VSub (drVec, mol[j1].r, mol[j2].r);
    if (strcmp(boundaryCond, "p") == 0)  VWrapAll (drVec); //boundary condition
    VSCopy(DR_A, (lUnit * 1.e10), drVec);
    rsq = VLenSq (DR_A);
//    printf("r = %f, cutforce = %f\n", sqrt(rsq), sqrt(cutforcesq));
    if (rsq < cutforcesq){
      itype = 1;
      jtype = 1;
      p = sqrt(rsq)*rdr + 1.0;
      m = (int) (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);
      coeff = rhor_spline[type2rhor[jtype][itype]][m];
      #pragma omp atomic
      rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      coeff = rhor_spline[type2rhor[itype][jtype]][m];
      #pragma omp atomic
      rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    }
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  frho_tot = 0.;
  #pragma omp parallel for reduction(+:frho_tot)private(p, m, coeff, phi) num_threads(Nthreads)
  for (i = 0; i < nMol; i ++){
    p = rho[i]*rdrho + 1.0;
    m = (int) (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[1]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
    phi *= scale[1][1];
    frho_tot += phi;
  }
//  printf("frho_tot_Avg = %feV\n", frho_tot / nMol);

  // compute forces on each atom
  // loop over neighbors of my atoms

//  DO_MOL VZero (mol[n].f);
  pairpoten_tot = 0.;
  uSum = 0.;
  virSum = 0.;
  AllocMem (fx, nMol, double);
  AllocMem (fy, nMol, double);
  AllocMem (fz, nMol, double);
  DO_MOL fx[n] = fy[n] = fz[n] = 0.;
  #pragma omp parallel for reduction(+:pairpoten_tot, virSum)private(j1, j2, i, j, drVec, DR_A, rsq, rr, itype, jtype, r, p, m, coeff, rhoip, rhojp, z2p, z2, recip, phi, phip, psip, fpair, fcVal) num_threads(Nthreads)
  for (n = 0; n < nebrTabLen; n ++) {
    j1 = nebrTab[2 * n];
    j2 = nebrTab[2 * n + 1];
    i = j1;
    j = j2;
    VSub (drVec, mol[j1].r, mol[j2].r); //dimensionless
    if (strcmp(boundaryCond, "p") == 0) VWrapAll (drVec);
    VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
    rsq = VLenSq (DR_A); //--A2
    if (rsq < cutforcesq){
        itype = 1;
        jtype = 1;
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = (int) (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r; //--A-1
        phi = z2*recip; //--eV
        phip = z2p*recip - phi*recip; //--eV/A
        psip = fp[i]*rhojp + fp[j]*rhoip + phip; //--eV/A
        fpair = -scale[itype][jtype]*psip*recip; //--eV/A2
        fcVal = fpair / ((eUnit / eleChar) / Sqr(lUnit * 1.e10)); //dimensionless
        #pragma omp atomic
	fx[j1] += fcVal * drVec.x;
        #pragma omp atomic
	fy[j1] += fcVal * drVec.y;
        #pragma omp atomic
	fz[j1] += fcVal * drVec.z;
        #pragma omp atomic
	fx[j2] += -fcVal * drVec.x;
        #pragma omp atomic
	fy[j2] += -fcVal * drVec.y;
        #pragma omp atomic
	fz[j2] += -fcVal * drVec.z;
        pairpoten_tot += phi; //--eV
	rr = rsq / Sqr(lUnit * 1.e10); //--dimensionless
	virSum += fabs(fcVal * rr); //--dimensionless
    }
  }
  
  DO_MOL VSet (mol[n].f, fx[n], fy[n], fz[n]); //--dimensionless

  uSum = frho_tot + pairpoten_tot; //--eV
//  printf("pairpoten_Avg = %f eV, uSum_Avg = %f eV, virSum_Avg = %f eV\n", pairpoten_tot / nMol, uSum / nMol, virSum * (eUnit / eleChar) / nMol);
  uSum = uSum / (eUnit / eleChar);//dimensionless

  free (rho);
  free (fp);
  free (fx);
  free (fy);
  free (fz);
} 
