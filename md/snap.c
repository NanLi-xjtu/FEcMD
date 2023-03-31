#include "snap.h"

  char snapFilename[128];
  int **sna_setflag, sna_nparas, sna_ncoeffq, sna_ncoeffall;;
  double **sna_cutsq;
  double sna_rcutmax;               // max cutoff for all sna_elements
  int sna_nelements;                // # of unique sna_elements
  char **sna_elements;              // names of unique sna_elements
  double *sna_radelem;              // element radii
  double *sna_wjelem;               // sna_elements weights
  double **sna_coeffelem;           // element bispectrum coefficients
  int *sna_map;                     // mapping from atom types to sna_elements
  int twojmax, diagonalstyle, switchflag, bzeroflag;
  double rfac0, rmin0, wj1, wj2;
  int rcutfacflag, twojmaxflag; // flags for required parameters
  double rcutfac;
  int quadraticflag; // declared public to workaround gcc 4.9
  int sna_ncoeff;                    //  compiler bug, manifest in KOKKOS package
  int sna_nmax = 0;

  double **sna_bvec, ***sna_dbvec;

  SNA* snaptr;
  double*** sna_uarraytot_r, *** sna_uarraytot_i;
  double***** sna_zarray_r, ***** sna_zarray_i;
  double*** sna_uarraytot_r_b, *** sna_uarraytot_i_b;
  double***** sna_zarray_r_b, ***** sna_zarray_i_b;
  double*** sna_uarray_r, *** sna_uarray_i;

  // Self-weight
  double sna_wself;
  // data for bispectrum coefficients

  double***** sna_cgarray;
  double** sna_rootpqarray;
  double*** sna_barray;

  // derivatives of data

  double**** sna_duarray_r, **** sna_duarray_i;
  double**** sna_dbarray;

  double *sna_bzero;  // array of B values for isolated atoms

  //use indexlist instead of loops, constructor generates these
  SNA_LOOPINDICES* sna_idxj;
  int sna_idxj_max;

  static const int sna_nmaxfactorial = 167;

  int flag_rij_memory = 0;
  int flag_sna_memory = 0;

  int sna_Nthreads = 7;
/* ----------------------------------------------------------------------
   factorial n table, size SNA::sna_nmaxfactorial+1
------------------------------------------------------------------------- */

const double nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // sna_nmaxfactorial = 167
};

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSNAP_allocate()
{
  int i, j, k;
  FILE *snap_file;
  char line[1024], filename[128];

  sprintf (filename, "%s.snapcoeff", snapFilename);
  snap_file = ReadFile(filename);
  for (i = 0; i < 4; i ++) fgets(line, 1024, snap_file);
  sscanf (line, "%d %d", &ntypes, &sna_nparas);
  puts (line);
  allocated = 1;
  int n = ntypes;
  sna_nelements = ntypes;

  AllocMem2 (sna_setflag, n+1, n+1, int);
  AllocMem2 (sna_cutsq, n+1, n+1, double);
  AllocMem (sna_map, n+1, int);

  fclose (snap_file);
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSNAP_coeff()
{
  int i, j, k;

  // read SNAP element names between 2 filenames
  // sna_nelements = # of SNAP sna_elements
  // sna_elements = list of unique element names

  if (!allocated) PairSNAP_allocate();

  if (sna_nelements < 1) {
    printf ("Incorrect sna_nelements for pair coefficients");
    exit (1);
  }

  AllocMem2 (sna_elements, sna_nelements, 10, char);

  char line[1024], filename[128];
  FILE *snap_file;

  sprintf (filename, "%s.snapcoeff", snapFilename);
  snap_file = ReadFile(filename);
  for (i = 0; i < 4; i ++) fgets(line, 1024, snap_file);
  for (i = 0; i < sna_nelements; i++) {
    fgets(line, 1024, snap_file);
    sscanf (line, "%s", sna_elements[i]);
    printf ("%s ", sna_elements[i]);
    for (j = 0; j < sna_nparas; j ++) fgets(line, 1024, snap_file);
  }
  printf ("\n");
  fclose (snap_file);

  // read snapcoeff and snapparam files
  char coefffilename[128], paramfilename[128];

  sprintf (coefffilename, "%s.snapcoeff", snapFilename);
  sprintf (paramfilename, "%s.snapparam", snapFilename);
  PairSNAP_read_files(coefffilename,paramfilename);

  if (!quadraticflag)
    sna_ncoeff = sna_ncoeffall - 1;
  else {

    // sna_ncoeffall should be (sna_ncoeff+2)*(sna_ncoeff+1)/2
    // so, sna_ncoeff = floor(sqrt(2*sna_ncoeffall))-1

    sna_ncoeff = sqrt(2*sna_ncoeffall)-1;
    sna_ncoeffq = (sna_ncoeff*(sna_ncoeff+1))/2;
    int ntmp = 1+sna_ncoeff+sna_ncoeffq;
    if (ntmp != sna_ncoeffall) {
      printf("sna_ncoeffall = %d ntmp = %d sna_ncoeff = %d \n",sna_ncoeffall,ntmp,sna_ncoeff);
      Error ("Incorrect SNAP coeff file");
    }
  }
  printf("sna_ncoeffall = %d sna_ncoeff = %d \n",sna_ncoeffall,sna_ncoeff);

  // read args that sna_map atom types to SNAP sna_elements
  // sna_map[i] = which element the Ith atom type is, -1 if not mapped
  // sna_map[0] is not used

  for (int i = 1; i <= ntypes; i++) {
    char* elemname = sna_elements[i-1];
    int jelem;
    for (jelem = 0; jelem < sna_nelements; jelem++)
      if (strcmp(elemname,sna_elements[jelem]) == 0)
        break;

    if (jelem < sna_nelements)
      sna_map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) sna_map[i] = -1;
  }

  int n = ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      sna_setflag[i][j] = 0;

  // set sna_setflag i,j for type pairs where both are mapped to sna_elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (sna_map[i] >= 0 && sna_map[j] >= 0) {
        sna_setflag[i][j] = 1;
        count++;
      }

  if (count == 0) Error("Incorrect args for pair coefficients");


  // Calculate maximum cutoff for all sna_elements

  sna_rcutmax = 0.0;
  for (int ielem = 0; ielem < sna_nelements; ielem++)
    sna_rcutmax = MAX(2.0*sna_radelem[ielem]*rcutfac,sna_rcutmax);

  printf ("\nrcutmax = %.3f\n", sna_rcutmax);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSNAP_init_style()
{
  int k;

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSNAP_init_one(int i, int j)
{
  if (sna_setflag[i][j] == 0) Error("All pair coeffs are not set");
  return (sna_radelem[sna_map[i]] +
          sna_radelem[sna_map[j]])*rcutfac;
}

/* ---------------------------------------------------------------------- */

void PairSNAP_read_files(char *coefffilename, char *paramfilename)
{
  int i, j, k;

  // open SNAP coefficient file on proc 0

  FILE *fpcoeff;

  fpcoeff = ReadFile (coefffilename);

  char line[1024],*ptr;
  int eof = 0;

  int n;

  // words = ptrs to all words in line
  // strip single and double quotes from words

  char* words[1024];
  int iword = 0;
  int nelemfile = sna_nelements;
  sna_ncoeffall = sna_nparas;

  // Set up element lists

  AllocMem (sna_radelem, sna_nelements, double);
  AllocMem (sna_wjelem, sna_nelements, double);
  AllocMem2 (sna_coeffelem, sna_nelements, sna_ncoeffall, double);

  int *found;
  AllocMem (found, sna_nelements, int);
  for (int ielem = 0; ielem < sna_nelements; ielem++)
    found[ielem] = 0;

  // Loop over sna_elements in the SNAP coefficient file
  for (i = 0; i < 4; i ++) fgets(line,1024,fpcoeff);
  for (int ielemfile = 0; ielemfile < nelemfile; ielemfile++) {

      ptr = fgets(line,1024,fpcoeff);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpcoeff);
      } else n = strlen(line) + 1;

    iword = 0;
    words[iword] = strtok(line,"' \t\n\r\f");
    iword = 1;
    words[iword] = strtok(NULL,"' \t\n\r\f");
    iword = 2;
    words[iword] = strtok(NULL,"' \t\n\r\f");

    char* elemtmp = words[0];
    double radtmp = atof(words[1]);
    double wjtmp = atof(words[2]);

    printf ("%s %.3f %.3f\n", elemtmp, radtmp, wjtmp);

    // skip if element name isn't in element list

    int ielem;
    for (ielem = 0; ielem < sna_nelements; ielem++)
      if (strcmp(elemtmp,sna_elements[ielem]) == 0) break;
    if (ielem == sna_nelements) {
        for (int icoeff = 0; icoeff < sna_ncoeffall; icoeff++)
          ptr = fgets(line,1024,fpcoeff);
      continue;
    }


    if (found[ielem]) {
        for (int icoeff = 0; icoeff < sna_ncoeffall; icoeff++)
          ptr = fgets(line,1024,fpcoeff);
      continue;
    }

    found[ielem] = 1;
    sna_radelem[ielem] = radtmp;
    sna_wjelem[ielem] = wjtmp;

    for (int icoeff = 0; icoeff < sna_ncoeffall; icoeff++) {
        ptr = fgets(line,1024,fpcoeff);
        if (ptr == NULL) {
          eof = 1;
          fclose(fpcoeff);
        } else n = strlen(line) + 1;
      iword = 0;
      words[iword] = strtok(line,"' \t\n\r\f");
      sna_coeffelem[ielem][icoeff] = atof(words[0]);

    }
  }

  // set flags for required keywords

  rcutfacflag = 0;
  twojmaxflag = 0;

  // Set defaults for optional keywords

  rfac0 = 0.99363;
  rmin0 = 0.0;
  diagonalstyle = 3;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  fpparam = ReadFile (paramfilename);
  eof = 0;
  for (i = 0; i < 3; i ++) fgets(line,1024,fpparam);
  while (1) {
      ptr = fgets(line,1024,fpparam);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpparam);
      } else n = strlen(line) + 1;
    if (eof) break;

    // words = ptrs to all words in line
    // strip single and double quotes from words

    char* keywd = strtok(line,"' \t\n\r\f");
    char* keyval = strtok(NULL,"' \t\n\r\f");
    if (strcmp(keywd,"rcutfac") == 0) {
      rcutfac = atof(keyval);
      rcutfacflag = 1;
    } else if (strcmp(keywd,"twojmax") == 0) {
      twojmax = atoi(keyval);
      twojmaxflag = 1;
    } else if (strcmp(keywd,"rfac0") == 0)
      rfac0 = atof(keyval);
    else if (strcmp(keywd,"rmin0") == 0)
      rmin0 = atof(keyval);
    else if (strcmp(keywd,"diagonalstyle") == 0)
      diagonalstyle = atoi(keyval);
    else if (strcmp(keywd,"switchflag") == 0)
      switchflag = atoi(keyval);
    else if (strcmp(keywd,"bzeroflag") == 0)
      bzeroflag = atoi(keyval);
    else if (strcmp(keywd,"quadraticflag") == 0)
      quadraticflag = atoi(keyval);
    else {
      printf ("Incorrect SNAP parameter file");
      exit (1);
    }
  }
  printf ("rcutfac = %.3f, twojmax = %d, rfac0 = %.3f, rmin0 = %.3f\n", rcutfac, twojmax, rfac0, rmin0);
  printf ("diagonalstyle = %d, switchflag = %d, bzeroflag = %d, quadraticflag = %d\n", \
           diagonalstyle, switchflag, bzeroflag, quadraticflag);

}

void SNA_free_rij()
{
  DestroyDouble2 (snaptr->rij);
  DestroyInt1 (snaptr->inside);
  DestroyDouble1 (snaptr->wj);
  DestroyDouble1 (snaptr->rcutij);
}

void SNA_grow_rij(int newnmax)
{
  int i, j, k;
  long n, nbytes;
  double *data;

  if(newnmax <= sna_nmax) return;

  sna_nmax = newnmax;
  if (flag_rij_memory == 1) SNA_free_rij();
  snaptr->rij = CreateDouble2 (sna_nmax, 3);
  snaptr->inside = CreateInt1 (sna_nmax);
  snaptr->wj = CreateDouble1 (sna_nmax);
  snaptr->rcutij = CreateDouble1 (sna_nmax);
  ZeroDouble2 (snaptr->rij, sna_nmax, 3);
  ZeroInt1 (snaptr->inside, sna_nmax);
  ZeroDouble1 (snaptr->wj, sna_nmax);
  ZeroDouble1 (snaptr->rcutij, sna_nmax);

  flag_rij_memory = 1;
}

void ComputeSNAAtom ()
{
  // construct sna_cutsq

  double cut;
  cutmax = 0.0;
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*sna_radelem[i-1]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    sna_cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (sna_radelem[i-1]+sna_radelem[j-1])*rcutfac;
      sna_cutsq[i][j] = sna_cutsq[j][i] = cut*cut;
    }
  }
  
  for(int i = 1; i <= ntypes; i++) {
    for(int j = 1; j <= ntypes; j++)
	printf ("cut(%d, %d) = %.3f A\n", i, j, sqrt (sna_cutsq[i][j]));
  }
}

void SNA_create_twojmax_arrays()
{
  int i, j, k, l, m, m1, m2, n, nbytes;
  int jdim = twojmax + 1;
  double *data, **level4, ***level3, ****level2, **cube, ***plane, **plane3;

  sna_uarraytot_r = CreateDouble3 (jdim, jdim, jdim);
  sna_uarraytot_i = CreateDouble3 (jdim, jdim, jdim);
  sna_uarray_r = CreateDouble3 (jdim, jdim, jdim);
  sna_uarray_i = CreateDouble3 (jdim, jdim, jdim);
  sna_barray = CreateDouble3 (jdim, jdim, jdim);
  sna_rootpqarray = CreateDouble2 (jdim+1, jdim+1);
  sna_zarray_r = CreateDouble5 (jdim, jdim, jdim, jdim, jdim);
  sna_zarray_i = CreateDouble5 (jdim, jdim, jdim, jdim, jdim);
  sna_cgarray = CreateDouble5 (jdim, jdim, jdim, jdim, jdim);
  if (bzeroflag)
    AllocMem (sna_bzero, jdim, double);
  else
    sna_bzero = NULL;
  sna_duarray_r = CreateDouble4 (jdim, jdim, jdim, 3);
  sna_duarray_i = CreateDouble4 (jdim, jdim, jdim, 3);
  sna_dbarray = CreateDouble4 (jdim, jdim, jdim, 3);
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

void SNA_init_sna_rootpqarray()
{
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      sna_rootpqarray[p][q] = sqrt((double)(p)/q);
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

double factorial(int n)
{
  if (n < 0 || n > sna_nmaxfactorial) {
    char str[128];
    sprintf(str, "Invalid argument to factorial %d", n);
//    error->all(FLERR, str);
    Error (str);
  }

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

double deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

void SNA_init_clebsch_gordan()
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= twojmax; j2++)
      for (int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
        for (int m1 = 0; m1 <= j1; m1 += 1) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2 += 1) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) continue;

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                   / 2, -(j - j1 - bb2) / 2));
                z <= MIN((j1 + j2 - j) / 2,
                         MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(factorial((j1 + aa2) / 2) *
                        factorial((j1 - aa2) / 2) *
                        factorial((j2 + bb2) / 2) *
                        factorial((j2 - bb2) / 2) *
                        factorial((j  + cc2) / 2) *
                        factorial((j  - cc2) / 2) *
                        (j + 1));

            sna_cgarray[j1][j2][j][m1][m2] = sum * dcg * sfaccg;
          }
        }
}

void SNA_build_indexlist()
{
  if(diagonalstyle == 0) {
    int sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
          sna_idxj_count++;

    // indexList can be changed here

    AllocMem (sna_idxj, sna_idxj_count, SNA_LOOPINDICES);
    sna_idxj_max = sna_idxj_count;

    sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
          sna_idxj[sna_idxj_count].j1 = j1;
          sna_idxj[sna_idxj_count].j2 = j2;
          sna_idxj[sna_idxj_count].j = j;
          sna_idxj_count++;
        }
  }

  if(diagonalstyle == 1) {
    int sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j = 0; j <= MIN(twojmax, 2 * j1); j += 2) {
        sna_idxj_count++;
      }

    // indexList can be changed here

    AllocMem (sna_idxj, sna_idxj_count, SNA_LOOPINDICES);
    sna_idxj_max = sna_idxj_count;

    sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j = 0; j <= MIN(twojmax, 2 * j1); j += 2) {
        sna_idxj[sna_idxj_count].j1 = j1;
        sna_idxj[sna_idxj_count].j2 = j1;
        sna_idxj[sna_idxj_count].j = j;
        sna_idxj_count++;
      }
  }

  if(diagonalstyle == 2) {
    int sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++) {
      sna_idxj_count++;
    }

    // indexList can be changed here

    AllocMem (sna_idxj, sna_idxj_count, SNA_LOOPINDICES);
    sna_idxj_max = sna_idxj_count;

    sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++) {
      sna_idxj[sna_idxj_count].j1 = j1;
      sna_idxj[sna_idxj_count].j2 = j1;
      sna_idxj[sna_idxj_count].j = j1;
      sna_idxj_count++;
    }
  }

  if(diagonalstyle == 3) {
    int sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) sna_idxj_count++;

    // indexList can be changed here

    AllocMem (sna_idxj, sna_idxj_count, SNA_LOOPINDICES);
    sna_idxj_max = sna_idxj_count;

    sna_idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) {
            sna_idxj[sna_idxj_count].j1 = j1;
            sna_idxj[sna_idxj_count].j2 = j2;
            sna_idxj[sna_idxj_count].j = j;
            sna_idxj_count++;
          }
  }

}

void SNA_free ()
{
SNA_free_rij();
free (snaptr->bvec);
free2 (snaptr->dbvec);
free (snaptr);
free3 (sna_uarraytot_r);
free3 (sna_uarraytot_i);
free3 (sna_uarray_r);
free3 (sna_uarray_i);
free2 (sna_rootpqarray);
free5 (sna_zarray_r);
free5 (sna_zarray_i);
free5 (sna_cgarray);
if (bzeroflag) {
    free (sna_bzero);
  } else 
    sna_bzero = NULL;
free4 (sna_duarray_r);
free4 (sna_duarray_i);
free4 (sna_dbarray);
free (sna_idxj);
}

int SNA_compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    if(diagonalstyle == 0) {
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          ncount++;
    } else if(diagonalstyle == 1) {
      int j2 = j1;

      for (int j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2)
        ncount++;
    } else if(diagonalstyle == 2) {
      ncount++;
    } else if(diagonalstyle == 3) {
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) ncount++;
    }

  return ncount;
}

void SNA_init ()
{
  int i,j,k;
  int jdim = twojmax + 1;

  sna_wself = 1.0;
  flag_rij_memory = 0;
  sna_ncoeff = SNA_compute_ncoeff();
  if (flag_sna_memory == 0) AllocMem (snaptr, 1, SNA);
  if (flag_sna_memory == 0) SNA_create_twojmax_arrays();
  ZeroDouble3 (sna_uarraytot_r, jdim, jdim, jdim);
  ZeroDouble3 (sna_uarraytot_i, jdim, jdim, jdim);
  ZeroDouble3 (sna_uarray_r, jdim, jdim, jdim);
  ZeroDouble3 (sna_uarray_i, jdim, jdim, jdim);
  ZeroDouble3 (sna_barray, jdim, jdim, jdim);
  ZeroDouble2 (sna_rootpqarray, jdim+1, jdim+1);
  ZeroDouble5 (sna_zarray_r, jdim, jdim, jdim, jdim, jdim);
  ZeroDouble5 (sna_zarray_i, jdim, jdim, jdim, jdim, jdim);
  ZeroDouble5 (sna_cgarray, jdim, jdim, jdim, jdim, jdim);
  if (bzeroflag) ZeroDouble1 (sna_bzero, jdim);
  ZeroDouble4 (sna_duarray_r, jdim, jdim, jdim, 3);
  ZeroDouble4 (sna_duarray_i, jdim, jdim, jdim, 3);
  ZeroDouble4 (sna_dbarray, jdim, jdim, jdim, 3);

  SNA_init_clebsch_gordan();
  SNA_init_sna_rootpqarray();
  if (flag_sna_memory == 0) SNA_build_indexlist();

  if (bzeroflag) {
    double www = sna_wself*sna_wself*sna_wself;
    for(int j = 0; j <= twojmax; j++)
      sna_bzero[j] = www*(j+1);
  }

  if (flag_sna_memory == 0) AllocMem (snaptr->bvec, sna_ncoeff, double);
  if (flag_sna_memory == 0) AllocMem2 (snaptr->dbvec, sna_ncoeff, 3, double);
  ZeroDouble1 (snaptr->bvec, sna_ncoeff);
  ZeroDouble2 (snaptr->dbvec, sna_ncoeff, 3);

  flag_sna_memory = 1;
}

void SNA_zero_uarraytot()
{
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        sna_uarraytot_r[j][ma][mb] = 0.0;
        sna_uarraytot_i[j][ma][mb] = 0.0;
      }
}

void SNA_addself_uarraytot(double sna_wself_in)
{
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++) {
      sna_uarraytot_r[j][ma][ma] = sna_wself_in;
      sna_uarraytot_i[j][ma][ma] = 0.0;
    }
}

void SNA_compute_uarray(double x, double y, double z,
                         double z0, double r)
{
  double r0inv;
  double a_r, b_r, a_i, b_i;
  double rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2

  sna_uarray_r[0][0][0] = 1.0;
  sna_uarray_i[0][0][0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {

    // fill in left side of matrix layer from previous layer

    for (int mb = 0; 2*mb <= j; mb++) {
      sna_uarray_r[j][0][mb] = 0.0;
      sna_uarray_i[j][0][mb] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = sna_rootpqarray[j - ma][j - mb];
        sna_uarray_r[j][ma][mb] +=
          rootpq *
          (a_r * sna_uarray_r[j - 1][ma][mb] +
           a_i * sna_uarray_i[j - 1][ma][mb]);
        sna_uarray_i[j][ma][mb] +=
          rootpq *
          (a_r * sna_uarray_i[j - 1][ma][mb] -
           a_i * sna_uarray_r[j - 1][ma][mb]);

        rootpq = sna_rootpqarray[ma + 1][j - mb];
        sna_uarray_r[j][ma + 1][mb] =
          -rootpq *
          (b_r * sna_uarray_r[j - 1][ma][mb] +
           b_i * sna_uarray_i[j - 1][ma][mb]);
        sna_uarray_i[j][ma + 1][mb] =
          -rootpq *
          (b_r * sna_uarray_i[j - 1][ma][mb] -
           b_i * sna_uarray_r[j - 1][ma][mb]);
      }
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    int mbpar = -1;
    for (int mb = 0; 2*mb <= j; mb++) {
      mbpar = -mbpar;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
        mapar = -mapar;
        if (mapar == 1) {
          sna_uarray_r[j][j-ma][j-mb] = sna_uarray_r[j][ma][mb];
          sna_uarray_i[j][j-ma][j-mb] = -sna_uarray_i[j][ma][mb];
        } else {
          sna_uarray_r[j][j-ma][j-mb] = -sna_uarray_r[j][ma][mb];
          sna_uarray_i[j][j-ma][j-mb] = sna_uarray_i[j][ma][mb];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double SNA_compute_sfac(double r, double rcut)
{
  if (switchflag == 0) return 1.0;
  if (switchflag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double SNA_compute_dsfac(double r, double rcut)
{
  if (switchflag == 0) return 0.0;
  if (switchflag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

void SNA_add_uarraytot(double r, double wj, double rcut)
{
  double sfac;

  sfac = SNA_compute_sfac(r, rcut);

  sfac *= wj;

  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        sna_uarraytot_r[j][ma][mb] +=
          sfac * sna_uarray_r[j][ma][mb];
        sna_uarraytot_i[j][ma][mb] +=
          sfac * sna_uarray_i[j][ma][mb];
      }
}

/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

void SNA_compute_ui(int jnum)
{
  double rsq, r, x, y, z, z0, theta0;


  SNA_zero_uarraytot();
  SNA_addself_uarraytot(sna_wself);
  for(int j = 0; j < jnum; j++) {
    x = snaptr->rij[j][0];
    y = snaptr->rij[j][1];
    z = snaptr->rij[j][2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    theta0 = (r - rmin0) * rfac0 * MY_PI / (snaptr->rcutij[j] - rmin0);
    z0 = r / tan(theta0);

    SNA_compute_uarray(x, y, z, z0, r);
    SNA_add_uarraytot(r, snaptr->wj[j], snaptr->rcutij[j]);
  }
}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

void SNA_compute_zi()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        for ma = 0,...,j
  //          for mb = 0,...,jmid
  //            z(j1,j2,j,ma,mb) = 0
  //            for ma1 = Max(0,ma+(j1-j2-j)/2),Min(j1,ma+(j1+j2-j)/2)
  //              sumb1 = 0
  //              ma2 = ma-ma1+(j1+j2-j)/2;
  //              for mb1 = Max(0,mb+(j1-j2-j)/2),Min(j1,mb+(j1+j2-j)/2)
  //                mb2 = mb-mb1+(j1+j2-j)/2;
  //                sumb1 += cg(j1,mb1,j2,mb2,j) *
  //                  u(j1,ma1,mb1) * u(j2,ma2,mb2)
  //              z(j1,j2,j,ma,mb) += sumb1*cg(j1,ma1,j2,ma2,j)
  // compute_dbidrj() requires full j1/j2/j chunk of z sna_elements
  // use zarray j1/j2 symmetry
  int j1, j2, j, ma, mb, ma1, ma2, mb1, mb2;
  double sumb1_r, sumb1_i;
  #pragma omp parallel for private (j1, j2, j, ma, mb, ma1, ma2, mb1, mb2, sumb1_r, sumb1_i) num_threads(sna_Nthreads)
  for(j1 = 0; j1 <= twojmax; j1++)
    for(j2 = 0; j2 <= j1; j2++) {
      for(j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for(mb = 0; 2*mb <= j; mb++)
          for(ma = 0; ma <= j; ma++) {
            sna_zarray_r[j1][j2][j][ma][mb] = 0.0;
            sna_zarray_i[j1][j2][j][ma][mb] = 0.0;

            for(ma1 = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma1 <= MIN(j1, (2 * ma - j + j2 + j1) / 2); ma1++) {
              sumb1_r = 0.0;
              sumb1_i = 0.0;

              ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;

              for(mb1 = MAX(0, (2 * mb - j - j2 + j1) / 2);
              mb1 <= MIN(j1, (2 * mb - j + j2 + j1) / 2); mb1++) {

                mb2 = (2 * mb - j - (2 * mb1 - j1) + j2) / 2;
                sumb1_r += sna_cgarray[j1][j2][j][mb1][mb2] *
                  (sna_uarraytot_r[j1][ma1][mb1] * sna_uarraytot_r[j2][ma2][mb2] -
                   sna_uarraytot_i[j1][ma1][mb1] * sna_uarraytot_i[j2][ma2][mb2]);
                sumb1_i += sna_cgarray[j1][j2][j][mb1][mb2] *
                  (sna_uarraytot_r[j1][ma1][mb1] * sna_uarraytot_i[j2][ma2][mb2] +
                   sna_uarraytot_i[j1][ma1][mb1] * sna_uarraytot_r[j2][ma2][mb2]);
              } // end loop over mb1

              sna_zarray_r[j1][j2][j][ma][mb] +=
                sumb1_r * sna_cgarray[j1][j2][j][ma1][ma2];
              sna_zarray_i[j1][j2][j][ma][mb] +=
                sumb1_i * sna_cgarray[j1][j2][j][ma1][ma2];
            } // end loop over ma1
          } // end loop over ma, mb
      } // end loop over j
    } // end loop over j1, j2
}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
------------------------------------------------------------------------- */

void SNA_compute_bi()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++) {
      for(int j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        sna_barray[j1][j2][j] = 0.0;

        for(int mb = 0; 2*mb < j; mb++)
          for(int ma = 0; ma <= j; ma++)
            sna_barray[j1][j2][j] +=
              sna_uarraytot_r[j][ma][mb] * sna_zarray_r[j1][j2][j][ma][mb] +
              sna_uarraytot_i[j][ma][mb] * sna_zarray_i[j1][j2][j][ma][mb];

        // For j even, special treatment for middle column

        if (j%2 == 0) {
          int mb = j/2;
          for(int ma = 0; ma < mb; ma++)
            sna_barray[j1][j2][j] +=
              sna_uarraytot_r[j][ma][mb] * sna_zarray_r[j1][j2][j][ma][mb] +
              sna_uarraytot_i[j][ma][mb] * sna_zarray_i[j1][j2][j][ma][mb];
          int ma = mb;
          sna_barray[j1][j2][j] +=
            (sna_uarraytot_r[j][ma][mb] * sna_zarray_r[j1][j2][j][ma][mb] +
             sna_uarraytot_i[j][ma][mb] * sna_zarray_i[j1][j2][j][ma][mb])*0.5;
        }

        sna_barray[j1][j2][j] *= 2.0;
        if (bzeroflag)
          sna_barray[j1][j2][j] -= sna_bzero[j];
      }
    }
}
/* ----------------------------------------------------------------------
   copy Bi array to a vector
------------------------------------------------------------------------- */

void SNA_copy_bi2bvec()
{
  int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++)
    if(diagonalstyle == 0) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2) {
          snaptr->bvec[ncount] = sna_barray[j1][j2][j];
          ncount++;
        }
    } else if(diagonalstyle == 1) {
      j2 = j1;
      for(j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        snaptr->bvec[ncount] = sna_barray[j1][j2][j];
        ncount++;
      }
    } else if(diagonalstyle == 2) {
      j = j2 = j1;
      snaptr->bvec[ncount] = sna_barray[j1][j2][j];
      ncount++;
    } else if(diagonalstyle == 3) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) {
            snaptr->bvec[ncount] = sna_barray[j1][j2][j];
            ncount++;
          }
    }
}

/* ----------------------------------------------------------------------
   compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

void SNA_compute_duarray(double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut)
{
  double r0inv;
  double a_r, a_i, b_r, b_i;
  double da_r[3], da_i[3], db_r[3], db_i[3];
  double dz0[3], dr0inv[3], dr0invdr;
  double rootpq;

  double rinv = 1.0 / r;
  double ux = x * rinv;
  double uy = y * rinv;
  double uz = z * rinv;

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

  dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  for (int k = 0; k < 3; k++) {
    da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }

  da_i[2] += -r0inv;

  for (int k = 0; k < 3; k++) {
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }

  db_i[0] += -r0inv;
  db_r[1] += r0inv;

  sna_uarray_r[0][0][0] = 1.0;
  sna_duarray_r[0][0][0][0] = 0.0;
  sna_duarray_r[0][0][0][1] = 0.0;
  sna_duarray_r[0][0][0][2] = 0.0;
  sna_uarray_i[0][0][0] = 0.0;
  sna_duarray_i[0][0][0][0] = 0.0;
  sna_duarray_i[0][0][0][1] = 0.0;
  sna_duarray_i[0][0][0][2] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    for (int mb = 0; 2*mb <= j; mb++) {
      sna_uarray_r[j][0][mb] = 0.0;
      sna_duarray_r[j][0][mb][0] = 0.0;
      sna_duarray_r[j][0][mb][1] = 0.0;
      sna_duarray_r[j][0][mb][2] = 0.0;
      sna_uarray_i[j][0][mb] = 0.0;
      sna_duarray_i[j][0][mb][0] = 0.0;
      sna_duarray_i[j][0][mb][1] = 0.0;
      sna_duarray_i[j][0][mb][2] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = sna_rootpqarray[j - ma][j - mb];
        sna_uarray_r[j][ma][mb] += rootpq *
                               (a_r *  sna_uarray_r[j - 1][ma][mb] +
                                a_i *  sna_uarray_i[j - 1][ma][mb]);
        sna_uarray_i[j][ma][mb] += rootpq *
                               (a_r *  sna_uarray_i[j - 1][ma][mb] -
                                a_i *  sna_uarray_r[j - 1][ma][mb]);

        for (int k = 0; k < 3; k++) {
          sna_duarray_r[j][ma][mb][k] +=
            rootpq * (da_r[k] * sna_uarray_r[j - 1][ma][mb] +
                      da_i[k] * sna_uarray_i[j - 1][ma][mb] +
                      a_r * sna_duarray_r[j - 1][ma][mb][k] +
                      a_i * sna_duarray_i[j - 1][ma][mb][k]);
          sna_duarray_i[j][ma][mb][k] +=
            rootpq * (da_r[k] * sna_uarray_i[j - 1][ma][mb] -
                      da_i[k] * sna_uarray_r[j - 1][ma][mb] +
                      a_r * sna_duarray_i[j - 1][ma][mb][k] -
                      a_i * sna_duarray_r[j - 1][ma][mb][k]);
        }

        rootpq = sna_rootpqarray[ma + 1][j - mb];
        sna_uarray_r[j][ma + 1][mb] =
          -rootpq * (b_r *  sna_uarray_r[j - 1][ma][mb] +
                     b_i *  sna_uarray_i[j - 1][ma][mb]);
        sna_uarray_i[j][ma + 1][mb] =
          -rootpq * (b_r *  sna_uarray_i[j - 1][ma][mb] -
                     b_i *  sna_uarray_r[j - 1][ma][mb]);

        for (int k = 0; k < 3; k++) {
          sna_duarray_r[j][ma + 1][mb][k] =
            -rootpq * (db_r[k] * sna_uarray_r[j - 1][ma][mb] +
                       db_i[k] * sna_uarray_i[j - 1][ma][mb] +
                       b_r * sna_duarray_r[j - 1][ma][mb][k] +
                       b_i * sna_duarray_i[j - 1][ma][mb][k]);
          sna_duarray_i[j][ma + 1][mb][k] =
            -rootpq * (db_r[k] * sna_uarray_i[j - 1][ma][mb] -
                       db_i[k] * sna_uarray_r[j - 1][ma][mb] +
                       b_r * sna_duarray_i[j - 1][ma][mb][k] -
                       b_i * sna_duarray_r[j - 1][ma][mb][k]);
        }
      }
    }

    int mbpar = -1;
    for (int mb = 0; 2*mb <= j; mb++) {
      mbpar = -mbpar;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
        mapar = -mapar;
        if (mapar == 1) {
          sna_uarray_r[j][j-ma][j-mb] = sna_uarray_r[j][ma][mb];
          sna_uarray_i[j][j-ma][j-mb] = -sna_uarray_i[j][ma][mb];
          for (int k = 0; k < 3; k++) {
            sna_duarray_r[j][j-ma][j-mb][k] = sna_duarray_r[j][ma][mb][k];
            sna_duarray_i[j][j-ma][j-mb][k] = -sna_duarray_i[j][ma][mb][k];
          }
        } else {
          sna_uarray_r[j][j-ma][j-mb] = -sna_uarray_r[j][ma][mb];
          sna_uarray_i[j][j-ma][j-mb] = sna_uarray_i[j][ma][mb];
          for (int k = 0; k < 3; k++) {
            sna_duarray_r[j][j-ma][j-mb][k] = -sna_duarray_r[j][ma][mb][k];
            sna_duarray_i[j][j-ma][j-mb][k] = sna_duarray_i[j][ma][mb][k];
          }
        }
      }
    }
  }

  double sfac = SNA_compute_sfac(r, rcut);
  double dsfac = SNA_compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;

  int j, ma, mb;
  #pragma omp parallel for private (j, ma, mb) num_threads(sna_Nthreads)
  for (j = 0; j <= twojmax; j++)
    for (ma = 0; ma <= j; ma++)
      for (mb = 0; mb <= j; mb++) {
        sna_duarray_r[j][ma][mb][0] = dsfac * sna_uarray_r[j][ma][mb] * ux +
                                  sfac * sna_duarray_r[j][ma][mb][0];
        sna_duarray_i[j][ma][mb][0] = dsfac * sna_uarray_i[j][ma][mb] * ux +
                                  sfac * sna_duarray_i[j][ma][mb][0];
        sna_duarray_r[j][ma][mb][1] = dsfac * sna_uarray_r[j][ma][mb] * uy +
                                  sfac * sna_duarray_r[j][ma][mb][1];
        sna_duarray_i[j][ma][mb][1] = dsfac * sna_uarray_i[j][ma][mb] * uy +
                                  sfac * sna_duarray_i[j][ma][mb][1];
        sna_duarray_r[j][ma][mb][2] = dsfac * sna_uarray_r[j][ma][mb] * uz +
                                  sfac * sna_duarray_r[j][ma][mb][2];
        sna_duarray_i[j][ma][mb][2] = dsfac * sna_uarray_i[j][ma][mb] * uz +
                                  sfac * sna_duarray_i[j][ma][mb][2];
      }
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

void SNA_compute_duidrj(double* rij, double wj, double rcut)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij[0];
  y = rij[1];
  z = rij[2];
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcut - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
}
/* ----------------------------------------------------------------------
   calculate derivative of Bi w.r.t. atom j
   variant using indexlist for j1,j2,j
   variant using symmetry relation
------------------------------------------------------------------------- */

void SNA_compute_dbidrj()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        zdb = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            zdb +=
  //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
  //        dbdr(j1,j2,j) += 2*zdb
  //        zdb = 0
  //        for mb1 = 0,...,j1mid
  //          for ma1 = 0,...,j1
  //            zdb +=
  //              Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j1+1)
  //        zdb = 0
  //        for mb2 = 0,...,j2mid
  //          for ma2 = 0,...,j2
  //            zdb +=
  //              Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j2+1)

  double* dbdr;
  double* dudr_r, *dudr_i;
  double sumzdu_r[3];
  double** jjjzarray_r;
  double** jjjzarray_i;
  double jjjmambzarray_r;
  double jjjmambzarray_i;
  int JJ, j1, j2, j, k, ma, mb, ma1, mb1, ma2, mb2;
  double j1fac, j2fac;
  #pragma omp parallel for private (JJ, j1, j2, j, dbdr, sumzdu_r, jjjzarray_r, jjjzarray_i, ma, mb, dudr_r, dudr_i, jjjmambzarray_r, jjjmambzarray_i, k, j1fac, ma1, mb1, j2fac, ma2, mb2) num_threads(sna_Nthreads)
  for(JJ = 0; JJ < sna_idxj_max; JJ++) {
    j1 = sna_idxj[JJ].j1;
    j2 = sna_idxj[JJ].j2;
    j = sna_idxj[JJ].j;

    dbdr = sna_dbarray[j1][j2][j];
    dbdr[0] = 0.0;
    dbdr[1] = 0.0;
    dbdr[2] = 0.0;

    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    for(k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j1 >= j2) {
      jjjzarray_r = sna_zarray_r[j1][j2][j];
      jjjzarray_i = sna_zarray_i[j1][j2][j];
    } else {
      jjjzarray_r = sna_zarray_r[j2][j1][j];
      jjjzarray_i = sna_zarray_i[j2][j1][j];
    }

    for(mb = 0; 2*mb < j; mb++)
      for(ma = 0; ma <= j; ma++) {

        dudr_r = sna_duarray_r[j][ma][mb];
        dudr_i = sna_duarray_i[j][ma][mb];
        jjjmambzarray_r = jjjzarray_r[ma][mb];
        jjjmambzarray_i = jjjzarray_i[ma][mb];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {
      mb = j/2;
      for(ma = 0; ma < mb; ma++) {
        dudr_r = sna_duarray_r[j][ma][mb];
        dudr_i = sna_duarray_i[j][ma][mb];
        jjjmambzarray_r = jjjzarray_r[ma][mb];
        jjjmambzarray_i = jjjzarray_i[ma][mb];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;
      }
      ma = mb;
      dudr_r = sna_duarray_r[j][ma][mb];
      dudr_i = sna_duarray_i[j][ma][mb];
      jjjmambzarray_r = jjjzarray_r[ma][mb];
      jjjmambzarray_i = jjjzarray_i[ma][mb];
      for(k = 0; k < 3; k++)
        sumzdu_r[k] +=
          (dudr_r[k] * jjjmambzarray_r +
           dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if jeven

    for(k = 0; k < 3; k++)
      dbdr[k] += 2.0*sumzdu_r[k];

    // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

    j1fac = (j+1)/(j1+1.0);

    for(k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j >= j2) {
      jjjzarray_r = sna_zarray_r[j][j2][j1];
      jjjzarray_i = sna_zarray_i[j][j2][j1];
    } else {
      jjjzarray_r = sna_zarray_r[j2][j][j1];
      jjjzarray_i = sna_zarray_i[j2][j][j1];
    }

    for(mb1 = 0; 2*mb1 < j1; mb1++)
      for(ma1 = 0; ma1 <= j1; ma1++) {

        dudr_r = sna_duarray_r[j1][ma1][mb1];
        dudr_i = sna_duarray_i[j1][ma1][mb1];
        jjjmambzarray_r = jjjzarray_r[ma1][mb1];
        jjjmambzarray_i = jjjzarray_i[ma1][mb1];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma1 mb1

    // For j1 even, handle middle column

    if (j1%2 == 0) {
      mb1 = j1/2;
      for(ma1 = 0; ma1 < mb1; ma1++) {
        dudr_r = sna_duarray_r[j1][ma1][mb1];
        dudr_i = sna_duarray_i[j1][ma1][mb1];
        jjjmambzarray_r = jjjzarray_r[ma1][mb1];
        jjjmambzarray_i = jjjzarray_i[ma1][mb1];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;
      }
      ma1 = mb1;
      dudr_r = sna_duarray_r[j1][ma1][mb1];
      dudr_i = sna_duarray_i[j1][ma1][mb1];
      jjjmambzarray_r = jjjzarray_r[ma1][mb1];
      jjjmambzarray_i = jjjzarray_i[ma1][mb1];
      for(k = 0; k < 3; k++)
        sumzdu_r[k] +=
          (dudr_r[k] * jjjmambzarray_r +
           dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if j1even

    for(k = 0; k < 3; k++)
      dbdr[k] += 2.0*sumzdu_r[k]*j1fac;

    // Sum over Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)

    j2fac = (j+1)/(j2+1.0);

    for(k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j1 >= j) {
      jjjzarray_r = sna_zarray_r[j1][j][j2];
      jjjzarray_i = sna_zarray_i[j1][j][j2];
    } else {
      jjjzarray_r = sna_zarray_r[j][j1][j2];
      jjjzarray_i = sna_zarray_i[j][j1][j2];
    }

    for(mb2 = 0; 2*mb2 < j2; mb2++)
      for(ma2 = 0; ma2 <= j2; ma2++) {

        dudr_r = sna_duarray_r[j2][ma2][mb2];
        dudr_i = sna_duarray_i[j2][ma2][mb2];
        jjjmambzarray_r = jjjzarray_r[ma2][mb2];
        jjjmambzarray_i = jjjzarray_i[ma2][mb2];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma2 mb2

    // For j2 even, handle middle column

    if (j2%2 == 0) {
      mb2 = j2/2;
      for(ma2 = 0; ma2 < mb2; ma2++) {
        dudr_r = sna_duarray_r[j2][ma2][mb2];
        dudr_i = sna_duarray_i[j2][ma2][mb2];
        jjjmambzarray_r = jjjzarray_r[ma2][mb2];
        jjjmambzarray_i = jjjzarray_i[ma2][mb2];
        for(k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
            dudr_i[k] * jjjmambzarray_i;
      }
      ma2 = mb2;
      dudr_r = sna_duarray_r[j2][ma2][mb2];
      dudr_i = sna_duarray_i[j2][ma2][mb2];
      jjjmambzarray_r = jjjzarray_r[ma2][mb2];
      jjjmambzarray_i = jjjzarray_i[ma2][mb2];
      for(k = 0; k < 3; k++)
        sumzdu_r[k] +=
          (dudr_r[k] * jjjmambzarray_r +
           dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if j2even

    for(k = 0; k < 3; k++) {
      dbdr[k] += 2.0*sumzdu_r[k]*j2fac;
    }

  } //end loop over j1 j2 j

}

/* ----------------------------------------------------------------------
   copy Bi derivatives into a vector
------------------------------------------------------------------------- */

void SNA_copy_dbi2dbvec()
{
  int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++) {
    if(diagonalstyle == 0) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2) {
          snaptr->dbvec[ncount][0] = sna_dbarray[j1][j2][j][0];
          snaptr->dbvec[ncount][1] = sna_dbarray[j1][j2][j][1];
          snaptr->dbvec[ncount][2] = sna_dbarray[j1][j2][j][2];
          ncount++;
        }
    } else if(diagonalstyle == 1) {
      j2 = j1;
      for(j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        snaptr->dbvec[ncount][0] = sna_dbarray[j1][j2][j][0];
        snaptr->dbvec[ncount][1] = sna_dbarray[j1][j2][j][1];
        snaptr->dbvec[ncount][2] = sna_dbarray[j1][j2][j][2];
        ncount++;
      }
    } else if(diagonalstyle == 2) {
      j = j2 = j1;
      snaptr->dbvec[ncount][0] = sna_dbarray[j1][j2][j][0];
      snaptr->dbvec[ncount][1] = sna_dbarray[j1][j2][j][1];
      snaptr->dbvec[ncount][2] = sna_dbarray[j1][j2][j][2];
      ncount++;
    } else if(diagonalstyle == 3) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) {
            snaptr->dbvec[ncount][0] = sna_dbarray[j1][j2][j][0];
            snaptr->dbvec[ncount][1] = sna_dbarray[j1][j2][j][1];
            snaptr->dbvec[ncount][2] = sna_dbarray[j1][j2][j][2];
            ncount++;
          }
    }
  }
}

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

void PairSNAP_compute_regular(int eflag, int vflag)
{
  int i,j,k,n,jnum,ninside;
  double delx,dely,delz,evdwl,rsq;
  double fij[3], fijx, fijy, fijz;
  int *jlist,*numneigh,**firstneigh;
  evdwl = 0.0;
  VecR drVec, DR_A, f;
  double start, end, time1, time2, time3;
  double bgb;

  SNA_init ();
  DO_MOL VZero(mol[n].f);
  time1 = time2 = time3 = 0.;
  for (int ii = 0; ii < nMol; ii ++) {
    i = ii;
    const double xtmp = mol[i].r.x * lUnit *1.e10; //--A
    const double ytmp = mol[i].r.y * lUnit *1.e10; //--A
    const double ztmp = mol[i].r.z * lUnit *1.e10; //--A

    int itype = 0;
    for (k = 1; k <= ntypes; k ++) {
      if (strcmp(mol[i].elem, sna_elements[k-1]) == 0) itype = k;
    }
    if (itype == 0) Error ("error(snap.c):itype, the snap potential file does not match\n");

    const int ielem = sna_map[itype];
    const double radi = sna_radelem[ielem];

    jlist = mol[i].id_nebr;
    jnum = mol[i].Nnebr;

    // insure rij, inside, wj, and rcutij are of size jnum

    SNA_grow_rij(jnum);
    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    start = omp_get_wtime ();
    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      VSub (drVec, mol[i].r, mol[j].r); //--dimensionless
      if (strcmp(boundaryCond, "p") == 0) VWrapAll (drVec); //--dimensionless
      VSCopy (DR_A, (lUnit * 1.e10), drVec); //--A
      delx = DR_A.x; //--A
      dely = DR_A.y;
      delz = DR_A.z;

      rsq = delx*delx + dely*dely + delz*delz;

      int jtype = 0;
      for (k = 1; k <= ntypes; k ++) {
        if (strcmp(mol[j].elem, sna_elements[k-1]) == 0) jtype = k;
      }
      if (jtype == 0) Error ("error(snap.c):jtype, the snap potential file does not match\n");

      int jelem = sna_map[jtype];

      if (rsq < sna_cutsq[itype][jtype]&&rsq>1e-20) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
        snaptr->wj[ninside] = sna_wjelem[jelem];
        snaptr->rcutij[ninside] = (radi + sna_radelem[jelem])*rcutfac;
        ninside++;
      }
    }
    time1 += (omp_get_wtime () - start);

    // compute Ui, Zi, and Bi for atom I
    start = omp_get_wtime ();
    SNA_compute_ui(ninside);
    SNA_compute_zi();
    if (quadraticflag) {
      SNA_compute_bi();
      SNA_copy_bi2bvec();
    }
    time2 += omp_get_wtime () - start;

    // for neighbors of I within cutoff:
    // compute dUi/drj and dBi/drj
    // Fij = dEi/dRj = -dEi/dRi => add to Fi, subtract from Fj

    double* coeffi = sna_coeffelem[ielem];
    start = omp_get_wtime ();
    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      SNA_compute_duidrj(snaptr->rij[jj],
                             snaptr->wj[jj],snaptr->rcutij[jj]);

      SNA_compute_dbidrj();
      SNA_copy_dbi2dbvec();

      fijx = 0.0;
      fijy = 0.0;
      fijz = 0.0;

      // linear contributions
      for (k = 1; k <= sna_ncoeff; k++) {
        bgb = coeffi[k];
        fijx += bgb*snaptr->dbvec[k-1][0];
        fijy += bgb*snaptr->dbvec[k-1][1];
        fijz += bgb*snaptr->dbvec[k-1][2];
      }

      // quadratic contributions

      if (quadraticflag) {
        int k = sna_ncoeff+1;
        for (int icoeff = 0; icoeff < sna_ncoeff; icoeff++) {
          double bveci = snaptr->bvec[icoeff];
          double fack = coeffi[k]*bveci;
          double* dbveci = snaptr->dbvec[icoeff];
          fijx += fack*dbveci[0];
          fijy += fack*dbveci[1];
          fijz += fack*dbveci[2];
          k++;
          for (int jcoeff = icoeff+1; jcoeff < sna_ncoeff; jcoeff++) {
            double facki = coeffi[k]*bveci;
            double fackj = coeffi[k]*snaptr->bvec[jcoeff];
            double* dbvecj = snaptr->dbvec[jcoeff];

            fijx += facki*dbvecj[0]+fackj*dbveci[0];
            fijy += facki*dbvecj[1]+fackj*dbveci[1];
            fijz += facki*dbvecj[2]+fackj*dbveci[2];
            k++;
          }
        }
      }
      VSet (f, fijx, fijy, fijz);
      VVSAdd (mol[i].f, lUnit * 1.e10 / (eUnit / eleChar), f); //--dimensionless
      VVSAdd (mol[j].f, - lUnit * 1.e10 / (eUnit / eleChar), f); //--dimensionless
/*
      // tally per-atom virial contribution

      if (vflag)
        ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                     fij[0],fij[1],fij[2],
                     -snaptr->rij[jj][0],-snaptr->rij[jj][1],
                     -snaptr->rij[jj][2]);
*/
    }
    time3 += omp_get_wtime () - start;
/*
    // tally energy contribution

    if (eflag) {

      // evdwl = energy of atom I, sum over coeffs_k * Bi_k

      evdwl = coeffi[0];
      if (!quadraticflag) {
        snaptr->compute_bi();
        snaptr->copy_bi2bvec();
      }

      // E = beta.B + 0.5*B^t.alpha.B
      // coeff[k] = beta[k-1] or
      // coeff[k] = alpha_ii or
      // coeff[k] = alpha_ij = alpha_ji, j != i

      // linear contributions

      for (int k = 1; k <= sna_ncoeff; k++)
        evdwl += coeffi[k]*snaptr->bvec[k-1];

      // quadratic contributions

      if (quadraticflag) {
        int k = sna_ncoeff+1;
        for (int icoeff = 0; icoeff < sna_ncoeff; icoeff++) {
          double bveci = snaptr->bvec[icoeff];
          evdwl += 0.5*coeffi[k++]*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < sna_ncoeff; jcoeff++) {
            evdwl += coeffi[k++]*bveci*snaptr->bvec[jcoeff];
          }
        }
      }
      ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
    }
*/
  }
  printf ("time: %f\n", time1+time2+time3);
}

