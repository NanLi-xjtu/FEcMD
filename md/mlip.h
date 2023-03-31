#ifndef MLIP_H
#define MLIP_H

extern "C" {
#include "main.h"
}
#include "init.h"

extern void MLIP_init(const char*, const char*, int, double&, int&);
extern void MLIP_calc_cfg(int, double*, double**, int*, int*, double&, double**, double*);
extern void MLIP_calc_nbh(char*, double*, int, int*, int*, int**, int, int, double**, int*, double**, double&, double*, double**);
extern void MLIP_finalize();

class PairMLIP {
 public:
  double cutoff;

  PairMLIP();
  virtual ~PairMLIP();
  virtual void compute(int, int);
  void settings();
  virtual void coeff();
  void init_style();
  double init_one(int, int);

 protected:
  int mode; // 0 - nbh mode (can't learn on the fly), 1 - cfg mode (typically for non-parallel lammps)
  bool inited;
  char MLIPsettings_filename[1000];
  char MLIPlog_filename[1000];
  double cutoffsq;
  void allocate();
};

class AtomLMP {
 public:
  int nlocal, nghost;    // # of owned and ghost atoms on this proc
  int *type;
  int ntypes;
  double **x, **v, **f;
};

class NeighList {
 public:
  int inum;            // # of I atoms neighbors are stored for
  int *ilist;          // local indices of I atoms
  int *numneigh;       // # of J neighbors for each I atom
  int **firstneigh;    // ptr to 1st J int value of each I atom
};

extern PairMLIP mlp;

void ComputeForcesMTP ();

#endif
