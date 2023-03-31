1. Introduction

Field emission coupled with molecular dynamics simulation (FEcMD) software package is a computational tool for studying the electron emission characteristics and the atomic structure evolution of micro- and nano-protrusions made of pure metals or multi-component alloys by means of multi-physics and multi-scale methodology. 

2. FEcMD Preparation

2.1 Software environment

Only Ubuntu 16 is supported. The g++, gfortran and cmake compilers need to be installed in advance. 

    $ sudo apt-get update
    $ sudo apt-get install g++
    $ sudo apt-get install gfortran
    $ sudo apt-get install cmake
	
vim and unzip sould be also installed with:

    $ sudo apt-get install vim
    $ sudo apt-get install unzip
	
The software OVITO can be used to visualize the results files "*.movie" or "*.xyz". The software PARAVIEW can be used to visualize the results files "*.vtk".

    $ sudo apt-get install ovito
    $ sudo apt-get install paraview

2.2 Installation
    Download the file FEcMD.zip.
    There are two approaches to install FEcMD (select one of them):

(1) Extract the executable without modifying the source code.

## Run the following commands to extract and access to all of the original files:
    $ unzip FEcMD.zip 
    $ sudo chmod -R 775 FEcMD
    $ cd FEcMD

## The following are documents and folders in FEcMD after extracting:
--fecmd: executable file for this program

--example: this directory contains the example in the paper, user can run the *.sh to test examples.

--in: this directory contains the input files

--obj: this directory contains the object files during compilation

--out: this directory contains the output files

--lib: this directory contains the libraries.

## Then configure, compile, and install the deal.II library with:
    $ cd lib/deal.II/dealii-9.2.0
    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_INSTALL_PREFIX=../../ ../
    $ make install    (alternatively $ make -j<N> install)
    $ make test
Add the dynamic library directory to the configuration file of the shared library:

    $ sudo vim /etc/ld.so.conf
	
Add ‘/PATH/TO/FEcMD/lib/deal.II/lib’ below the ‘include ld.so.conf.d/*.conf’, save and quit.

    $ sudo /sbin/ldconfig -v
    $ sudo ldconfig
	
## Next, extract and run the executable
    $ cd ../../../../
	$ unzip executable_file.zip
	$ cp executable_file/fecmd ./
	$ sudo chmod 775 fecmd
    $ ./fecmd.
## How to run examples:
    $ cd example/the/folder/contains/*.sh
    $ bash *.sh

(2) Install the program by compiling the FEcMD source code.

This method is relative complicated, but if you need to modify the source code, use this method to compile and install FEcMD.

Compile deal.II, libxc, FEMOCS, and MLIP library, respectively. Then, place the static libraries in the lib directory and link the dynamic libraries. The detailed steps are in the README or Makefile in the corresponding directories.

## Compile and install deal.II library. 

The method refer to (1).

## Compile and install libxc library with:
    $ cd lib/libxc
    $ cmake -H. -Bobjdir
    $ cd objdir && make
    $ make test
    $ sudo make install 
## Compile and install FEMOCS library and GETELE library with:
    $ cd ../../femocs/GETELEC
    $ make clean
    $ make
    $ cd ../
    $ make clean
    $ make
## Compile and install MLIP library with:
    $ cd ../mlip-2-master
    $ ./configure
    $ make mlp
    $ make libinterface
## Install the following essential libraries:
    $ sudo apt-get install libtbb-dev
    $ sudo apt-get install libboost-dev
    $ sudo apt-get install liblapack-dev
    $ sudo apt-get install libz-dev
## After that, link all library functions to compile and install FEcMD 
    $ cd ../../
    $ make

3. Run FEcMD

3.1 Commands

There are various commands to perform the corresponding function.

Perform ED-MD simulation:

    $ ./fecmd
	
Perform ED simulation:

    $ ./fecmd ED
	
Perform MD simulation:

    $ ./fecmd MD
	
Perform Maxwell stress MD simulation:

    $ ./fecmd Maxwell
	
Perform heat transport simulation:

    $ ./fecmd heattransport
	
Generate Pyramidal-hemisphere tip:

    $ ./fecmd PH
	
Generate Mushroom-head tip:

    $ ./fecmd MH
	
Generate Prolate-spheroidal tip:

    $ ./fecmd PS
	
Generate Hemi-ellipsoidal tip:

    $ ./fecmd HE
	
Generate Cylinder tip:

    $ ./fecmd Cylinder
	
Generate Cone tip:

    $ ./fecmd Cone

3.2 Examples

Five examples in the corresponding paper were attached in FEcMD/example. Change to the example directory and run the bash script:

    $ cd the/folder/contains/*.sh/in/example
    $ bash *.sh

3.3 input and output files

The following is the input file and output file documents.

## Input file:
Definition of variable

===============Example for md.in:==============
#==========================================
# MD parameters
# Commands and arguments should be separated with " = " sign
#------------------------------------------
deltaT = 4           # The time interval of each step [fs]
temperature = 300      # Initial temperature [K]
Trate = 0              # Rate of temperature rise per step [K]
vMax_set = 1000000000  # Maximim atomic velocity [K]
vMin_set = 0           # Minimim atomic velocity [K]
Nelems = 1             # The number of elements
elems = Cu             # Elements types
mass = 63.546          # Relative atom mass for each elements
stepAvg = 10           # The average steps of statistical energy and temperature
stepMovie = 100        # The average steps of the output atomic trajectory
stepEquil = 1000       # The steps for inital equilibrium
stepLimit = 100000     # Total number of simulation steps
randSeed = 7           # Random seed
boundary = p           # Boundary condition: periodic (p) or nonperiodic (n)
boundary_fly = 0       # Whether to delete the atoms flying from the surface: 
                       # yes (1) or no (0)
regionScale = 1 1 1    # MD box scaling factor
structure_type = FCC   # Structure type: FCC, BCC or cubic 
                       # (consistent with the input atomic structure file)
lattice = 3.6147       # Lattice parameters [Angstrom]
interact_method = nebr # Statistical method of interacting atom pairs: 
                       # allpairs, cell or nebr
force_type = MTP       # Types of atomic interactions: 
                       # lj, metal, alloy, snap, eamfs, MTP
force_file = in/potential/Cu.mtp  # Atomic interaction potential file name
lj_strength = 0.6      # L-J strength length parameter [eV]
lj_length = 3.4        # L-J potential length parameter [Angstrom]
ensemble = NVT         # Ensemble: NVE, NVT, NPT or nonEquil
Nthreads = 8           # Parallel threads for MD
pedestal_thick = 0     # Fixed atomic thickness at the bottom
#tip Maker
#------------------------------------------
cell_order = disorder  # Atomic ordering of the alloys constructed:
                       # order or disorder
cell_file = in/mdlat.in.xyz  # Field emitter atomic strcture file name
initUcell = 10 10 10   # The size of cell expansion
tip_file = in/mdlat.in.xyz   # The generated field emitter atomic strcture file name
tip_R = 30             # The radius of the mushroom cap R [Angstrom]
tip_r = 15            # The radius of cone top rt [Angstrom]
tip_h = 500            # Total hight [Angstrom]
tip_theta = 3          # Half-angle [degree]
#neighbor list
#------------------------------------------
nebrTabFac = 200       # The Maximum number of adjacent atoms in the cutoff radius
rNebrShell = 0.3       # The radius of shell for neighbor list method [Angstrom]
#adjust temperature
#------------------------------------------
stepAdjustTemp = 10    # Average steps for adjusting the temperature
stepInitlzTemp = 1     # Average steps for adjusting the inital temperature
#ligancy analysis && mark surface atoms
#------------------------------------------
rSurfCut = 2.951       # The cut-off radius for extracting surface atoms [Angstrom]
nSurfCut = 8           # The atomic number threshold for extracting surface atoms
rFlyCut = 3.818        # The cut-off radius for extracting atoms flying out of the surface [Angstrom]
nFlyCut = 2            # The atomic number threshold for extracting atoms flying out of the surface
stepSurf = 8           # Average steps for extracting surface atoms
stepLigancy = 8        # Average steps for calculating the coordination number
#velocity distribution
#------------------------------------------
rangeVel = 10.         # The range of velocity distribution [Angstrom/ps]
limitVel = 100         # The limit steps of velocity distribution
stepVel = 5            # The average steps of velocity distribution
sizeHistVel = 200      # The size of histogram
#RDF
#------------------------------------------
limitRdf = 100         # The limit steps of RDF
stepRdf = 10           # The average steps of RDF
rangeRdf = 4.          # The range of RDF [Angstrom]
sizeHistRdf = 200      # The size of histogram
#NPT
#------------------------------------------
extPressure = 26.2     # External pressure [GPa]
massS = 100            # mass parameter whose value must be determined empirically
massV = 0.000001       # mass parameter whose value must be determined empirically
#Maxwell stress 
#------------------------------------------
Maxwell_rate = 0.      # The rise rate of Maxwell stress [GPa/ps]
Maxwell_max = 10       # The Maximum of Maxwell stress [GPa]
Maxwell_begin = 0      # The initial Maxwell stress [GPa]
#lattice conductivity 
#------------------------------------------
gravField = 0.
wallTempHi = 800.      # The upper temperature limit of heat conduction
wallTempLo = 300.      # The lower temperature limit of heat conduction
#T/V-profile
#------------------------------------------
limitGrid = 100        # The limit steps of T/V-profile
stepGrid = 50          # The average steps of T/V-profile
sizeHistGrid = 1 1 25  # The size of histogram grids


===============Example for femocs.in: ===============
#==========================================
# FEMOCS parameters
#------------------------------------------
infile = out/md/femocs.in.xyz
# General
#---------------
md_timestep = 4               # MD time step [fs]
time_limit = 200                 # simulation time limit [ps]
movie_timestep = 1000		# time between frames of movie files; 0 turns writing off [fs]
latconst = 3.615                # lattice constant [A]
radius = 70                     # inner radius of coarsening cylinder [A]
distance_tol = 0.45             # max RMS distance atoms are allowed to move between runs before the solution is recalculated; 0 forces to recalculate every time-step
n_read_conf = 0                 # nr of timesteps between re-reading conf file; 0 reads file only once
femocs_verbose_mode = verbose    # control verbosity in console; mute, silent or verbose
clear_output = true             # clear out folder before the run
coarse_theta = 10               # apex angle of coarsening cylinder [deg]
coarse_factor = 0.5 12 2         # coarsening factor; bigger number gives coarser surface
                                # 1st - coarsening factor for atoms outside the coarsening cylinder
                                # 2nd - minimum distance between atoms in nanotip below apex [latconst/4]
                                # 3rd - minimum distance between atoms in nanotip apex [latconst/4]

# Atom processing
#---------------
nnn = 12                       # nr of nearest neighbours of bulk material within coord_cutoff radius; needs to be adjusted if coord_cutoff is changed
coord_cutoff = 3.1              # coordination analysis cut-off radius [A]
cluster_cutoff = 4.2             # cluster anal. cut-off radius [A]; if 0, cluster anal. uses coord_cutoff instead
mesh_quality = 1.8              # minimum mesh quality Tetgen is allowed to make; smaller nr gives more symmetric elements
coplanarity = 1.0               # parameter defining the max flatness of tetrahedron
surface_smooth_factor = 1.1     # surface smoothing factor; bigger number gives smoother but more distorted surface
smooth_steps = 3                # number of surface mesh smoothing iterations; 0 turns smoothing off
cluster_anal = ture            # enable or disable cluster analysis
clean_surface = ture           # activate extra effort to clean surface atoms
charge_smooth_factor = 100.0    # parameter controlling face-charge distribution

# Parameters for extending simulation domain
#---------------
extended_atoms = in/extension.xyz   # file with atoms of extended surface
femocs_periodic = false          # imported atoms have periodic boundaries in x- & y-direction; must be false in systems without slab
box_width = 10                  # minimal simulation box width [tip height]
box_height = 10                  # simulation box height [tip height]
bulk_height = 20                # bulk substrate height [latconst]

# Field emission parameters
#---------------
work_function = 4.59             # work function [eV]
emitter_blunt = false            # if true blunt emitter SN barrier approximation used
emitter_cold = false
emitter_qe = true        #if true quantum effect of space charge will be considered
func_X_id = 1            #exchange functional id
func_C_id = 9            #correlation functional id

# Heating parameters
#---------------
heat_mode = transient           # method to calculate heating effects in material; none, stationary, transient, converge
force_mode = all                # forces to be calculated; all, lorentz, none
vscale_tau = 1500.0             # time constant in temperature scaler
#t_ambient =  300.0             # temperature at the bottom of simulation cell [K]
heat_limit = 0 1.e5             # Minimum amd Maxmum allowed temperature [K]
lorentz = 2.e-8                 # Lorentz constant [W Ohm K-2]
heat_cp = 3.491e-24             # Volumetric heat capacity [J/(K*Ang^3)]
Ncell = 3 3 100                 # the size of temperature cells
temperature_mode = double       #single (electron heat conduction) or double (two-temperature model)
heat_mapping_R = 200.;          # the radius of the heat mapping box [Angstrom]
heat_mapping_H = 1000.;         # the height of the heat mapping box [Angstrom]
heat_dt = 40                    # heat calculation time step [fs]

# Field parameters
#---------------
charge_tolerance 0.8 1.2        # min and max tolerance of charge calculation
field_tolerance 0.1 100.0         # min and max allowed deviation of maximum field from the semianalytical one
field_mode = transient          # method to calculate field; laplace or transient
potential_limit = -2 1.e4       # the minimum and maxmum value of potential
field_cgtol = 1e-3              # tolerance of field solver
elfield_mode = DC               # Method to control the applied electric field; DC, AC, pulse, or index E*(1-exp(-t/tau))
elfield = -0.0320000             # DC--value of applied elctric field; pulse--initial value  [V/Angstrom][10000 MV/m]
anode_bc = neumann              # boundary condition type at anode; dirichlet or neumann

# PIC parameters
#---------------
pic_dtmax = 0.5                 # maximum time interval between particle push in PIC [fs]
pic_collide_ee = true           # turn on electron superparticle collisions
electron_weight = 0.01          # electron superparticle weight [nr of SPs]
pic_periodic = true

## Atomic structure unit cell file CONTCAR for building a specified geometry tip:
Cu
1.0
        3.6147000790         0.0000000000         0.0000000000
        0.0000000000         3.6147000790         0.0000000000
        0.0000000000         0.0000000000         3.6147000790
   Cu
    4
Direct
     0.000000000         0.000000000         0.000000000
     0.000000000         0.500000000         0.500000000
     0.500000000         0.000000000         0.500000000
     0.500000000         0.500000000         0.000000000

## Potential function requires input parameter or file format:
(1) The LJ potential parameters
Strength parameter ε
Length parameter σ 
(2) The EAM potential file
Single-element: DYNAMO funcfl format [47]
line 1: comment (ignored)
line 2: atomic number, mass, lattice constant, lattice type (e.g. FCC)
line 3: Nrho, drho, Nr, dr, cutoff
Following the three header lines are three arrays of tabulated values:
embedding function F(rho) (Nrho values)
effective charge function Z(r) (Nr values)
density function rho(r) (Nr values)
Multi-element: DYNAMO setfl format
lines 1,2,3 = comments (ignored)
line 4: Nelements Element1 Element2 … ElementN
line 5: Nrho, drho, Nr, dr, cutoff
Following the 5 header lines are Nelements sections, one for each element, each with the following format:
line 1 = atomic number, mass, lattice constant, lattice type (e.g. FCC)
embedding function F(rho) (Nrho values)
density function rho(r) (Nr values)
(3) MLP
Training: 
MTP initial file of a certain level
Training mode (active training or passive learning)
Atomic structure in cfg format [35]
DFT calculation input file
MTP file format:
species_count is the number of components (atomic types) in the system investigated;
radial_basis_type is the type of radial basis used (the most frequently used radial basis is the Chebyshev basis);
min_dist is the minimal distance between atoms in the training set (in Angstroms);
max_dist is the cutoff radius (in Angstroms);
radial_basis_size is the size of radial basis (the Chebyshev basis for the given MTP).

## output file:
The output files are placed in the FEcMD/out directory.
*.xyz and *.movie files can be visualized by ‘ovito’.
*.vtk files can be visualized by ‘paraview’.
== finite-element mesh generation results
atomreader.ckx: the atomic coordinates read by FECMOS library
surface_dense.xyz: the coordinates of the extracted tip surface points
surface_coarse.xyz: the coordinates of the coarse-grained tip surface points
mesh/*vtk: meshes generated for finite element calculation
== ED results:
ch_solver.movie: current and heat solver results during the evolution
ch_solver.xyz: current and heat solver results currently being simulated
E0.dat: the value of applied electric field during the evolution
electrons.movie: the electron trajectory simulated by the PIC during the evolution
electrons.xyz: the electron trajectory simulated by the PIC currently being simulated
emission.dat: the electron filed emission data
exchange-correlation.movie: exchange-correlation potential during the evolution
exchange-correlation.xyz: exchange-correlation potential currently being simulated
fields.movie: the electric field distribution on the tip during the evolution
fields.xyz: the electric field distribution on the tip currently being simulated
surface_fields.movie: the electric field distribution on the tip surface during the evolution
surface_fields.xyz: the electric field distribution on the tip surface currently being simulated
FJI.dat: the electric field, current density and current data
forces.movie: the force induced by electric field during the evolution
forces.xyz: the force induced by electric field currently being simulated
surface_temperature.movie:  the electron temperature distribution on the tip surface during the evolution
surface_temperature.xyz:  the electron temperature distribution on the tip surface currently being simulated
temperature_phonon.movie: the electron and phonon temperature distribution of the tip during the evolution
temperature_phonon.xyz: the electron and phonon temperature distribution of the tip currently being simulated
temperature.dat: the electron and phonon temperature data.
mesh/
=== MD results: (out/MD)
md.movie: the trajectory file of the tip atoms during the evolution
summary.dat: atomic average energy, average temperature, tip height. etc
T_profile.dat: the temperature profile obtained from MD
V_profile.dat: the atomic velocity profile
rdf.dat:  the radial distribution function
