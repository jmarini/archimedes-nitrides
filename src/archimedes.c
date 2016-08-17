/* archimedes.c -- This belongs to GNU archimedes.

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>
   <jsellier@purdue.edu>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


// =================================================================
// File Name : archimedes.c
// Version   : release 2.0.0
// Date of Creation : 07 Dec.2003, Paris, France, Jean Michel Sellier
// Last Revision : 15 Sept. 2011, Carry le Rouet, France, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

#include <getopt.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#ifdef	HAVE_STRING_H
    #include <string.h>
#else
    #include <strings.h>
#endif

// Preprocessor Definitions
#define real double
#define ON 1
#define OFF 0
#define INVALID -1             // invalid value for most enum-like variables
#define KANE 0                 // conduction band model, kane model
#define PARABOLIC 1            // conduction band model, parabolic approximation
#define FULL 2                 // conduction band model, full band model
#define QEP_BOHM 0             // quantum effective potential, bohm potential
#define QEP_CALIBRATED_BOHM 1  // quantum effective potential, calibrated bohm potential
#define QEP_FULL 2             // quantum effective potential, full effective potential
#define QEP_DENSITY_GRADIENT 3 // quantum effective potential, density gradient
#define MN3 4                  // number of summary values to save per cell
#define NXM 308                // maximum number of cells in x-direction
#define NYM 308                // maximum number of cells in y-direction
#define DIME 3003              // maximum number of points in energy mesh
#define ITMAX 10000000         // maximum number of monte carlo iterations
#define POISSONITMAX 1500      // maximum number of poisson iterations
#define SMALL 1.e-5            // defines what is a "small" number/delta
#define VMAX 1000000
#define NPMAX 10000000         // maximum number of super-particles
#define MCE 0                  // MCE stands for MC for electrons only
#define MCH 1                  // MCH stands for MC for holes only
#define MCEH 2                 // MCEH stands for MC both for electrons and holes
#define MEPE 3                 // MEPE stands for MEP model for electrons only
#define MEPH 4                 // MEPH stands for MEP model for holes only
#define MEPEH 5                // MEPEH stands for MEP model for electrons and holes
#define GNUPLOTFORMAT 0        // output file in GNUPLOT format
#define MESHFORMAT 1           // output file in Mesh format
#define MAX_VALLEYS 4          // maximum number of valleys

// definition of the material reference table
#define NOAMTIA 17   // Number Of All Material Taken Into Account (excluding SiO2)
// ********************
#define iSIO2 -1     // reference number to Silicon Oxide
#define SILICON 0    // reference number to Silicon
#define GAAS 1       // reference number to GaAs
#define GERMANIUM 2  // reference number to Germanium
#define INSB 3       // reference number to InSb
#define ALSB 4       // reference number to AlSb
#define ALXINXSB 5   // reference number to Al_x In_x Sb
#define ALXIN1XSB 6  // reference number to Al_x In_(1-x) Sb
#define ALAS 7       // reference number to AlAs
#define ALP 8        // reference number to AlP
#define GAP 9        // reference number to GaP
#define GASB 10      // reference number to GaSb
#define INAS 11      // reference number to InAs
#define INP 12       // reference number to InP
#define INXGA1XAS 13 // reference number to In_x Ga_(1-x) As
#define INXAL1XAS 14 // reference number to In_x Al_(1-x) As
#define INXGAXXAS 15 // reference number to In_x Ga_(1-x) As (second zone)
#define GAN 16       // reference number to GaN
// ********************

#define NUMSIO2 2 // maximum number of SiO2 interfaces
// ===============================
#include "constants.h"
#include "extrema.h"
#include "sign.h"
#include "mm.h"
#include "mm2.h"
#include "particle.h"
extern inline int mc_does_particle_exist(particle_t *particle);
extern inline void mc_remove_particle(particle_t *particle);
extern inline real mc_particle_ksquared(particle_t *particle);
extern inline real mc_particle_k(particle_t *particle);
// ===============================

// All integers here...
int NUM_VERT;            // number of vertices in the meshing
int NUM_EXAHEDRA;        // number of quadrilaterals in the meshing
int MEDIA;               // number of time steps macroscopic variables will be averaged/computed over, defaults to 500
int MAXIMINI;            // boolean, whether to save max & min values of macroscopic variables during simulation, defaults to 0
int SAVEALWAYS;          // boolean, whether to save information at each step, defaults to 0
int SCATTERING_OUTPUT;   // boolean, whether to output scattering rates, defaults to 0
int nx, ny;              // number of cells in x- and y-directions
int ISEED;               // seed for random number generator, starts at 38467
int NP1;                 // number of particles in n+ cell, defaults to 2500
int INUM;                // number of electrons sumulated
int c;                   // iteration number
int Model_Number;        // enum, controls which simulation model is used, values include MCE, MCH, MCEH, MEPE, MEPH, MEPEH, defaults to MCE
int File_Format;         // enum, controls which file format to output to, values include GNUPLOTFORMAT and MESHFORMAT, defaults to GNUPLOTFORMAT
int Quantum_Flag;        // boolean, controls whether quantum effects are simulated using effective potential method, defaults to 0
int NG;                  // number of mesh nodes
int NE;                  // number of mesh triangles
int leid_flag;           // boolean, controls whether starting point used previously saved results, defaults to 0
int SIO2_UP_FLAG;        // boolean, controls whether SIO2 is above the devide, defaults to 0
int SIO2_DOWN_FLAG;      // boolean, controls whether SIO2 is below the devide, defaults to 0
int FARADAYFLAG;         // boolean, controls whether evolution of magnetic field will be calculated, defaults to 0
int i_dom[NXM+1][NYM+1]; // material at each mesh node, array indexed by node (i, j)
int NOVALLEY[NOAMTIA+1]; // number of valleys to simulate, array indexed by material
int ZSCATTER[NOAMTIA+1][6][6]; // number of equivalent valleys for scattering, array indexed by material, starting valley and ending valley
int ACOUSTICPHONONS;     // boolean, controls whether acoustic phonon scattering is used, defaults to 1
int OPTICALPHONONS;      // boolean, controls whether optical phonon scattering is used, defaults to 1
int IMPURITYPHONONS;     // boolean, controls whether impurity scattering is used, defaults to 1
int PIEZOELECTRIC;       // boolean, controls whether piezoelectric scattering is used, defaults to 0
int CONDUCTION_BAND;     // enum, controls which band structure model is used, values include PARABOLIC, KANE, FULL, defaults to KANE
int SAVE_MESH;           // boolean, controls whether mesh is saved, defaults to 0
int NODE_GEO[3][(NXM+1)*(NYM+1)];  // stores the node geometry, nxm+1 x nym+1 x 3 array
long long int PARTICLE_ID;         // tracker for next particle id

// All "real"'s here...
real moving_average[NXM+1][NYM+1][MN3+1]; // Holds moving average of calculated values, array indexed by mesh node and value type:
                                          //  type = 0: unused
                                          //  type = 1: unused
                                          //  type = 2: particle x-velocity
                                          //  type = 3: particle y-velocity
                                          //  type = 4: particle energy
real moving_alpha;                  // for calculating exponential moving average - represents the degree of weighting decrease for older observations
particle_info_t particle_info[NPMAX+1];   // Holds summary information for particles to be output, array indexed by particle index
real u2d[NXM+1][NYM+1][MN3+1];      // Hold summary values for electrons per cell, array indexed by mesh node and value type:
                                    //  type = 0: quantum effective potential
                                    //  type = 1: electron density
                                    //  type = 2: running sum of electron x-velocity (divide by MEDIA to get average)
                                    //  type = 3: running sum of electron y-velocity (divide by MEDIA to get average)
                                    //  type = 4: running sum of electron energy     (divide by MEDIA to get average)
real h2d[NXM+1][NYM+1][MN3+1];      // Hold summary values for holes per cell, array indexed by mesh node and value type:
                                    //  type = 0: quantum effective potential
                                    //  type = 1: hole density
                                    //  type = 2: running sum of hole x-velocity (divide by MEDIA to get average)
                                    //  type = 3: running sum of hole y-velocity (divide by MEDIA to get average)
                                    //  type = 4: running sum of hole energy     (divide by MEDIA to get average)
real PSI[NXM+1][NYM+1];             // Potential, indexed by mesh node
real E[NXM+1][NYM+1][2];            // E-field, indexed by mesh node
real N_D[NXM+1][NYM+1];             // Donor concentration, indexed by mesh node, defaults to NI
real N_H[NXM+1][NYM+1];             // Acceptor concentration, indexed by mesh node, defaults to NI
real dx, dy;                        // length of cells in x & y directions
real TEMPO=0.;                      // current time in simulation, starts at 0
real TF;                            // final time, defaults to 5e-12
real LX, LY;                        // length of device in x & y directions
real TL;                            // lattice temperature
real DT;                            // time step, defaults to 0.001e-12
real BKTQ;                          // precomputed constant, k * T_lattice / Q [eV]
real QH;                            // precomputed constant, q / hbar
real SMH[NOAMTIA+1][3];             // precomputed constant, sqrt(2 * m* * m_e * q) / hbar, array indexed by material and valley number
real HHM[NOAMTIA+1][3];             // precomputed constant, hbar^2 / (2 * m* * m_e * q), array indexed by material and valley number
real HM[NOAMTIA+1][3];              // precomputed constant, hbar / (m* * m_e), array indexed by material and valley number
real GM[NOAMTIA+1];                 // total scattering rate, Gamma=1/t0, array indexed by material
real SWK[NOAMTIA+1][3][14][DIME+1]; // scattering rate, indexed by material, valley, phonon mode/scattering type, energy step (i*DE)
particle_t P[NPMAX+1];              // particle information, array indexed by particle
real EPP;                           // number of carriers per particle (?)
real DDmax;                         // maximum donor density (?)
real EDGE[4][NXM+NYM+1][4];         // stores information on edges, array indexed by edge type (0=bottom, 1=right, 2=top, 3=left),
                                    //                                               cell index (i or j),
                                    //                                               information type (0=boundary type (0=insulator, 1=schottky, 2=ohmic),
                                    //                                                                 1=potential,
                                    //                                                                 2=contact electron density,
                                    //                                                                 3=contact hole density)
real CIMP;                          // impurity concentration
real QD2;                           // precomputed constant, qd^2, qd=sqrt(q * cimp / ktq / eps)
real TAUW;                          // energy relaxation time, defaults to 0.4e-12
real EPSRSIO2;
real bufx2d[NXM+1][NYM+1];          // MEP
real bufy2d[NXM+1][NYM+1];          // MEP
real ux2d[NXM+1][NYM+1][MN3+1];     // MEP
real uy2d[NXM+1][NYM+1][MN3+1];     // MEP
real f2d[NXM+1][NYM+1][MN3+1];      // MEP
real g2d[NXM+1][NYM+1][MN3+1];      // MEP
real fx2d[NXM+1][NYM+1][MN3+1];     // MEP
real gy2d[NXM+1][NYM+1][MN3+1];     // MEP
real c11[7],c12[7],c21[7],c22[7];   // MEP
real u[7],f[7],g[7],cw[7];          // MEP
real SIO2_INI[NUMSIO2];
real SIO2_FIN[NUMSIO2];
real SIO2_POT[NUMSIO2];
real SIO2_THICKNESS[NUMSIO2];
real SIO2[NUMSIO2][NXM+1][NYM+1];
real B[NXM+1][NYM+1];         // magnetic field, indexed by mesh node
real EPSR[NOAMTIA+1];         // static dielectric constant, array indexed by material
real EPF[NOAMTIA+1];          // high frequency dielectric constant, array indexed by material
real MSTAR[NOAMTIA+1][6];     // effective mass, array indexed by material and by valley number
real MSTAR_VB[NOAMTIA+1][3];  // effective valence band mass, array indexed by material and by valley number
                              //   valley = 0: heavy hole
                              //   valley = 1: light hole
                              //   valley = 2: split-off
real DELTAE_VB[NOAMTIA+1][3]; // energy difference between VBM and valence band, array indexed by material and by valley number
                              //   valley = 0: heavy hole
                              //   valley = 1: light hole
                              //   valley = 2: split-off
real alphaK[NOAMTIA+1][4];    // valley non-parabolicity, array indexed by material annd by valley number
real EG[NOAMTIA+1];           // band gap, array indexed by material
real HWO[NOAMTIA+1][6];       // optical phonon scattering energy, array indexed by material and up to 6 different values
real DTK[NOAMTIA+1][6];       // optical coupling constant, array indexed by material and by up to 6 different values
real ZF[NOAMTIA+1][6];        // optical phonon z-factor, array indexed by material and by up to 6 different values
real RHO[NOAMTIA+1];          // crystal density, array indexed by material
real DA[NOAMTIA+1];           // acoustic deformation potential, array indexed by material
real UL[NOAMTIA+1];           // longitudinal sound velocity, array indexed by material
real EMIN[NOAMTIA+1][4];      // energy difference between valley minimum and CB minimum, array indexed by material and by valley number
real XVAL[NOAMTIA+1];         // x-mole fraction, array indexed by material
real LATTCONST[NOAMTIA+1];    // lattice constant, array indexed by material
real CB_FULL[NOAMTIA+1][11];  // polynomial coefficients (up to 9th order) for full band structure, array indexed by material
real KAV[NOAMTIA+1];          // electro-mechanical coupling constants, array indexed by material
real QEP_ALPHA;
real QEP_GAMMA;
real QEP_MODEL;
real COORD[2][(NXM+1)*(NYM+1)];

// All structures here...
time_t binarytime;
struct tm *nowtm;
struct option longopts[] =
{
  { "version", no_argument, NULL, 'v' },
  { "help", no_argument, NULL, 'h' }
};
// All files here...
FILE *fp;
// All strings here...
static char *progname;

#include "utility.h"
#include "mesher.h"
#include "poissonbcs.h"
#include "faradaybcs.h"
#include "media.h"
#include "saveoutput2dmeshformat.h"
#include "saveoutput2dgnuplot.h"
#include "saveoutput2dholegnuplot.h"
#include "saveoutput2dholemeshformat.h"
#include "saveoutputfiles.h"
#include "quantumeffectivepotential.h"
#include "electric_field.h"
#include "faraday.h"
#include "random.h"
#include "scattering_rates.h"
#include "optical_absorption.h"
#include "deviceconfig.h"
#include "particlecreation.h"
#include "drift.h"
#include "scattering.h"
#include "ensemblemontecarlo.h"
#include "charge.h"
#include "computecurrents.h"
#include "MEP_interpolation.h"
#include "electron_relaxation.h"
#include "HMEPbcs.h"
#include "ParabMEP2D.h"
#include "Hole_bcs.h"
#include "holemep2d.h"
#include "hole_relaxation.h"
#include "updating.h"
#include "readinputfile.h"
//#include "SaveRappture.h"

// provide extern declarations of functions to fix compiler error
extern inline real rnd(void);
extern inline particle_t creation(int i, real t, int edge);
extern inline real MM(real a, real b);
extern inline real MM2(real x, real a, real b);
extern inline real sign(real a, real b);
extern inline real minimus(real x, real y);
extern inline real maximus(real x, real y);
extern inline int mc_is_boundary_insulator(int direction, int index);
extern inline int mc_is_boundary_schottky(int direction, int index);
extern inline int mc_is_boundary_ohmic(int direction, int index);
extern inline int mc_is_boundary_contact(int direction, int index);
extern inline void mc_particle_coords(particle_t *particle, int *i, int *j);
extern inline char* mc_material_name(int material);
extern inline char* mc_band_model_name(int model);
particle_info_t mc_calculate_particle_info(particle_t *p);


int main(int argc, char *argv[]) {

    int optc;
    int h = 0,
        v = 0,
        z = 0,
        lose = 0;
    progname = argv[0];

    while((optc = getopt_long(argc, argv, "hv", longopts, (int *) 0)) != EOF) {
        switch (optc) {
            case 'v':
                v = 1;
                break;
            case 'h':
                h = 1;
                break;
            default:
                lose = 1;
                break;
        }
    }

    if(optind == argc - 1) {
        z = 1;
    }
    else if(lose || optind < argc) {
        /* Print error message and exit.  */
        if(optind < argc) {
            printf("Too many arguments\n");
        }
        printf("Try `%s --help' for more information.\n",progname);
        exit(1);
    }

    /* `help' should come first.  If `help' is requested, ignore the other
       options. */
    if(h) {
        /* Print help info and exit. */
        printf("GNU archimedes, a simulator for submicron and nanoscale "
               "semiconductor devices.\n"
               "Copyright (C) 2004-2011 Jean Michel D. Sellier.\n\n");
        printf ("Usage: %s [OPTION] file...\n\n", progname);

        printf("-h, --help          display this help and exit\n"
               "-v, --version       display version information and exit\n"
               "\n");

        printf ("Report bugs to jeanmichel.sellier@gmail.com "
                "or jsellier@purdue.edu\n");
        exit(EXIT_SUCCESS);
    }

    if(v) {
        /* Print version number and exit. */
        printf("archimedes - GNU archimedes 2.0.1\n\n");
        printf("Copyright (C) 2004 - 2011 Jean Michel D. Sellier.\n\n"
               "There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A\n"
               "PARTICULAR PURPOSE.\n"
               "You may redistribute copies of GNU Archimedes under the terms\n"
               "of the GNU General Public License.\n"
               "For more information about these matters, see the file named COPYING.\n");
        exit(EXIT_SUCCESS);
    }

    if(!z) {
        // No filename has been specified
        // ==============================
        printf("%s: no input file\n", progname);
        exit(EXIT_FAILURE);
    }

    // In case of filename specified
    // =============================
    fp = fopen(argv[1], "r"); // here we open the input file...
    // File Control, just in case the file does not exist...
    if(fp == NULL) {
        printf("%s: fatal error in opening the input file %s\n",
               progname, argv[1]);
        exit(EXIT_FAILURE);
    }

    // material constants definition
    #include "material_parameters.h"


    // We reset the some arrays
    // ========================
    memset(&u2d, 0, sizeof(u2d));
    memset(&E, 0, sizeof(E));
    memset(&EDGE, 0, sizeof(EDGE));
    memset(&SIO2, 0, sizeof(SIO2));

    // Read the geometrical and physical description of the MESFET
    // ===========================================================
    Read_Input_File();
    // ===========================================================

    // Construction of the mesh for the electrostatic potential
    // (to properly take into account the oxyde layers)
    mesher();

    // III-V Semiconductor materials energy gap (depending on the lattice temperature)
    printf("\n");
    EG[SILICON]=1.21-3.333e-4*TL;
    printf("EG[SILICON]      = %g\n",EG[SILICON]);
    EG[GERMANIUM]=0.747-3.587e-4*TL;
    printf("EG[GERMANIUM]    = %g\n",EG[GERMANIUM]);
    EG[GAAS]=1.54-4.036e-4*TL;
    printf("EG[GAAS]         = %g\n",EG[GAAS]);
    EG[INSB]=0.2446-2.153e-4*TL;
    printf("EG[INSB]         = %g\n",EG[INSB]);
    EG[ALSB]=1.696-2.20e-4*TL;
    printf("EG[ALSB]         = %g\n",EG[ALSB]);
    EG[ALAS]=2.314-3.0e-4*TL;
    printf("EG[ALAS]         = %g\n",EG[ALAS]);
    EG[ALP]=2.51-3.333e-4*TL;
    printf("EG[ALP]         = %g\n",EG[ALP]);
    EG[GAP]=2.35-2.667e-4*TL;
    printf("EG[GAP]         = %g\n",EG[GAP]);
    EG[GASB]=0.81-3.667e-4*TL;
    printf("EG[GASB]         = %g\n",EG[GASB]);
    EG[INAS]=0.434-2.601e-4*TL;
    printf("EG[INAS]         = %g\n",EG[INAS]);
    EG[INP]=1.445-3.296e-4*TL;
    printf("EG[INP]         = %g\n",EG[INP]);
    EG[GAN]=3.47 - 7.7e-4 * TL * TL / (TL + 600);
    printf("EG[GAN]         = %g\n", EG[GAN]);
    printf("\n");

    if(CONDUCTION_BAND == KANE ||
       CONDUCTION_BAND == PARABOLIC ||
       CONDUCTION_BAND == FULL) {
        // USED WHATEVER IS THE CONDUCTION BAND FOR THE INITIAL PSUEDO WAVE VECTOR
        // OF THE PSEUDO PARTICLES
        // all the following non-parabolicity coefficients depend on lattice temperature
        // non-parabolicity coefficient for GaAs in the GAMMA-valley
        alphaK[GAAS][1]=pow(1.-MSTAR[GAAS][1],2.)/(EG[GAAS]+EMIN[GAAS][1]);//expected value = 0.611
        printf("alphaK_gamma[GaAs] = %g\n",alphaK[GAAS][1]);
        // non-parabolicity coefficient for GaAs in the L-valley
        alphaK[GAAS][2]=pow(1.-MSTAR[GAAS][2],2.)/(EG[GAAS]+EMIN[GAAS][2]);//expected value = 0.242;
        printf("alphaK_L[GaAs]     = %g\n",alphaK[GAAS][2]);
        // non-parabolicity coefficient for InSb in the GAMMA-valley
        alphaK[INSB][1]=pow(1.-MSTAR[INSB][1],2.)/(EG[INSB]+EMIN[INSB][1]);//5.59;
        printf("alphaK_gamma[InSb] = %g\n",alphaK[INSB][1]);
        // non-parabolicity coefficient for AlSb in the GAMMA-valley
        alphaK[ALSB][1]=pow(1.-MSTAR[ALSB][1],2.)/(EG[ALSB]+EMIN[ALSB][1]);//0.321;
        printf("alphaK_gamma[AlSb] = %g\n",alphaK[ALSB][1]);
        // non-parabolicity coefficient for AlAs in the GAMMA-valley
        alphaK[ALAS][1]=pow(1.-MSTAR[ALAS][1],2.)/(EG[ALAS]+EMIN[ALAS][1]);
        printf("alphaK_gamma[AlAs] = %g\n",alphaK[ALAS][1]);
        // non-parabolicity coefficient for AlP in the GAMMA-valley
        alphaK[ALP][1]=pow(1.-MSTAR[ALP][1],2.)/(EG[ALP]+EMIN[ALP][1]);
        printf("alphaK_gamma[AlP] = %g\n",alphaK[ALP][1]);
        // non-parabolicity coefficient for GaP in the GAMMA-valley
        alphaK[GAP][1]=pow(1.-MSTAR[GAP][1],2.)/(EG[GAP]+EMIN[GAP][1]);
        printf("alphaK_gamma[GaP] = %g\n",alphaK[GAP][1]);
        // non-parabolicity coefficient for GaSb in the GAMMA-valley
        alphaK[GASB][1]=pow(1.-MSTAR[GASB][1],2.)/(EG[GASB]+EMIN[GASB][1]);
        printf("alphaK_gamma[GaSb] = %g\n",alphaK[GASB][1]);
        // non-parabolicity coefficient for InAs in the GAMMA-valley
        alphaK[INAS][1]=pow(1.-MSTAR[INAS][1],2.)/(EG[INAS]+EMIN[INAS][1]);
        printf("alphaK_gamma[InAs] = %g\n",alphaK[INAS][1]);
        // non-parabolicity coefficient for InP in the GAMMA-valley
        alphaK[INP][1]=pow(1.-MSTAR[INP][1],2.)/(EG[INP]+EMIN[INP][1]);
        printf("alphaK_gamma[InP] = %g\n",alphaK[INP][1]);

        alphaK[GAN][1] = pow(1. - MSTAR[GAN][1], 2.) / (EG[GAN] + EMIN[GAN][1]); // expected value = 0.189
        printf("alphaK_gamma1[GAN] = %g\n", alphaK[GAN][1]);
    }

    // Semiconductor compounds
    // ***
    // Relative dielectric constant for semiconductor compounds
    EPSR[ALXINXSB]=XVAL[ALXINXSB]*EPSR[ALSB]+XVAL[ALXINXSB]*EPSR[INSB];
    EPSR[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPSR[ALSB]+(1.-XVAL[ALXIN1XSB])*EPSR[INSB];
    EPSR[INXGA1XAS]=XVAL[INXGA1XAS]*EPSR[INAS]+(1.-XVAL[INXGA1XAS])*EPSR[GAAS];
    EPSR[INXAL1XAS]=XVAL[INXAL1XAS]*EPSR[INAS]+(1.-XVAL[INXAL1XAS])*EPSR[ALAS];
    EPSR[INXGAXXAS]=XVAL[INXGAXXAS]*EPSR[INAS]+XVAL[INXGAXXAS]*EPSR[GAAS];
    // semiconductor compounds high frequency dieletric constant
    EPF[ALXINXSB]=XVAL[ALXINXSB]*(EPF[ALSB]+EPF[INSB]);
    EPF[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPF[ALSB]+(1.-XVAL[ALXIN1XSB])*EPF[INSB];
    EPF[INXGA1XAS]=XVAL[INXGA1XAS]*EPF[INAS]+(1.-XVAL[INXGA1XAS])*EPF[GAAS];
    EPF[INXAL1XAS]=XVAL[INXAL1XAS]*EPF[INAS]+(1.-XVAL[INXAL1XAS])*EPF[ALAS];
    EPF[INXGAXXAS]=XVAL[INXGAXXAS]*EPF[INAS]+XVAL[INXGAXXAS]*EPF[GAAS];
    // semiconductor compounds optical phonon scattering energy (eV)
    HWO[ALXINXSB][0]=XVAL[ALXINXSB]*(HWO[ALSB][0]+HWO[INSB][0]);
    HWO[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*HWO[ALSB][0]+(1.-XVAL[ALXIN1XSB])*HWO[INSB][0];
    HWO[INXGA1XAS][0]=XVAL[INXGA1XAS]*HWO[INAS][0]+(1.-XVAL[INXGA1XAS])*HWO[GAAS][0];
    HWO[INXAL1XAS][0]=XVAL[INXAL1XAS]*HWO[INAS][0]+(1.-XVAL[INXAL1XAS])*HWO[ALAS][0];
    HWO[INXGAXXAS][0]=XVAL[INXGAXXAS]*HWO[INAS][0]+XVAL[INXGAXXAS]*HWO[GAAS][0];
    // semiconductor compounds optical coupling constants (eV/m)
    DTK[ALXINXSB][0]=XVAL[ALXINXSB]*(DTK[ALSB][0]+DTK[INSB][0]);
    DTK[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*DTK[ALSB][0]+(1.-XVAL[ALXIN1XSB])*DTK[INSB][0];
    DTK[INXGA1XAS][0]=XVAL[INXGA1XAS]*DTK[INAS][0]+(1.-XVAL[INXGA1XAS])*DTK[GAAS][0];
    DTK[INXAL1XAS][0]=XVAL[INXAL1XAS]*DTK[INAS][0]+(1.-XVAL[INXAL1XAS])*DTK[ALAS][0];
    DTK[INXGAXXAS][0]=XVAL[INXGAXXAS]*DTK[INAS][0]+XVAL[INXGAXXAS]*DTK[GAAS][0];
    // semiconductor compounds optical phonon Z-factor
    ZF[ALXINXSB][0]=XVAL[ALXINXSB]*(ZF[ALSB][0]+ZF[INSB][0]);
    ZF[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*ZF[ALSB][0]+(1.-XVAL[ALXIN1XSB])*ZF[INSB][0];
    ZF[INXGA1XAS][0]=XVAL[INXGA1XAS]*ZF[INAS][0]+(1.-XVAL[INXGA1XAS])*ZF[GAAS][0];
    ZF[INXAL1XAS][0]=XVAL[INXAL1XAS]*ZF[INAS][0]+(1.-XVAL[INXAL1XAS])*ZF[ALAS][0];
    ZF[INXGAXXAS][0]=XVAL[INXGAXXAS]*ZF[INAS][0]+XVAL[INXGAXXAS]*ZF[GAAS][0];
    // semiconductor compounds Crystal Density (Kg/m^3)
    RHO[ALXINXSB]=XVAL[ALXINXSB]*(RHO[ALSB]+RHO[INSB]);
    RHO[ALXIN1XSB]=XVAL[ALXIN1XSB]*RHO[ALSB]+(1.-XVAL[ALXIN1XSB])*RHO[INSB];
    RHO[INXGA1XAS]=XVAL[INXGA1XAS]*RHO[INAS]+(1.-XVAL[INXGA1XAS])*RHO[GAAS];
    RHO[INXAL1XAS]=XVAL[INXAL1XAS]*RHO[INAS]+(1.-XVAL[INXAL1XAS])*RHO[ALAS];
    RHO[INXGAXXAS]=XVAL[INXGAXXAS]*RHO[INAS]+XVAL[INXGAXXAS]*RHO[GAAS];
    // semiconductor compounds acoustic deformation potential (Joule)
    DA[ALXINXSB]=XVAL[ALXINXSB]*(DA[ALSB]+DA[INSB]);
    DA[ALXIN1XSB]=XVAL[ALXIN1XSB]*DA[ALSB]+(1.-XVAL[ALXIN1XSB])*DA[INSB];
    DA[INXGA1XAS]=XVAL[INXGA1XAS]*DA[INAS]+(1.-XVAL[INXGA1XAS])*DA[GAAS];
    DA[INXAL1XAS]=XVAL[INXAL1XAS]*DA[INAS]+(1.-XVAL[INXAL1XAS])*DA[ALAS];
    DA[INXGAXXAS]=XVAL[INXGAXXAS]*DA[INAS]+XVAL[INXGAXXAS]*DA[GAAS];
    // semiconductor compounds longitudinal sound velocity (m/sec)
    UL[ALXINXSB]=XVAL[ALXINXSB]*(UL[ALSB]+UL[INSB]);
    UL[ALXIN1XSB]=XVAL[ALXIN1XSB]*UL[ALSB]+(1.-XVAL[ALXIN1XSB])*UL[INSB];
    UL[INXGA1XAS]=XVAL[INXGA1XAS]*UL[INAS]+(1.-XVAL[INXGA1XAS])*UL[GAAS];
    UL[INXAL1XAS]=XVAL[INXAL1XAS]*UL[INAS]+(1.-XVAL[INXAL1XAS])*UL[ALAS];
    UL[INXGAXXAS]=XVAL[INXGAXXAS]*UL[INAS]+XVAL[INXGAXXAS]*UL[GAAS];
    // semiconductor compounds energy gap
    EG[ALXINXSB]=XVAL[ALXINXSB]*(EG[ALSB]+EG[INSB]);
    EG[ALXIN1XSB]=XVAL[ALXIN1XSB]*EG[ALSB]+(1.-XVAL[ALXIN1XSB])*EG[INSB];
    EG[INXGA1XAS]=XVAL[INXGA1XAS]*EG[INAS]+(1.-XVAL[INXGA1XAS])*EG[GAAS];
    EG[INXAL1XAS]=XVAL[INXAL1XAS]*EG[INAS]+(1.-XVAL[INXAL1XAS])*EG[ALAS];
    EG[INXGAXXAS]=XVAL[INXGAXXAS]*EG[INAS]+XVAL[INXGAXXAS]*EG[GAAS];
    // semiconductor compounds energy minimum of GAMMA-valley
    EMIN[ALXINXSB][1]=XVAL[ALXINXSB]*(EMIN[ALSB][1]+EMIN[INSB][1]);
    EMIN[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*EMIN[ALSB][1]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][1];
    EMIN[INXGA1XAS][1]=XVAL[INXGA1XAS]*EMIN[INAS][1]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][1];
    EMIN[INXAL1XAS][1]=XVAL[INXAL1XAS]*EMIN[INAS][1]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][1];
    EMIN[INXGAXXAS][1]=XVAL[INXGAXXAS]*EMIN[INAS][1]+XVAL[INXGAXXAS]*EMIN[GAAS][1];
    // semiconductor compounds energy minimum 0f L-valley
    EMIN[ALXINXSB][2]=XVAL[ALXINXSB]*(EMIN[ALSB][2]+EMIN[INSB][2]);
    EMIN[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*EMIN[ALSB][2]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][2];
    EMIN[INXGA1XAS][2]=XVAL[INXGA1XAS]*EMIN[INAS][2]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][2];
    EMIN[INXAL1XAS][2]=XVAL[INXAL1XAS]*EMIN[INAS][2]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][2];
    EMIN[INXGAXXAS][2]=XVAL[INXGAXXAS]*EMIN[INAS][2]+XVAL[INXGAXXAS]*EMIN[GAAS][2];
    // GAMMA-valley effective mass
    MSTAR[ALXINXSB][1]=XVAL[ALXINXSB]*(MSTAR[ALSB][1]+MSTAR[INSB][1]);
    MSTAR[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*MSTAR[ALSB][1]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][1];
    MSTAR[INXGA1XAS][1]=XVAL[INXGA1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][1];
    MSTAR[INXAL1XAS][1]=XVAL[INXAL1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][1];
    MSTAR[INXGAXXAS][1]=XVAL[INXGAXXAS]*MSTAR[INAS][1]+XVAL[INXGAXXAS]*MSTAR[GAAS][1];
    // L-valley effective mass
    MSTAR[ALXINXSB][2]=XVAL[ALXINXSB]*(MSTAR[ALSB][2]+MSTAR[INSB][2]);
    MSTAR[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*MSTAR[ALSB][2]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][2];
    MSTAR[INXGA1XAS][2]=XVAL[INXGA1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][2];
    MSTAR[INXAL1XAS][2]=XVAL[INXAL1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][2];
    MSTAR[INXGAXXAS][2]=XVAL[INXGAXXAS]*MSTAR[INAS][2]+XVAL[INXGAXXAS]*MSTAR[GAAS][2];
    // non-parabolicity coefficient for semiconductor compounds in the GAMMA-valley
    alphaK[ALXINXSB][1]=XVAL[ALXINXSB]*(alphaK[ALSB][1]+alphaK[INSB][1]);
    alphaK[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*alphaK[ALSB][1]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][1];
    alphaK[INXGA1XAS][1]=XVAL[INXGA1XAS]*alphaK[INAS][1]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][1];
    alphaK[INXAL1XAS][1]=XVAL[INXAL1XAS]*alphaK[INAS][1]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][1];
    alphaK[INXGAXXAS][1]=XVAL[INXGAXXAS]*alphaK[INAS][1]+XVAL[INXGAXXAS]*alphaK[GAAS][1];
    // non-parabolicity coefficient for Al_x In_(1-x) Sb in the L-valley
    alphaK[ALXINXSB][2]=XVAL[ALXINXSB]*(alphaK[ALSB][2]+alphaK[INSB][2]);
    alphaK[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*alphaK[ALSB][2]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][2];
    alphaK[INXGA1XAS][2]=XVAL[INXGA1XAS]*alphaK[INAS][2]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][2];
    alphaK[INXAL1XAS][2]=XVAL[INXAL1XAS]*alphaK[INAS][2]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][2];
    alphaK[INXGAXXAS][2]=XVAL[INXGAXXAS]*alphaK[INAS][2]+XVAL[INXGAXXAS]*alphaK[GAAS][2];
    // ***

    // Closure of the input file
    // =========================
    fclose(fp);
    printf("\nInput file read...\n");

    // Read all the coefficients for MEP simulation
    // ============================================
    MEP_coefficients();

    // Loading of the initial runtime
    // ==============================
    binarytime = time(NULL);
    nowtm = localtime(&binarytime);

    printf("\n\nComputation Started at %s\n", asctime(nowtm));

    // Boundary conditions for the model simulated
    // ===========================================
    PoissonBCs();
    if(FARADAYFLAG) {
        FaradayBCs();
    }
    printf("Boundary conditions calculated...\n");

    // Initialization for Monte Carlo
    // ==============================
    if(Model_Number == MCE || Model_Number == MCEH) {
        int i;
        for(i = 0; i < NOAMTIA; i++) {
            calc_scattering_rates(i);
        }
        printf("Scattering rates calculated...\n");
        MCdevice_config( );
        printf("Device configuration complete...\n");
    }
    printf("\n");


    // HERE IS THE SIMULATION
    // ======================
    for(c = 1; c <= ITMAX; c++) {
        updating(Model_Number);
    }

    // Here we save the outputs
    // ========================
    SaveOutputFiles(File_Format, 0);
    printf("\nFinal Output has been saved\n");

    binarytime=time(NULL);
    nowtm=localtime(&binarytime);
    printf("Computation Finished at %s\n",asctime(nowtm));
    return(EXIT_SUCCESS);
}
