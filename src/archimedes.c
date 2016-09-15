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
#include "global_defines.h"
// ===============================
#include "configuration.h"
#include "mesh.h"
#include "constants.h"
#include "mep/constants.h"
#include "mep/extrema.h"
#include "mep/sign.h"
#include "mep/mm.h"
#include "mep/mm2.h"
#include "particle.h"
#include "material.h"

Configuration *g_config;
Mesh *g_mesh;
Material g_materials[NOAMTIA];

extern inline int mc_does_particle_exist(Particle *particle);
extern inline void mc_remove_particle(Particle *particle);
extern inline real mc_particle_ksquared(Particle *particle);
extern inline real mc_particle_k(Particle *particle);
// ===============================


// All integers here...
int NUM_VERT;            // number of vertices in the meshing
int NUM_EXAHEDRA;        // number of quadrilaterals in the meshing
int ISEED;               // seed for random number generator, starts at 38467
int c;                   // iteration number
int SIO2_UP_FLAG;        // boolean, controls whether SIO2 is above the devide, defaults to 0
int SIO2_DOWN_FLAG;      // boolean, controls whether SIO2 is below the devide, defaults to 0
int NOVALLEY[NOAMTIA+1]; // number of valleys to simulate, array indexed by material
int ZSCATTER[NOAMTIA+1][6][6]; // number of equivalent valleys for scattering, array indexed by material, starting valley and ending valley

// All "real"'s here...
real moving_average[NXM+1][NYM+1][MN3+1]; // Holds moving average of calculated values, array indexed by mesh node and value type:
                                          //  type = 0: unused
                                          //  type = 1: unused
                                          //  type = 2: particle x-velocity
                                          //  type = 3: particle y-velocity
                                          //  type = 4: particle energy
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
real BKTQ;                          // precomputed constant, k * T_lattice / Q [eV]
real QH;                            // precomputed constant, q / hbar
real SMH[NOAMTIA+1][3];             // precomputed constant, sqrt(2 * m* * m_e * q) / hbar, array indexed by material and valley number
real HHM[NOAMTIA+1][3];             // precomputed constant, hbar^2 / (2 * m* * m_e * q), array indexed by material and valley number
real HM[NOAMTIA+1][3];              // precomputed constant, hbar / (m* * m_e), array indexed by material and valley number
real GM[NOAMTIA+1];                 // total scattering rate, Gamma=1/t0, array indexed by material
real SWK[NOAMTIA+1][3][14][DIME+1]; // scattering rate, indexed by material, valley, phonon mode/scattering type, energy step (i*DE)
Particle P[NPMAX+1];              // particle information, array indexed by particle
real EDGE[4][NXM+NYM+1][4];         // stores information on edges, array indexed by edge type (0=bottom, 1=right, 2=top, 3=left),
                                    //                                               cell index (i or j),
                                    //                                               information type (0=boundary type (0=insulator, 1=schottky, 2=ohmic),
                                    //                                                                 1=potential,
                                    //                                                                 2=contact electron density,
                                    //                                                                 3=contact hole density)
real QD2;                           // precomputed constant, qd^2, qd=sqrt(q * cimp / ktq / eps)
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

#include "mep/HMEPbcs.h"
#include "mep/Hole_bcs.h"
#include "mep/MEP_interpolation.h"
#include "mep/ParabMEP2D.h"
#include "mep/electron_relaxation.h"
#include "mep/hole_relaxation.h"
#include "mep/holemep2d.h"

#include "utility.h"
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
#include "particles_per_cell.h"
#include "computecurrents.h"
#include "readinputfile.h"
#include "updating.h"

// provide extern declarations of functions to fix compiler error
extern inline Particle creation(int i, real t, int edge);
extern inline real MM(real a, real b);
extern inline real MM2(real x, real a, real b);
extern inline real sign(real a, real b);
extern inline real minimus(real x, real y);
extern inline real maximus(real x, real y);
extern inline int mc_is_boundary_insulator(int direction, int index);
extern inline int mc_is_boundary_schottky(int direction, int index);
extern inline int mc_is_boundary_ohmic(int direction, int index);
extern inline int mc_is_boundary_contact(int direction, int index);
extern inline char* mc_material_name(int material);
extern inline char* mc_band_model_name(int model);
particle_info_t mc_calculate_particle_info(Particle *p);


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

    g_config = malloc(sizeof *g_config);
    g_mesh = malloc(sizeof *g_mesh);

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
    if(mc_build_mesh(g_mesh) != 0) {
        printf("Error: Unexpected error while building mesh.\n");
        exit(EXIT_FAILURE);
    }
    if(g_config->save_mesh) {
        if(mc_save_mesh(g_mesh, "device.mesh") != 0) {
            printf("Error: Unexpected error while saving mesh.\n");
            exit(EXIT_FAILURE);
        }
    }

    // Closure of the input file
    // =========================
    fclose(fp);
    printf("\nInput file read...\n");

    // material constants definition
    #include "material_parameters.h"

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
    if(poisson_boundary_conditions( ) != 0) {
        printf("Error: Unknown error calculating Poisson boundary conditions.\n");
        return 1;
    }
    if(g_config->faraday_flag) {
        FaradayBCs();
    }
    printf("Boundary conditions calculated...\n");

    // Initialization for Monte Carlo
    // ==============================
    double transistion_rate[NOAMTIA][DIME][3];
    if(g_config->simulation_model == MCE || g_config->simulation_model == MCEH) {
        for(int i = 0; i < NOAMTIA; i++) {
            calc_scattering_rates(i);
            calc_absorption_rates(g_materials[i], transistion_rate);
        }
        printf("Scattering rates calculated...\n");
        MCdevice_config( );
        // photoexcite_carriers(g_mesh, 4.0, transistion_rate, GM, P);
        printf("Device configuration complete...\n");
    }
    printf("\n");


    // FILE *fp = fopen("photoexcited_particles.csv", "w");
    // fprintf(fp, "index x y kx ky kz energy\n");
    // for(int n = 0; n < g_config->num_particles; ++n) {
    //     Particle *p = &P[n];
    //     fprintf(fp, "%d %g %g %g %g %g %g\n", n, p->x, p->y, p->kx, p->ky, p->kz, mc_particle_energy(p));
    // }
    // fclose(fp);


    // HERE IS THE SIMULATION
    // ======================
    for(c = 1; c <= ITMAX; c++) {
        updating(g_config->simulation_model);
    }

    // Here we save the outputs
    // ========================
    SaveOutputFiles(g_config->output_format, 0);
    printf("\nFinal Output has been saved\n");

    binarytime=time(NULL);
    nowtm=localtime(&binarytime);
    printf("Computation Finished at %s\n",asctime(nowtm));
    return(EXIT_SUCCESS);
}
