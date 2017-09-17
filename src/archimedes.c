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

#include "global_defines.h"
#include "configuration.h"
#include "mesh.h"
#include "constants.h"
#include "particle.h"
#include "material.h"

// Extern variables
Configuration *g_config;
Mesh *g_mesh;
Material g_materials[NOAMTIA];
Direction direction_t = {.BOTTOM=0, .RIGHT=1, .TOP=2, .LEFT=3};
Boundary boundary_t = {.INSULATOR=0, .SCHOTTKY=1, .OHMIC=2, .VACUUM=3};


real moving_average[NXM+1][NYM+1][MN3+1]; // Holds moving average of calculated values, array indexed by mesh node and value type:
                                          //  type = 0: unused
                                          //  type = 1: unused
                                          //  type = 2: particle x-velocity
                                          //  type = 3: particle y-velocity
                                          //  type = 4: particle energy
real BKTQ;                          // precomputed constant, k * T_lattice / Q [eV]
real GM[NOAMTIA+1];                 // total scattering rate, Gamma=1/t0, array indexed by material
real SWK[NOAMTIA+1][4][14][DIME+1]; // scattering rate, indexed by material, valley, phonon mode/scattering type, energy step (i*DE)
real QD2;                           // precomputed constant, qd^2, qd=sqrt(q * cimp / ktq / eps)
real XVAL[NOAMTIA+1];         // x-mole fraction, array indexed by material
real CB_FULL[NOAMTIA+1][11];  // polynomial coefficients (up to 9th order) for full band structure, array indexed by material

FILE *input_fp;
FILE *emitted_fp;
FILE *tracking_fp;
FILE *valley_occupation_fp;
FILE *velocity_fp;

static char *progname;


#include "mep/mep.h"

#include "utility.h"
#include "electrostatics.h"
#include "media.h"
#include "saveoutput2dmeshformat.h"
#include "saveoutput2dgnuplot.h"
#include "saveoutput2dholegnuplot.h"
#include "saveoutput2dholemeshformat.h"
#include "saveoutputfiles.h"
#include "random.h"
#include "scattering_rates.h"
#include "optical_absorption.h"
#include "particle_creation.h"
#include "drift.h"
#include "scattering.h"
#include "ensemblemontecarlo.h"
#include "particles_per_cell.h"
#include "computecurrents.h"
#include "readinputfile.h"
#include "updating.h"

// extern declarations of functions
extern inline int mc_does_particle_exist(Particle *particle);
extern inline void mc_remove_particle(Particle *particle);
extern inline char* mc_band_model_name(int model);


int main(int argc, char *argv[]) {

    int optc;
    int h = 0,
        v = 0,
        z = 0,
        lose = 0;
    progname = argv[0];

    struct option longopts[] = {
        {"version", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
    };

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
    input_fp = fopen(argv[1], "r"); // here we open the input file...
    // File Control, just in case the file does not exist...
    if(input_fp == NULL) {
        printf("%s: fatal error in opening the input file %s\n",
               progname, argv[1]);
        exit(EXIT_FAILURE);
    }

    g_config = malloc(sizeof *g_config);
    g_mesh = malloc(sizeof *g_mesh);

    // Read the geometrical and physical description of the MESFET
    // ===========================================================
    read_input_file(input_fp);
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
    fclose(input_fp);
    printf("\nInput file read...\n");

    // material constants definition
    #include "material_parameters.h"

    // Read all the coefficients for MEP simulation
    // ============================================
    MEP_coefficients();

    // Loading of the initial runtime
    // ==============================
    time_t binarytime = time(NULL);
    struct tm *nowtm = localtime(&binarytime);

    printf("\n\nComputation Started at %s\n", asctime(nowtm));

    if(g_config->surface_bb_flag == ON) {
        for(int i = 1; i <= g_mesh->nx + 1; ++i) {
            for(int j = 1; j <= g_mesh->ny + 1; ++j) {
                Node *node = mc_node(i, j);
                surface_band_bending(g_mesh, node, g_config->surface_bb_delV, g_config->surface_bb_direction);
            }
        }

        electric_field(g_mesh);
    }
    else if(g_config->constant_efield_flag == ON) {
        constant_efield(g_mesh, g_mesh->edges[direction_t.TOP][1].potential);
    }

    if(g_config->poisson_flag == ON) {
        // Boundary conditions for the model simulated
        // ===========================================
        if(poisson_boundary_conditions(g_mesh) != 0) {
            printf("Error: Unknown error calculating Poisson boundary conditions.\n");
            return 1;
        }
        if(g_config->faraday_flag) {
            faraday_boundary_conditions(g_mesh);
        }
        printf("Boundary conditions calculated...\n");
    }

    // Initialization for Monte Carlo
    // ==============================
    double transistion_rate[NOAMTIA][DIME][3];
    if(g_config->simulation_model == MCE || g_config->simulation_model == MCEH) {
        for(int i = 0; i < NOAMTIA; i++) {
            calculate_scattering_rates(&g_materials[i]);
            calc_absorption_rates(g_materials[i], transistion_rate);
        }
        printf("Scattering rates calculated...\n");

        if(g_config->photoexcitation_flag == ON) {
            int num = photoexcite_carriers(g_mesh, g_config->photon_energy, transistion_rate, GM);
            printf("Photoexcited %d carriers\n", num);
        }
        else {
            populate_superparticles(g_mesh, 0.8, GM);
        }
        printf("Device configuration complete...\n");
    }
    printf("\n");

    int before = g_config->num_particles;


    if(g_config->photoexcitation_flag == ON) {
        FILE *excited_fp = fopen("photoexcited_particles.csv", "w");
        fprintf(excited_fp, "id x y energy\n");
        for(int n = 0; n < g_config->num_particles; ++n) {
            Particle *p = &g_mesh->particles[n];
            fprintf(excited_fp, "%lld %g %g %g\n", p->id, p->x, p->y, mc_particle_energy(p));
        }
        fclose(excited_fp);
    }

    emitted_fp = fopen("emitted.csv", "w");
    fprintf(emitted_fp, "id time energy\n");

    if(g_config->tracking_output == ON) {
      tracking_fp = fopen("tracking.csv", "w");
      fprintf(tracking_fp, "id time x y E valley\n");
    }

    FILE *particles_fp = fopen("particles.csv", "w");
    fprintf(particles_fp, "timestep time count\n");

    valley_occupation_fp = fopen("valley_occupation.csv", "w");
    fprintf(valley_occupation_fp, "timestep time c1 c2 c3\n");

    velocity_fp = fopen("velocity.csv", "w");
    fprintf(velocity_fp, "timestep x y\n");

    // HERE IS THE SIMULATION
    // ======================
    int valley_occupation[10];
    for(int it = 1; it <= ITMAX; it++) {
        memset(&valley_occupation, 0, sizeof(valley_occupation));
        for(int n = 0; n < g_config->num_particles; ++n) {
            valley_occupation[g_mesh->particles[n].valley] += 1;
        }
        fprintf(valley_occupation_fp, "%d %g %d %d %d\n", it, g_config->time, valley_occupation[1], valley_occupation[2], valley_occupation[3]);

        fprintf(particles_fp, "%d %g %lld\n", it, g_config->time, g_config->num_particles);
        if(it % 10 == 0) {
            fflush(valley_occupation_fp);
        }

        if(updating(it, g_config->simulation_model) != 0) {
            break;
        }
    }

    fclose(particles_fp);
    fclose(emitted_fp);
    fclose(valley_occupation_fp);
    fclose(velocity_fp);
    if(g_config->tracking_output == ON) {
        fclose(tracking_fp);
    }

    // Here we save the outputs
    // ========================
    SaveOutputFiles(g_config->output_format, 0);
    printf("\nFinal Output has been saved\n");

    int after = g_config->num_particles;
    if(g_config->photoexcitation_flag == ON) {
        printf("%d -> %d\n", before, after);
        printf("%d (%lf%%) Particles emitted\n", before - after, 100.0 * (double)(before - after) / (double)before);
    }

    binarytime=time(NULL);
    nowtm=localtime(&binarytime);
    printf("Computation Finished at %s\n",asctime(nowtm));
    return(EXIT_SUCCESS);
}
