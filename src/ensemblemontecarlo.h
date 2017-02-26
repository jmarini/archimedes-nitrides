/* ensemblemontecarlo.h -- This file is part of Archimedes release 1.2.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It includes some quantum effects by means
   of effective potential method. It is also able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier
   <jeanmichel.sellier@gmail.com>
   <jsellier@purdue.edu>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANparticle->y WARRANTparticle->y; without even the implied warranty of
   MERCHANTABILITparticle->y or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   particle->you should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "particle.h"
#include "mesh.h"


// Ensemble Monte Carlo method
void EMC(int iteration) {
    int i = 0,
        j = 0;
    real tdt = g_config->time + g_config->dt,
         tau = 0.;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;
    Node *node = NULL;

    int npt[NXM+NYM+1][4];
    memset(&npt, 0, sizeof(npt));

    long int n = 1;  // index of current particle
    do {
        Particle *particle = &P[n];
        // information about particle n is set up in easy access variables
        real ti = g_config->time;

        // while the particle's time is less than the time for the step...
        while(particle->t <= tdt) {
            tau = particle->t - ti;                // the dt for the current step
            drift(particle, tau);                  // drift for dt
            node = mc_get_particle_node(particle);

            if(g_config->tracking_output == ON
               && particle->id % g_config->tracking_mod == 0) {
                mc_print_tracking(iteration, particle);
            }
            int s = scatter(particle, node->mat);        // scatter particle


            if(g_config->tracking_output == ON
               && particle->id % g_config->tracking_mod == 0
               && s) {
                mc_print_tracking(iteration, particle);
            }

            ti = particle->t;                      // update the time
            particle->t = ti - log(rnd()) / GM[node->material]; // update particle time
        }
        tau = tdt - ti;              // calculate unused time in step
        drift(particle, tau);        // drift for unused time in step


        // check if a particle is going out from the right edge of the device
        int direction = direction_t.RIGHT;
        if(mc_does_particle_exist(particle)) {
            Index index = mc_particle_edge_coords(particle);
            int i = index.i,
                j = index.j;
            if(i >= nx + 1 && mc_is_boundary_contact(direction, j)) {
                mc_remove_particle(particle);
                if(npt[j][direction] < (g_config->particles_per_cell/2) && j > 1 && j < ny+1){
                    npt[j][direction]++;
                    particle->valley = 1;
                }
                else if(npt[j][direction] < (g_config->particles_per_cell/4) && (j <= 1 || j >= ny+1)){
                    npt[j][direction]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the left edge of the device
        direction = direction_t.LEFT;
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(i<=1 && mc_is_boundary_contact(direction, j)) {
                mc_remove_particle(particle);
                if(npt[j][direction]<(g_config->particles_per_cell/2) && j>1 && j<ny+1){
                    npt[j][direction]++;
                    particle->valley = 1;
                }
                else if(npt[j][direction]<(g_config->particles_per_cell/4) && (j<=1 || j>=ny+1)){
                    npt[j][direction]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the bottom edge of the device
        direction = direction_t.BOTTOM;
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(j<=1 && mc_is_boundary_contact(direction, i)) {
                mc_remove_particle(particle);
                if(npt[i][direction]<(g_config->particles_per_cell/2) && (i>1 || i<nx+1)){
                    npt[i][direction]++;
                    particle->valley = 1;
                }
                if(npt[i][direction]<(g_config->particles_per_cell/4) && (i<=1 || i>=nx+1)){
                    npt[i][direction]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the upper edge of the device
        direction = direction_t.TOP;
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(j>=ny+1 && mc_is_boundary_contact(direction, i)) {
                mc_remove_particle(particle);
                if(npt[i][direction]<(g_config->particles_per_cell/2) && (i>1 || i<nx+1)){
                    npt[i][direction]++;
                    particle->valley = 1;
                }
                if(npt[i][direction]<(g_config->particles_per_cell/4) && (i<=1 || i>=nx+1)){
                    npt[i][direction]++;
                    particle->valley = 1;
                }
            }
        }

        if(mc_does_particle_exist(particle)) { ++n; }
        else {
            // move the last particle to the now empty spot to keep
            // a dense array in constant time. Entries beyond num_particles
            // can be assumed to be "empty"
            P[n] = P[g_config->num_particles - 1];
            --g_config->num_particles;
        }
    } while(n < g_config->num_particles);

    // if(iteration <= 100) {
    //     char s[150];
    //     sprintf(s, "tracking%03d.bin", iteration);
    //     FILE *fp = fopen(s, "wb");
    //     fwrite(PTRACK, sizeof(particle_tracking_t) * n, 1, fp);
    //     fclose(fp);
    // }
    // else {
    //     exit(0);
    // }


    // create particles at ohmic contacts of the bottom edge
    for(i=1; i<=nx+1; i++) {
        int direction = direction_t.BOTTOM;
        if(mc_is_boundary_ohmic(direction, i)) {
            int ni=(g_config->particles_per_cell/2)-npt[i][direction];
            if(i==1 || i==nx+1) {
                ni=g_config->particles_per_cell/4-npt[i][direction];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=g_config->num_particles+j;
                    P[n] = creation(i, g_config->time, direction);
                }
            g_config->num_particles += ni;
            }
        }
    }

    // create particles at ohmic contacts of the upper edge
    for(i=1; i<=nx+1; i++) {
        int direction = direction_t.TOP;
        if(mc_is_boundary_ohmic(direction, i)) {
            int ni=(g_config->particles_per_cell/2)-npt[i][direction];
            if(i==1 || i==nx+1) {
                ni=g_config->particles_per_cell/4-npt[i][direction];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=g_config->num_particles+j;
                    P[n] = creation(i, g_config->time, direction);
                }
            g_config->num_particles += ni;
            }
        }
    }

    // create particles at ohmic contacts of the right edge
    for(i=1; i<=ny+1; i++) {
        int direction = direction_t.RIGHT;
        if(mc_is_boundary_ohmic(direction, i)) {
            int ni=(g_config->particles_per_cell/2)-npt[i][direction];
            if(i==1 || i==ny+1) {
                ni=g_config->particles_per_cell/4-npt[i][direction];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=g_config->num_particles+j;
                    P[n] = creation(i, g_config->time, direction);
                }
            g_config->num_particles += ni;
            }
        }
    }

    // create particles at ohmic contacts of the left edge
    for(i=1; i<=ny+1; i++) {
        int direction = direction_t.LEFT;
        if(mc_is_boundary_ohmic(direction, i)) {
            int ni=(g_config->particles_per_cell/2)-npt[i][direction];
            if(i==1 || i==ny+1) {
                ni=g_config->particles_per_cell/4-npt[i][direction];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++){
                    n=g_config->num_particles+j;
                    P[n] = creation(i, g_config->time, direction);
                }
                g_config->num_particles += ni;
            }
        }
    }

    printf("\nActual number of electron super-particles = %lld\n", g_config->num_particles);
    if(g_config->num_particles > NPMAX) {
        printf("%s: too big actual number of particles\n", progname);
        exit(EXIT_FAILURE);
    }
}

// ============================================================
