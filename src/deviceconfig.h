/* deviceconfig.h -- This file is part of Archimedes release 0.1.2.
Archimedes is a simulator for Submicron 2D III-V semiconductor
Devices. It implements both the Monte Carlo method
for the simulation of the semiclassical Boltzmann equation for both
electrons and holes. It also includes the quantum effects by means 
of effective potential method. It is now able to simulate applied
magnetic fields along with self consistent Faraday equation.

Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>

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

#ifndef ARCHIMEDES_DEVICE_CONFIG_H
#define ARCHIMEDES_DEVICE_CONFIG_H


#include "mesh.h"
#include "particle_creation.h"
#include "vec.h"


// Initial Configuration of the Particles (electrons)
// simulated in the devices.
void MCdevice_config(Mesh *mesh) {
    long int n = 0;
    int nx  = mesh->nx,
        ny  = mesh->ny;
    real dx = mesh->dx,
         dy = mesh->dy;

    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 1; j <= ny + 1; ++j){
            Node *node = mc_node(i, j);

            int np = 0; // number of superparticles in the (i,j)-th cell
            if(g_config->load_initial_data == ON) {
                np = (int)(u2d[i][j][1] * dx * dy / g_config->carriers_per_superparticle + 0.5);
            }
            else if(g_config->tcad_data == ON) {
                np = (int)(node->e.density * dx * dy / g_config->carriers_per_superparticle + 0.5);
            }
            else {
                np = (int)(node->donor_conc * dx * dy / g_config->carriers_per_superparticle + 0.5);
            }
            if((i == 1) || (i == nx + 1)) { np /= 2; }
            if((j == 1) || (j == ny + 1)) { np /= 2; }
            if(np > 0) {
                for(int m = 1; m <= np; m++) {
                    n++;
                    if(n > NPMAX) {
                        printf("%s: too big number of particles\n", progname);
                        exit(EXIT_FAILURE);
                    }

                    mesh->particles[n] = create_particle(mesh, node, 0.8, GM);
                }
            }
        }
    }

    g_config->num_particles = n;

    // compatability
    for(int i=1;i<=nx+1;i++) {
        for(int j=1;j<=ny+1;j++) {
            u2d[i][j][2]=0.;
            u2d[i][j][3]=0.;
            u2d[i][j][4]=0.;
        }
    }

    printf("Initial Number of Electron Super-particles = %lld\n", g_config->num_particles);
}

#endif
