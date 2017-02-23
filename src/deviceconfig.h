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
            int material = node->material;
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

                    // We assume that the particles are initially
                    // at near thermal equilibrium

                    // In the case of two-valleys materials, 80% of the electrons are considered in the Gamma
                    // valley in the starting simulation time, while the other 20% are
                    // in the L-valley.
                    int valley = 1;
                    double kf = 0.;

                    // The following is in case of precendtly loaded initial data.
                    if(g_config->load_initial_data == ON) {
                        valley = 1;
                        // In this case c1 represents the mean electron energy
                        // loaded from precedent simulations and have nothing to
                        // do with the lattice energy.
                        double c1=-u2d[i][j][4]/u2d[i][j][1]/Q;
                        if(g_materials[material].cb.num_valleys == 1) {
                            kf=g_materials[material].cb.smh[0]*sqrt(-1.5*c1*(1.-g_materials[material].cb.alpha[1]*1.5*c1));
                        }
                        else if(g_materials[material].cb.num_valleys >= 2) {
                            valley=1;
                            kf=g_materials[material].cb.smh[valley]*sqrt(-1.5*c1*(1.-g_materials[material].cb.alpha[1]*1.5*c1));
                            if(rnd()>0.8){
                                valley=2;
                                kf=g_materials[material].cb.smh[valley]*sqrt(-1.5*c1*(1.-g_materials[material].cb.alpha[2]*1.5*c1));
                            }
                        }
                    }
                    else {
                        double tmp = 1.5 * BKTQ * log(rnd());;
                        valley = 1;
                        if(node->mat->cb.num_valleys >= 2 && rnd() > 0.8) {
                            valley = 2; // 80% of particles in 1st valley, 20% in 2nd valley
                        }
                        kf = node->mat->cb.smh[valley]
                           * sqrt(-tmp * (1. - node->mat->cb.alpha[valley] * tmp));
                    }

                    double c3 = 1. - 2. * rnd();
                    double c4 = sqrt(1. - c3 * c3);
                    double c5 = 2. * PI * rnd();
                    double c6 = sin(c5);
                    double c7 = cos(c5);
                    P[n].id = mc_next_particle_id( );
                    P[n].valley = valley;
                    P[n].kx = kf * c3 * c6;
                    P[n].ky = kf * c4 * c6;
                    P[n].kz = kf * c7;
                    P[n].t  = -log(rnd()) / GM[material];
                    P[n].x  = dx * (rnd() + (double)(i) - 1.5);
                    P[n].y  = dy * (rnd() + (double)(j) - 1.5);
                    P[n].photoemission_flag = 0;

                    if(i == 1) { P[n].x = dx * 0.5 * rnd(); }
                    if(j == 1) { P[n].y = dy * 0.5 * rnd(); }
                    if(i == nx + 1) { P[n].x = mesh->width  - dx * 0.5 * rnd(); }
                    if(j == ny + 1) { P[n].y = mesh->height - dy * 0.5 * rnd(); }
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
