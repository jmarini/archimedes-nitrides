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


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 13 Sep.2007, Siracusa, J.M.Sellier
// ######################################################

// Initial Configuration of the Particles (electrons)
// simulated in the devices.

void MCdevice_config(void) {
    long int n = 0;
    int i  = 0,
        j  = 0,
        np = 0,
        m  = 0;
    int valley   = 0,
        material = 0;
    real c1 = 0.,
         c2 = 0.,
         c3 = 0.,
         c4 = 0.,
         c5 = 0.,
         c6 = 0.,
         c7 = 0.;
    int nx  = g_mesh->nx,
        ny  = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;

    // Number of carriers per particle
    g_config->carriers_per_particle = g_config->max_doping * dx * dy / g_config->particles_per_cell;

    for(i=1;i<=nx+1;i++) {
        for(j=1;j<=ny+1;j++){
            material = g_mesh->info[i][j].material;
            // np=number of particles in the (i,j)-th cell
            if(g_config->load_initial_data == 1) {
                np = (int)(u2d[i][j][1] * dx * dy / g_config->carriers_per_particle + 0.5);
            }
            else {
                np = (int)(g_mesh->info[i][j].donor_conc * dx * dy / g_config->carriers_per_particle + 0.5);
            }
            if((i == 1) || (i == nx + 1)) { np /= 2; }
            if((j == 1) || (j == ny + 1)) { np /= 2; }
            if(np > 0) {
                for(m = 1; m <= np; m++) {
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
                    if(g_config->load_initial_data == 0) {
                        valley=1;
                        c1=log(rnd());
                        if(NOVALLEY[material]==1) {
                            c2=SMH[material][0]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[material][1]*1.5*BKTQ*c1));
                        }
                        if(NOVALLEY[material]>=2){
                            // 80% of the created particles goes in the first valley
                            valley=1;
                            c2=SMH[material][valley]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[material][1]*1.5*BKTQ*c1));
                            // 20% of the created particles goes in the second valley
                            if(rnd()>0.8){
                                valley=2;
                                c2=SMH[material][valley]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[material][2]*1.5*BKTQ*c1));
                            }
                        }
                    }
                    // The following is in case of precendtly loaded initial data.
                    else if(g_config->load_initial_data == 1) {
                        valley=1;
                        // In this case c1 represents the mean electron energy
                        // loaded from precedent simulations and have nothing to
                        // do with the lattice energy.
                        c1=-u2d[i][j][4]/u2d[i][j][1]/Q;
                        if(NOVALLEY[material]==1) {
                            c2=SMH[material][0]*sqrt(-1.5*c1*(1.-alphaK[material][1]*1.5*c1));
                        }
                        if(NOVALLEY[material]>=2){
                            valley=1;
                            c2=SMH[material][valley]*sqrt(-1.5*c1*(1.-alphaK[material][1]*1.5*c1));
                            if(rnd()>0.8){
                                valley=2;
                                c2=SMH[material][valley]*sqrt(-1.5*c1*(1.-alphaK[material][2]*1.5*c1));
                            }
                        }
                    }
                    c3=1.-2.*rnd();
                    c4=sqrt(1.-c3*c3);
                    c5=2.*PI*rnd();
                    c6=sin(c5);
                    c7=cos(c5);
                    P[n].id = g_config->next_particle_id++;
                    P[n].valley = valley;
                    P[n].kx = c2 * c3 * c6;
                    P[n].ky = c2 * c4 * c6;
                    P[n].kz = c2 * c7;
                    P[n].t  = -log(rnd())/GM[material];
                    P[n].x  = dx*(rnd()+(real)(i)-1.5);
                    P[n].y  = dy*(rnd()+(real)(j)-1.5);
                    if(i==1) P[n].x=dx*0.5*rnd();
                    if(j==1) P[n].y=dy*0.5*rnd();
                    if(i==nx+1) P[n].x=g_mesh->width-dx*0.5*rnd();
                    if(j==ny+1) P[n].y=g_mesh->height-dy*0.5*rnd();
                }
            }
        }
    }

    g_config->num_particles = n;
    for(i=1;i<=nx+1;i++) {
        for(j=1;j<=ny+1;j++) {
            u2d[i][j][2]=0.;
            u2d[i][j][3]=0.;
            u2d[i][j][4]=0.;
        }
    }

    printf("Initial Number of Electron Super-particles = %lld\n", g_config->num_particles);
}
