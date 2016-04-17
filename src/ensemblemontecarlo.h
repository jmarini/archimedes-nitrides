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


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 25 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Ensemble Monte Carlo method

void EMC(void)
{
    long int n=1;  // index of current particle
    int i, ni, j, npt[NXM+NYM+1][4];
    real tdt, ti, tau;

    memset(&npt, 0, sizeof(npt));
    tdt = TEMPO + DT;

    do {
        particle_t *particle = &P[n];
        // information about particle n is set up in easy access variables
        ti = TEMPO;

        // while the particle's time is less than the time for the step...
        while(particle->t <= tdt) {
            tau = particle->t - ti;                // the dt for the current step
            drift(particle, tau);                  // drift for dt
            mc_particle_coords(particle, &i, &j);  // get updated particle coords
            scatter(particle, i_dom[i][j]);        // scatter particle
            ti = particle->t;                      // update the time
            mc_particle_coords(particle, &i, &j);  // get updated particle coords
            particle->t = ti - log(rnd()) / GM[i_dom[i][j]]; // update particle time
        }
        tau = tdt - ti;              // calculate unused time in step
        drift(particle, tau);        // drift for unused time in step

        // check if a particle is going out from the right edge of the device
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(i >= nx + 1 && mc_is_boundary_contact(direction_t.RIGHT, j)) {
                mc_remove_particle(particle);
                if(npt[j][1] < (NP1/2) && j > 1 && j < ny+1){
                    npt[j][1]++;
                    particle->valley = 1;
                }
                else if(npt[j][1] < (NP1/4) && (j <= 1 || j >= ny+1)){
                    npt[j][1]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the left edge of the device
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(i<=1 && mc_is_boundary_contact(direction_t.LEFT, j)) {
                mc_remove_particle(particle);
                if(npt[j][3]<(NP1/2) && j>1 && j<ny+1){
                    npt[j][3]++;
                    particle->valley = 1;
                }
                else if(npt[j][3]<(NP1/4) && (j<=1 || j>=ny+1)){
                    npt[j][3]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the bottom edge of the device
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(j<=1 && mc_is_boundary_contact(direction_t.BOTTOM, i)) {
                mc_remove_particle(particle);
                if(npt[i][0]<(NP1/2) && (i>1 || i<nx+1)){
                    npt[i][0]++;
                    particle->valley = 1;
                }
                if(npt[i][0]<(NP1/4) && (i<=1 || i>=nx+1)){
                    npt[i][0]++;
                    particle->valley = 1;
                }
            }
        }

        // check if a particle is going out from the upper edge of the device
        if(mc_does_particle_exist(particle)) {
            i = (int)(particle->x / dx + 1.5);
            j = (int)(particle->y / dy + 1.5);
            if(j>=ny+1 && mc_is_boundary_contact(direction_t.TOP, i)) {
                mc_remove_particle(particle);
                if(npt[i][2]<(NP1/2) && (i>1 || i<nx+1)){
                    npt[i][2]++;
                    particle->valley = 1;
                }
                if(npt[i][2]<(NP1/4) && (i<=1 || i>=nx+1)){
                    npt[i][2]++;
                    particle->valley = 1;
                }
            }
        }

        // ==============================================
        if(mc_does_particle_exist(particle)) {
            n++;
        }
        else {
            P[n] = P[INUM];
            INUM--;
        }
    } while(n < INUM);

    // create particles at ohmic contacts of the bottom edge
    for(i=1; i<=nx+1; i++) {
        if(mc_is_boundary_ohmic(direction_t.BOTTOM, i)) {
            ni=(NP1/2)-npt[i][0];
            if(i==1 || i==nx+1) {
                ni=NP1/4-npt[i][0];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=INUM+j;
                    P[n] = creation(i,TEMPO,direction_t.BOTTOM);
                }
            INUM += ni;
            }
        }
    }

    // create particles at ohmic contacts of the upper edge
    for(i=1; i<=nx+1; i++) {
        if(mc_is_boundary_ohmic(direction_t.TOP, i)) {
            ni=(NP1/2)-npt[i][2];
            if(i==1 || i==nx+1) {
                ni=NP1/4-npt[i][2];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=INUM+j;
                    P[n] = creation(i,TEMPO,direction_t.TOP);
                }
            INUM += ni;
            }
        }
    }

    // create particles at ohmic contacts of the right edge
    for(i=1; i<=ny+1; i++) {
        if(mc_is_boundary_ohmic(direction_t.RIGHT, i)) {
            ni=(NP1/2)-npt[i][1];
            if(i==1 || i==ny+1) {
                ni=NP1/4-npt[i][1];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++) {
                    n=INUM+j;
                    P[n] = creation(i,TEMPO,direction_t.RIGHT);
                }
            INUM += ni;
            }
        }
    }

    // create particles at ohmic contacts of the left edge
    for(i=1; i<=ny+1; i++) {
        if(mc_is_boundary_ohmic(direction_t.LEFT, i)) {
            ni=(NP1/2)-npt[i][3];
            if(i==1 || i==ny+1) {
                ni=NP1/4-npt[i][3];
            }
            if(ni > 0) {
                for(j=1;j<=ni;j++){
                    n=INUM+j;
                    P[n] = creation(i,TEMPO,direction_t.LEFT);
                }
                INUM += ni;
            }
        }
    }

    printf("\nActual number of electron super-particles = %d\n", INUM);
    if(INUM > NPMAX) {
        printf("%s: too big actual number of particles\n", progname);
        exit(EXIT_FAILURE);
    }
}

// ============================================================
