/* emc_threaded.h -- This file is part of Archimedes release 1.2.0.
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
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


// ######################################################
// Created on 30 oct.2015, J. Marini
// Last modif. : 30 oct.2015, J. Marini
// ######################################################


typedef struct {
    int start;
    int tid;
} thread_data_t;


thread_data_t thread_data[NUMBER_THREADS];


void EMC_worker(thread_data_t *data)
{
    int start = data->start;
    int tid = data->tid;
    int stop = start + SEGMENT_SIZE;
    real ti  = TEMPO,
         tdt = TEMPO + DT,
         tau = 0.0;
    int i = 0,
        j = 0,
        n = 0;

    for(n = start; n < stop && n < INUM; n++) {
        particle_t *particle = &P[n];

        // while the particle's time is less than the time for the step...
        while(particle->t <= tdt) {
            tau = particle->t - ti;                // the dt for the current step
            drift(particle, tau);                  // drift for dt
            mc_particle_coords(particle, &i, &j);  // get updated particle coords
            scatter(particle, i_dom[i][j], tid);   // scatter particle
            ti = particle->t;                      // update the time
            particle->t += -log(trnd(tid)) / GM[i_dom[i][j]]; // update particle time
        }
        tau = tdt - ti;              // calculate unused time in step
        drift(particle, tau);        // drift for unused time in step
    }


}



void EMC_cleanup(void)
{
    int i = 0,
        j = 0,
        n = 1,
        ni;
    int npt[NXM+NYM+1][NUMBER_DIRECTIONS];
    memset(&npt, 0, sizeof(npt));

    for(n = 1; n < INUM; n++) {
        particle_t *particle = &P[n];


        i = (int)(particle->x / dx + 1.5);
        j = (int)(particle->y / dy + 1.5);
        mc_check_particle_leaving(particle, i, j, direction_t.BOTTOM, npt);
        mc_check_particle_leaving(particle, i, j, direction_t.RIGHT, npt);
        mc_check_particle_leaving(particle, i, j, direction_t.TOP, npt);
        mc_check_particle_leaving(particle, i, j, direction_t.LEFT, npt);

    }


    for(n = INUM; n > 0; n--) {
        if(!mc_does_particle_exist(&P[n])) {
            P[n] = P[INUM];
            INUM--;
        }
    }


    // create particles at ohmic contacts of the bottom edge
    for(i=1; i<=nx+1; i++) {
        if(mc_boundary_type(direction_t.BOTTOM, i) == boundary_t.OHMIC) {
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
        if(mc_boundary_type(direction_t.TOP, i) == boundary_t.OHMIC) {
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
        if(mc_boundary_type(direction_t.RIGHT, i) == boundary_t.OHMIC) {
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
        if(mc_boundary_type(direction_t.LEFT, i) == boundary_t.OHMIC) {
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

}


void EMC_threaded(void)
{
    int n = 1,  // particle index
        t = 0;  // thread index

    printf("\n"
           "  NUMBER PARTICLES = %'d\n"
           "  NUMBER THREADS   = %'d\n"
           "  SEGMENT SIZE     = %'d\n"
           "  TOTAL QUEUED     = %.2f\n",
           INUM, NUMBER_THREADS, SEGMENT_SIZE, (float)INUM / (float)SEGMENT_SIZE);

    while(n < INUM) {
        thread_data[t] = (thread_data_t){.start=n, .tid=t};

        thpool_add_work(thread_pool, (void *)EMC_worker, &thread_data[t]);

        n += SEGMENT_SIZE;
        t++;
    }

    thpool_wait(thread_pool);

    EMC_cleanup();
}
