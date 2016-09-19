/* particle.h -- This file is part of Archimedes release 1.2.0.
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

#ifndef ARCHIMEDES_PARTICLE_H
#define ARCHIMEDES_PARTICLE_H


#include <math.h>

#include "global_defines.h"
#include "mesh.h"
#include "vec.h"


typedef struct {
    long long int id;      // unique identifier used to track particle
    int valley;  // number id of the valley the particle is in
    real kx;     // momentum coordinates - relative to valley minimum
    real ky;
    real kz;
    real t;      // time
    real x;      // real coordinates
    real y;

    int photoemission_flag; // flag tracking state of electron through photoemission
                            //   0: normal electron (default)
                            //   1: photoexcited carrier
                            //   2: photoemitted carrier
} Particle;


typedef struct {
    long long int id;         // unique identifier used to track particle
    int valley;     // valley the particle is in
    real kx;        // momentum coordinates - relative to valley minimum
    real ky;
    real kz;
    real energy;    // particle energy - relative to valley minimum
    real t;         // time
    real x;         // real coordinates
    real y;
    int i;          // mesh cell
    int j;
    real vx;        // real-space velocity
    real vy;

    int photoemission_flag; // flag tracking state of electron through photoemission
                            //   0: normal electron (default)
                            //   1: photoexcited carrier
                            //   2: photoemitted carrier
} particle_info_t;


inline int mc_does_particle_exist(Particle *particle) {
    return particle->valley != 9;
}


inline void mc_remove_particle(Particle *particle) {
    particle->valley = 9;
}


inline real mc_particle_ksquared(Particle *particle) {
    return (particle->kx * particle->kx +
            particle->ky * particle->ky +
            particle->kz * particle->kz);
}


inline real mc_particle_k(Particle *particle) {
    return sqrt(mc_particle_ksquared(particle));
}


Index mc_particle_coords(Particle *particle);
Index mc_particle_edge_coords(Particle *particle);

Node * mc_get_particle_node(Particle *particle);

long long int mc_next_particle_id( );

#endif
