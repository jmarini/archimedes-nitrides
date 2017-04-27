/* optical_absorption.h -- This file is part of Archimedes release 1.2.0.
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

#ifndef ARCHIMEDES_ABSORPTION_H
#define ARCHIMEDES_ABSORPTION_H


#include "material.h"
#include "mesh.h"
#include "particle.h"


// Returns the optical joint density of states for the given conduction and valence bands
// evaluated at the given energy
double optical_joint_DOS(Material material, int conduction_band, int valence_band, double energy);


// Returns the transition rate (1/s) between the given conduction and valence bands at the given
// energy
double optical_transition_rate(Material material, int conduction_band,
                               int valence_band, double photon_energy);


// Returns the absorption coefficient at the given energy. Takes into account direct interband
// transitions only from HH, LH, SO valence bands to G-1 conduction band
double absorption_coefficient(Material material, double photon_energy);


// Calculates the number of photoexcited electrons in the given node at the specified energy.
// Photons are assumed to enter device at x=0
// TODO: number of carriers needs to depend on energy - relative W @ E vs max W
int electrons_in_cell(Mesh *mesh, Node *node, double photon_energy);


// Calculated relative absorption rates for the different possible transitions at the given
// energy.
int calc_absorption_rates(Material material, double transistion_rate[NOAMTIA][DIME][3]);


Particle create_photoexcited_carrier(Node *node, double photon_energy,
                                     double total_scattering_rate[NOAMTIA+1],
                                     int conduction_band, int valence_band);


int photoexcite_carriers(Mesh *mesh, double photon_energy,
                         double transistion_rate[NOAMTIA][DIME][3],
                         double total_scattering_rate[NOAMTIA+1]);


#endif
