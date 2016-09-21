/* configuration.h -- This file is part of Archimedes release 1.2.0.
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

#ifndef ARCHIMEDES_CONFIGURATION_H
#define ARCHIMEDES_CONFIGURATION_H


typedef struct {
    int simulation_model;

    // particle info
    long long int num_particles;
    long long int next_particle_id;
    double carriers_per_superparticle;

    // scattering mechanism flags
    int optical_phonon_scattering;
    int acoustic_phonon_scattering;
    int impurity_scattering;
    int piezoelectric_scattering;

    // band structure & quantum correction models
    int conduction_band;
    int quantum_flag;
    int qep_model;
    double qep_alpha;
    double qep_gamma;

    int photoexcitation_flag;
    double photon_energy;

    int faraday_flag;
    double screening_length;

    double lattice_temp;
    double impurity_conc;

    // averaging
    int particles_per_cell;
    int avg_steps; // media
    double avg_alpha;

    // intput / output
    int save_mesh;
    int max_min_output;
    int save_step_output;
    int scattering_output;
    int output_format;
    int load_initial_data;
    int tcad_data;

    // simulation timing parameters
    double time;
    double tf;
    double dt;
    double tauw; // MEP

    double max_doping;
} Configuration;


void mc_initialize_config(Configuration *config);


extern Configuration *g_config;


#endif
