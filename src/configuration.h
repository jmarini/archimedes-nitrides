/* mesh.h -- This file is part of Archimedes release 1.2.0.
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
// Created on 12 aug.2016, J. Marini
// Last modif. : 12 aug.2016, J. Marini
// ######################################################

#ifndef ARCHIMEDES_CONFIGURATION_H
#define ARCHIMEDES_CONFIGURATION_H


typedef struct {
    real qep;             // quantum effective potential
    real potential;
    real efield_x;        // electric field - x
    real efield_y;        //                - y
    real magnetic_field;
} mc_poisson_t;


typedef struct {
    real density;
    real xvel;    // running sum velocity - x
    real yvel;    //                      - y
    real energy;  // running sum energy
} mc_carrier_t;


typedef struct {
    int material;
    real donor_conc;
    real acceptor_conc;

    mc_carrier_t e;     // electrons
    mc_carrier_t h;     // holes

    mc_poisson_t poisson;
} mc_node_information_t;


typedef struct {
    int simulation_model;

    // scattering mechanism flags
    int optical_phonon_scattering;
    int acoustic_phonon_scattering;
    int impurity_scattering;
    int piezoelectric_scattering;

    // band structure & quantum correction models
    int conduction_band;
    int quantum_flag;
    int qep_model;
    real qep_alpha;
    real qep_gamma;

    int faraday_flag;

    real lattice_temp;
    real impurity_conc;

    // averaging
    int particles_per_cell;
    int avg_steps; // media
    real avg_alpha;

    // intput / output
    int save_mesh;
    int max_min_output;
    int save_step_output;
    int scattering_output;
    int output_format;
    int load_initial_data;

    // simulation timing parameters
    real tf;
    real dt;
    real tauw;

    // mesh & geometry parameters
    int nx;
    int ny;
    real width;  // x
    real height; // y
    real dx;
    real dy;

    real max_doping;

    mc_node_information_t nodes[NXM + 1][NYM + 1];
} mc_configuration_t;


extern mc_configuration_t *mc_config;


#endif
