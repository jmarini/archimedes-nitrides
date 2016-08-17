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
// Created on 20 oct.2015, J. Marini
// Last modif. : 20 oct.2015, J. Marini
// ######################################################

#ifndef ARCHIMEDES_MESH_H
#define ARCHIMEDES_MESH_H

#define NXM 308
#define NXY 308


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
    int nx; // number of cells in x-direction
    int ny; //                    y-direction

    real dx; // cell size in x-direction
    real dy; //              y-direction

    real width;  // device width
    real height; //        height

    int num_nodes;
    int num_triangles;

    mc_node_information_t info[NXM + 1][NYM + 1];

    real nodes[NXM * NYM][2];
    int triangles[NXM * NYM][3];
} mc_mesh_t;


extern mc_mesh_t *g_mesh;


int mc_build_mesh(mc_mesh_t *mesh);

int mc_save_mesh(mc_mesh_t *mesh, char *filename);


#endif
