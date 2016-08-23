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


#ifndef ARCHIMEDES_MESH_H
#define ARCHIMEDES_MESH_H


#include "global_defines.h"
#include "material.h"
#include "vec.h"


// maximum number of mesh cells in x and y
#define NXM 308
#define NXY 308


typedef struct {
    real density;
    vec2d velocity; // running sum velocity
    real energy;  // running sum energy
} mc_carrier_t;


typedef struct {
    union {     // index of the node, provides i & j
        struct  index_s;
        index_s index;
    };

    int material;
    real donor_conc;
    real acceptor_conc;

    mc_carrier_t e;     // electrons
    mc_carrier_t h;     // holes

    real qep;             // quantum effective potential
    real potential;
    vec2d efield;         // electric field
    real magnetic_field;
} mc_node_t;


typedef struct {
    int nx; // number of cells in x-direction
    int ny; //                    y-direction

    real dx; // cell size in x-direction
    real dy; //              y-direction

    union { // size of the mesh, provides width & height
        struct  size_s;
        size_s size;
    };

    int num_nodes;
    int num_triangles;

    mc_node_t info[NXM + 1][NYM + 1];

    real coordinates[NXM * NYM][2];
    int triangles[NXM * NYM][3];
} mc_mesh_t;


int mc_build_mesh(mc_mesh_t *mesh);


int mc_save_mesh(mc_mesh_t *mesh, char *filename);


mc_node_t * mc_node(int i, int j);
mc_node_t * mc_node_s(index_s index);


// define global extern variable
extern mc_mesh_t *g_mesh;


#endif
