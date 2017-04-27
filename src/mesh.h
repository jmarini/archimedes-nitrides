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


#include "material.h"
#include "particle.h"
#include "vec.h"


// maximum number of mesh cells in x and y
#define NXM 308
#define NYM 308


typedef struct {
    double density;
    Vec2   velocity; // running sum velocity
    double energy;   // running sum energy
} Carrier_Info;


typedef struct {
    union {     // index of the node, provides i & j
        struct Index;
        Index index;
    };

    Material *material;
    int material_id;

    double donor_conc;
    double acceptor_conc;
    double fixed_charge;

    Carrier_Info e;     // electrons
    Carrier_Info h;     // holes

    double qep;         // quantum effective potential
    double potential;
    Vec2 efield;         // electric field
    double magnetic_field;
} Node;


typedef struct {
    int boundary;
    double potential;
    double n;
    double p;
} Edge;


typedef struct {
    int nx; // number of cells in x-direction
    int ny; //                    y-direction

    double dx; // cell size in x-direction
    double dy; //              y-direction

    union { // size of the mesh, provides width & height
        struct Dimensions;
        Dimensions size;
    };

    int num_nodes;
    int num_triangles;

    Node nodes[NXM + 1][NYM + 1];
    Edge edges[4][NXM + 1]; // edges, indexed by direction and index (i or j)

    Vec2 coordinates[NXM * NYM];
    int triangles[NXM * NYM][3];

    Particle particles[NPMAX + 1];
} Mesh;


typedef struct {
    int BOTTOM;
    int RIGHT;
    int TOP;
    int LEFT;
} Direction;


typedef struct {
    int INSULATOR;
    int SCHOTTKY;
    int OHMIC;
    int VACUUM;
} Boundary;


int mc_build_mesh(Mesh *mesh);
int mc_save_mesh(Mesh *mesh, char *filename);


int mc_is_boundary_insulator(int direction, int index);
int mc_is_boundary_schottky(int direction, int index);
int mc_is_boundary_ohmic(int direction, int index);
int mc_is_boundary_vacuum(int direction, int index);
int mc_is_boundary_contact(int direction, int index);


Node * mc_node(int i, int j);
Node * mc_node_s(Index index);


Node * mc_get_particle_node(Particle *p);
Vec2 mc_random_location_in_node(Node *node);


// define global extern variable
extern Mesh *g_mesh;
extern Direction direction_t;
extern Boundary boundary_t;


#endif
