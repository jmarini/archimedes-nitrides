/* poissonbcs.h -- This file is part of Archimedes release 0.0.5.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and a simplified
   MEP model for the simulation of the semiclassical Boltzmann
   equation for both electrons and holes. It also includes the
   quantum effects by means of effective potential method.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.  */


// Boundary Conditions for the non-stationary Poisson equation
// For more informations about this equation
// see the manual of GNU Archimedes release 0.0.5.


#include "global_defines.h"
#include "mesh.h"
#include "utility.h"


int poisson_boundary_conditions(void) {
    if(SIO2_UP_FLAG || SIO2_DOWN_FLAG) {
        printf("Error: SIO2 flag is deprecated.\n");
        return 1;
    }

    int nx = g_mesh->nx,
        ny = g_mesh->ny;

    // These are completely generic boundary conditions
    for(int i = 1; i <= nx + 1; ++i) {
        // Bottom Edge
        // ===========

        // insulator
        if(mc_is_boundary_insulator(direction_t.BOTTOM, i) || mc_is_boundary_vacuum(direction_t.BOTTOM, i)) {
            // if there is no potential on insulator
            g_mesh->nodes[i][1].potential = g_mesh->nodes[i][2].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(EDGE[direction_t.BOTTOM][i][1] != 0.0) {
                g_mesh->nodes[i][1].potential = EDGE[direction_t.BOTTOM][i][1];
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.BOTTOM, i)) {
            g_mesh->nodes[i][1].potential = EDGE[direction_t.BOTTOM][i][1];
        }

        // Upper Edge
        // ==========

        // insulator
        if(mc_is_boundary_insulator(direction_t.TOP, i) || mc_is_boundary_vacuum(direction_t.TOP, i)) {
            // if there is no potential on insulator
            g_mesh->nodes[i][ny + 1].potential = g_mesh->nodes[i][ny    ].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(EDGE[direction_t.TOP][i][1] != 0.0) {
                g_mesh->nodes[i][ny + 1].potential = EDGE[direction_t.TOP][i][1];
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.TOP, i)) {
            g_mesh->nodes[i][ny + 1].potential = EDGE[direction_t.TOP][i][1];
        }

    }


    for(int j = 1; j <= ny + 1; ++j) {
        // Left Edge
        // =========

        // insulator
        if(mc_is_boundary_insulator(direction_t.LEFT, j) || mc_is_boundary_vacuum(direction_t.LEFT, j)) {
            g_mesh->nodes[1][j].potential = g_mesh->nodes[2][j].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(EDGE[direction_t.LEFT][j][1] != 0.0) {
                g_mesh->nodes[1][j].potential = EDGE[direction_t.LEFT][j][1];
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.LEFT, j)) {
            g_mesh->nodes[1][j].potential = EDGE[direction_t.LEFT][j][1];
        }

        // Right Edge
        // ==========

        // insulator
        if(mc_is_boundary_insulator(direction_t.RIGHT, j) || mc_is_boundary_vacuum(direction_t.RIGHT, j)) {
            g_mesh->nodes[nx + 1][j].potential = g_mesh->nodes[nx - 1][j].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(EDGE[direction_t.RIGHT][j][1] != 0.0) {
                g_mesh->nodes[nx + 1][j].potential = EDGE[direction_t.RIGHT][j][1];
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.RIGHT, j)) {
            g_mesh->nodes[nx + 1][j].potential = EDGE[direction_t.RIGHT][j][1];
        }

    }

    return 0;

}
