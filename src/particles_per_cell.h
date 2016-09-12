/* particles_per_cell.h -- This file is part of Archimedes release 0.1.2
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for
   electrons. It also includes the quantum effects by means
   of Bohm effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>

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


#include "mesh.h"


// calculate electron density per cell using particle in cell method
int calculate_particles_per_cell(void) {
    int nx  = g_mesh->nx,
        ny  = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;


    // resetting of the electronic density
    // a simple way to avoid NaN propagation...
    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 1; j <= ny + 1; ++j) {
            mc_node(i, j)->e.density = 0;
        }
    }

    // cloud in cell method
    for(int n = 1; n <= g_config->num_particles; ++n) {
        real x = P[n].x / dx;
        real y = P[n].y / dy;
        int i = (int)(x + 1.);
        int j = (int)(y + 1.);

        real x1 = (real)i - x;
        real y1 = (real)j - y;
        real x2  = x - (real)(i - 1);
        real y2  = y - (real)(j - 1);

        mc_node(i, j)->e.density += x1 * y1;
        if(i <= nx) {
            mc_node(i + 1, j)->e.density += x2 * y1;
        }
        if(j <= ny) {
            mc_node(i, j + 1)->e.density += x1 * y2;
        }
        if(i <= nx && j <= ny) {
            mc_node(i + 1, j + 1)->e.density += x2 * y2;
        }
    }

    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 1; j <= ny + 1; ++j) {
            mc_node(i, j)->e.density *= g_config->carriers_per_superparticle / (dx * dy);
            if(i == 1 || i == nx + 1) { mc_node(i, j)->e.density *= 2.; }
            if(j == 1 || j == ny + 1) { mc_node(i, j)->e.density *= 2.; }
        }
    }
    mc_node(nx + 1, ny + 1)->e.density = mc_node(nx, ny + 1)->e.density;
    mc_node(     1, ny + 1)->e.density = mc_node( 1, ny    )->e.density;

    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 1; j <= ny + 1; ++j) {
            u2d[i][j][1] = mc_node(i, j)->e.density;
        }
    }

    return 0;
}
