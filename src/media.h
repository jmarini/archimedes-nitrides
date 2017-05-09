/* media.h -- This file is part of Archimedes release 1.5.0
   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

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

#include "mesh.h"


void media(Mesh *mesh, int iteration) {
    printf("Computation of macroscopic observables\n");

    int i = 0,
        j = 0,
        n = 0;
    int density[NXM+1][NYM+1];
    real xvel[NXM+1][NYM+1],
         yvel[NXM+1][NYM+1],
         ener[NXM+1][NYM+1];

    // resetting of the electronic density
    // a simple way to avoid NaN propagation...
    memset(density, 0, sizeof(density[0][0]) * (NXM + 1) * (NYM + 1));
    memset(xvel,    0, sizeof(xvel[0][0])    * (NXM + 1) * (NYM + 1));
    memset(yvel,    0, sizeof(yvel[0][0])    * (NXM + 1) * (NYM + 1));
    memset(ener,    0, sizeof(ener[0][0])    * (NXM + 1) * (NYM + 1));

    // calculate info for each particle
    for(n = 1; n <= g_config->num_particles; n++) {
        particle_info_t info = mc_calculate_particle_info(&(mesh->particles[n]));
        i = info.i;
        j = info.j;

        density[i][j]++;
        ener[i][j] += info.energy;
        ener[i][j] += material_node(i, j).cb.emin[info.valley];
        xvel[i][j] += info.vx;
        yvel[i][j] += info.vy;
    }

    // Mean Value of the macroscopic variables
    // =======================================

    // Average velocity and energy over grid cell
    for(i = 1; i <= g_mesh->nx + 1; i++) {
        for(j = 1; j <= g_mesh->ny + 1; j++) {
            if(density[i][j] != 0) {
                xvel[i][j] /= (real)density[i][j];
                yvel[i][j] /= (real)density[i][j];
                ener[i][j] /= (real)density[i][j];
            }
        }
    }

    for(i = 1; i <= g_mesh->nx + 1; i++) {
        for(j = 1; j <= g_mesh->ny + 1; j++){
            u2d[i][j][2] += xvel[i][j];
            u2d[i][j][3] += yvel[i][j];
            u2d[i][j][4] += ener[i][j];
            if(g_config->time < g_config->dt) {
                moving_average[i][j][2] = xvel[i][j];
                moving_average[i][j][3] = yvel[i][j];
                moving_average[i][j][4] = ener[i][j];
            }
            else {
                moving_average[i][j][2] = g_config->avg_alpha * xvel[i][j] + (1. - g_config->avg_alpha) * moving_average[i][j][2];
                moving_average[i][j][3] = g_config->avg_alpha * yvel[i][j] + (1. - g_config->avg_alpha) * moving_average[i][j][3];
                moving_average[i][j][4] = g_config->avg_alpha * ener[i][j] + (1. - g_config->avg_alpha) * moving_average[i][j][4];
            }
        }
    }
}
