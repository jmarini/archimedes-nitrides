/* hole_relaxation.c -- This file is part of GNU/Archimedes 0.0.3
   This code is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means
   of effective potential method.

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
   USA.
*/


#ifndef ARCHIMEDES_MEP_HOLE_RELAXATION_H
#define ARCHIMEDES_MEP_HOLE_RELAXATION_H


#include "mesh.h"


void MEP_hole_relaxation(Mesh *mesh) {
    int ND = 2;
    int nx = mesh->nx,
        ny = mesh->ny;

    // Euler Relaxation step
    for(int c = 1; c <= 50; ++c) {
        for(int i = ND; i <= nx + ND; ++i) {
            for(int j = ND; j <= ny + ND; ++j) {
                double t = (2. / 3.) * (h2d[i][j][4] / h2d[i][j][1]) / KB;
                double ktaup = M * mstarhole * MIU0hole * g_config->lattice_temp / Q;
                double ktauw = MIU0hole * g_config->lattice_temp * KB
                             / (Q * VShole * VShole);
                double taup = ktaup / t;
                double tauw = 0.5 * ktaup / t
                            + 1.5 * ktauw * t / (t + g_config->lattice_temp);
                h2d[i][j][4] += -g_config->dt / 50.
                              * (-Q * (h2d[i][j][2] * g_mesh->nodes[i-1][j-1].efield.x +
                                       h2d[i][j][3] * g_mesh->nodes[i-1][j-1].efield.y)
                                 + h2d[i][j][1] * 1.5 * KB * (t - g_config->lattice_temp) / tauw);
                h2d[i][j][2] += -g_config->dt / 50.
                              * (-h2d[i][j][1] * Q * g_mesh->nodes[i-1][j-1].efield.x / (M * mstarhole)
                                 + h2d[i][j][2] / taup);
                h2d[i][j][3] += -g_config->dt / 50
                              * (-h2d[i][j][1] * Q * g_mesh->nodes[i-1][j-1].efield.y / (M * mstarhole)
                                 + h2d[i][j][3] / taup);
            }
        }
    }
}

#endif
