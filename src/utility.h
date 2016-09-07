/* utility.h -- This file is part of Archimedes release 1.2.0.
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
// Created on 21 oct.2015, J. Marini
// Last modif. : 21 oct.2015, J. Marini
// ######################################################

#ifndef ARCHIMEDES_UTILITY_H
#define ARCHIMEDES_UTILITY_H

#include "particle.h"
#include "mesh.h"


struct {
    int BOTTOM;
    int RIGHT;
    int TOP;
    int LEFT;
} direction_t = {.BOTTOM=0, .RIGHT=1, .TOP=2, .LEFT=3};


struct {
    int INSULATOR;
    int SCHOTTKY;
    int OHMIC;
} boundary_t = {.INSULATOR=0, .SCHOTTKY=1, .OHMIC=2};


inline int mc_is_boundary_insulator(int direction, int index) {
    return EDGE[direction][index][0] == boundary_t.INSULATOR;
}


inline int mc_is_boundary_schottky(int direction, int index) {
    return EDGE[direction][index][0] == boundary_t.SCHOTTKY;
}


inline int mc_is_boundary_ohmic(int direction, int index) {
    return EDGE[direction][index][0] == boundary_t.OHMIC;
}


inline int mc_is_boundary_contact(int direction, int index) {
    return mc_is_boundary_schottky(direction, index)
        || mc_is_boundary_ohmic(direction, index);
}


real mc_particle_energy(Particle *particle) {
    Node *node = mc_get_particle_node(particle);

    if(g_config->conduction_band == PARABOLIC) {
        return HHM[node->material][particle->valley] * mc_particle_ksquared(particle);
    }
    else if(g_config->conduction_band == KANE) {
        real alpha = alphaK[node->material][particle->valley];
        real gamma = HHM[node->material][particle->valley] * mc_particle_ksquared(particle);
        return (-1.0 + sqrt(1.0 + 4.0 * alpha * gamma)) / (2.0 * alpha);
    }
    else {
        return -1.0;
    }
}


inline char* mc_material_name(int material) {
    switch(material) {
        case iSIO2: return "iSIO2";
        case SILICON: return "Si";
        case GAAS: return "GaAs";
        case GERMANIUM: return "Ge";
        case INSB: return "InSb";
        case ALSB: return "AlSb";
        case ALXINXSB: return "AlxInxSb";
        case ALXIN1XSB: return "AlxIn1xSb";
        case ALAS: return "AlAs";
        case ALP: return "AlP";
        case GAP: return "GaP";
        case GASB: return "GaSb";
        case INAS: return "InAs";
        case INP: return "InP";
        case INXGA1XAS: return "InxGa1xAs";
        case INXAL1XAS: return "InxAl1xAs";
        case INXGAXXAS: return "InxGaXxAs";
        case GAN: return "GaN";
        default: return "Unknown Material";
    }
}


inline char* mc_band_model_name(int model) {
    switch(model) {
        case PARABOLIC: return "Parabolic";
        case KANE: return "Kane";
        case FULL: return "Full";
        default: return "Unknown Band Model";
    }
}


particle_info_t mc_calculate_particle_info(Particle *p) {
    // calculate particle coordinates
    int i = 0,
        j = 0;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;

    i = (int)(p->x / dx + 1.5);
    if(i <= 1) { i = 1; }
    if(i >= nx + 1) { i = nx + 1; }

    j = (int)(p->y / dy + 1.5);
    if(j <= 1) { j = 1; }
    if(j >= ny + 1) { j = ny + 1; }

    int material = g_mesh->info[i][j].material;

    // calculate particle energy and velocity
    real ksquared = mc_particle_ksquared(p);
    real energy = 0.;
    real xvelocity = 0.,
         yvelocity = 0.;

    if(g_config->conduction_band == PARABOLIC) {
        energy = HHM[material][p->valley] * ksquared;
        xvelocity = p->kx * HM[material][p->valley];
        yvelocity = p->ky * HM[material][p->valley];
    }
    else if(g_config->conduction_band == KANE) {
        real sq = sqrt(1. + 4. * alphaK[material][p->valley] * HHM[material][p->valley] * ksquared);
        energy = (sq - 1.) / (2. * alphaK[material][p->valley]);
        xvelocity = p->kx * HM[material][p->valley] / sq;
        yvelocity = p->ky * HM[material][p->valley] / sq;
    }
    else if(g_config->conduction_band == FULL) {
        real k, k2, k4;
        real d;
        k = sqrt(ksquared) * 0.5 / PI * 1.e-12;
        // periodicity on reciprocal lattice
        k2 = k * k;
        k4 = k2 * k2;
        energy = CB_FULL[material][0] * k4 * k4 * k2
               + CB_FULL[material][1] * k4 * k4 * k
               + CB_FULL[material][2] * k4 * k4
               + CB_FULL[material][3] * k4 * k2 * k
               + CB_FULL[material][4] * k4 * k2
               + CB_FULL[material][5] * k4 * k
               + CB_FULL[material][6] * k4
               + CB_FULL[material][7] * k2 * k
               + CB_FULL[material][8] * k2
               + CB_FULL[material][9] * k
               + CB_FULL[material][10]; // in eV

        d = 10. * CB_FULL[material][0] * k4 * k4 * k
          +  9. * CB_FULL[material][1] * k4 * k4
          +  8. * CB_FULL[material][2] * k4 * k2 * k
          +  7. * CB_FULL[material][3] * k4 * k2
          +  6. * CB_FULL[material][4] * k4 * k
          +  5. * CB_FULL[material][5] * k4
          +  4. * CB_FULL[material][6] * k2 * k
          +  3. * CB_FULL[material][7] * k2
          +  2. * CB_FULL[material][8] * k
          +       CB_FULL[material][9];
        k *= 1.e+12 * 2. * PI;
        d *= 1.e-12 * 0.5 / PI;
        xvelocity = QH * d * p->kx / k;
        yvelocity = QH * d * p->ky / k;
    }


    return (particle_info_t){
        .id=p->id,
        .valley=p->valley,
        .kx=p->kx,
        .ky=p->ky,
        .kz=p->kz,
        .energy=energy,
        .t=p->t,
        .x=p->x,
        .y=p->y,
        .i=i,
        .j=j,
        .vx=xvelocity,
        .vy=yvelocity
    };
}


#endif
