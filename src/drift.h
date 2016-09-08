/* drift.h -- This file is part of Archimedes release 1.2.0
   Archimedes is a simulator for Submicron 2D III-V Semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means
   of effective potential method. It is now able to simulate applied
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


#include "mesh.h"
#include "particle.h"
#include "vec.h"


// calculation of drift process over time tau
void drift(Particle *particle, real tau) {
    Vec2 dk = {0., 0.};
    Vec2 v = {0., 0.};

    if(!mc_does_particle_exist(particle)) { return; }

    Node *node = mc_get_particle_node(particle);
    int material = node->material;

    // Electron drift process
    // second order Runge-Kutta method
    real hmt = HM[material][particle->valley] * tau;
    real ksquared = mc_particle_ksquared(particle);

    if(g_config->conduction_band == KANE) {
        real thesquareroot, gk;
        gk = HHM[material][particle->valley] * ksquared;
        thesquareroot = sqrt(1. + 4. * node->mat->cb.alpha[particle->valley] * gk);
        v.x = particle->kx * HM[material][particle->valley] / thesquareroot;
        v.y = particle->ky * HM[material][particle->valley] / thesquareroot;
        dk.x = -QH * (node->efield.x + v.y * node->magnetic_field) * tau;
        dk.y = -QH * (node->efield.y - v.x * node->magnetic_field) * tau;
        particle->x += hmt * (particle->kx + 0.5 * dk.x) / thesquareroot;
        particle->y += hmt * (particle->ky + 0.5 * dk.y) / thesquareroot;
        particle->kx += dk.x;
        particle->ky += dk.y;
    }
    else if(g_config->conduction_band == PARABOLIC) {
        v.x = particle->kx * HM[material][particle->valley];
        v.y = particle->ky * HM[material][particle->valley];
        dk.x = -QH * (node->efield.x + v.y * node->magnetic_field) * tau;
        dk.y = -QH * (node->efield.y - v.x * node->magnetic_field) * tau;
        particle->x += hmt * (particle->kx + 0.5 * dk.x);
        particle->y += hmt * (particle->ky + 0.5 * dk.y);
        particle->kx += dk.x;
        particle->ky += dk.y;
    }
    else if(g_config->conduction_band == FULL) {
        real k4, k2, ks;
        real dx, dy, d;
        v.x = particle->kx * HM[material][particle->valley];
        v.y = particle->ky * HM[material][particle->valley];
        dk.x = -QH * (node->efield.x + v.y * node->magnetic_field) * tau;
        dk.y = -QH * (node->efield.y - v.x * node->magnetic_field) * tau;
        k2 = (particle->kx + 0.5 * dk.x) * (particle->kx + 0.5 * dk.x)
           + (particle->ky + 0.5 * dk.y) * (particle->ky + 0.5 * dk.y)
           +  particle->kz               *  particle->kz;
        ks = sqrt(k2) * 1.e-12 * 0.5 / PI;
        k2 = ks * ks;
        k4 = k2 * k2;
        d = 10. * CB_FULL[material][0] * k4 * k4 * ks
          +  9. * CB_FULL[material][1] * k4 * k4
          +  8. * CB_FULL[material][2] * k4 * k2 * ks
          +  7. * CB_FULL[material][3] * k4 * k2
          +  6. * CB_FULL[material][4] * k4 * ks
          +  5. * CB_FULL[material][5] * k4
          +  4. * CB_FULL[material][6] * k2 * ks
          +  3. * CB_FULL[material][7] * k2
          +  2. * CB_FULL[material][8] * ks
          +       CB_FULL[material][9];
        ks *= 1.e+12 * 2.  * PI;
        d  *= 1.e-12 * 0.5 / PI;
        dx = QH * d * tau * (particle->kx + 0.5 * dk.x) / ks;
        dy = QH * d * tau * (particle->ky + 0.5 * dk.y) / ks;
        particle->kx += dk.x;
        particle->ky += dk.y;
        particle->x += dx;
        particle->y += dy;
    }

    // check if some particles are out of the device
    node = mc_get_particle_node(particle);


    // Generic boundary conditions for the super-particles
    // ===================================================

    // left edge
    // =========
    // ---Insulator---
    if(particle->x <= 0. && mc_is_boundary_insulator(direction_t.LEFT, node->j)) {
        particle->x  *= -1.;
        particle->kx *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->x <= 0. && mc_is_boundary_contact(direction_t.LEFT, node->j)) {
        mc_remove_particle(particle);
        return;
    }

    // right edge
    // ==========
    // ---Insulator---
    if(particle->x >= g_mesh->width && mc_is_boundary_insulator(direction_t.RIGHT, node->j)) {
        particle->x = g_mesh->width - (particle->x - g_mesh->width);
        particle->kx *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->x >= g_mesh->width && mc_is_boundary_contact(direction_t.RIGHT, node->j)) {
        mc_remove_particle(particle);
        return;
    }

    // bottom edge
    // ===========
    // ---Insulator---
    if(particle->y <= 0. && mc_is_boundary_insulator(direction_t.BOTTOM, node->i)) {
        particle->y  *= -1.;
        particle->ky *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->y <= 0. && mc_is_boundary_contact(direction_t.BOTTOM, node->i)) {
        mc_remove_particle(particle);
        return;
    }

    // upper edge
    // ==========
    // ---Insulator---
    if(particle->y >= g_mesh->height && mc_is_boundary_insulator(direction_t.TOP, node->i)) {
        particle->y = g_mesh->height - (particle->y - g_mesh->height);
        particle->ky *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->y >= g_mesh->height && mc_is_boundary_contact(direction_t.TOP, node->i)) {
        mc_remove_particle(particle);
        return;
    }

}
