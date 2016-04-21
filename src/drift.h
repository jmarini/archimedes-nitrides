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


// ######################################################
// Created on 06 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 31 Aug.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// calculation of drift process

void drift(particle_t *particle, real tau)
{
    int iaux;
    int i,j;
    real dkx,dky,hmt,ksquared;
    real vx,vy;

    if(!mc_does_particle_exist(particle)) { return; }

    mc_particle_coords(particle, &i, &j);

    if(NOVALLEY[i_dom[i][j]] == 1) { iaux = 0; }
    if(NOVALLEY[i_dom[i][j]] >= 2) { iaux = particle->valley; }

    // Electron drift process
    // second order Runge-Kutta method
    hmt = HM[i_dom[i][j]][iaux] * tau;
    ksquared = mc_particle_ksquared(particle);

    if(CONDUCTION_BAND == KANE) {
        real thesquareroot, gk;
        gk = HHM[i_dom[i][j]][iaux] * ksquared;
        thesquareroot = sqrt(1. + 4. * alphaK[i_dom[i][j]][particle->valley] * gk);
        vx = particle->kx * HM[i_dom[i][j]][iaux] / thesquareroot;
        vy = particle->ky * HM[i_dom[i][j]][iaux] / thesquareroot;
        dkx = -QH * (E[i][j][0] + vy * B[i][j]) * tau;
        dky = -QH * (E[i][j][1] - vx * B[i][j]) * tau;
        particle->x += hmt * (particle->kx + 0.5 * dkx) / thesquareroot;
        particle->y += hmt * (particle->ky + 0.5 * dky) / thesquareroot;
        particle->kx += dkx;
        particle->ky += dky;
    }
    else if(CONDUCTION_BAND == PARABOLIC) {
        vx = particle->kx * HM[i_dom[i][j]][iaux];
        vy = particle->ky * HM[i_dom[i][j]][iaux];
        dkx = -QH * (E[i][j][0] + vy * B[i][j]) * tau;
        dky = -QH * (E[i][j][1] - vx * B[i][j]) * tau;
        particle->x += hmt * (particle->kx + 0.5 * dkx);
        particle->y += hmt * (particle->ky + 0.5 * dky);
        particle->kx += dkx;
        particle->ky += dky;
    }
    else if(CONDUCTION_BAND == FULL) {
        real k4, k2, ks;
        real dx, dy, d;
        vx = particle->kx * HM[i_dom[i][j]][iaux];
        vy = particle->ky * HM[i_dom[i][j]][iaux];
        dkx = -QH * (E[i][j][0] + vy * B[i][j]) * tau;
        dky = -QH * (E[i][j][1] - vx * B[i][j]) * tau;
        k2 = (particle->kx + 0.5 * dkx) * (particle->kx + 0.5 * dkx)
           + (particle->ky + 0.5 * dky) * (particle->ky + 0.5 * dky)
           +  particle->kz              *  particle->kz;
        ks = sqrt(k2) * 1.e-12 * 0.5 / PI;
        k2 = ks * ks;
        k4 = k2 * k2;
        d = 10. * CB_FULL[i_dom[i][j]][0] * k4 * k4 * ks
          +  9. * CB_FULL[i_dom[i][j]][1] * k4 * k4
          +  8. * CB_FULL[i_dom[i][j]][2] * k4 * k2 * ks
          +  7. * CB_FULL[i_dom[i][j]][3] * k4 * k2
          +  6. * CB_FULL[i_dom[i][j]][4] * k4 * ks
          +  5. * CB_FULL[i_dom[i][j]][5] * k4
          +  4. * CB_FULL[i_dom[i][j]][6] * k2 * ks
          +  3. * CB_FULL[i_dom[i][j]][7] * k2
          +  2. * CB_FULL[i_dom[i][j]][8] * ks
          +       CB_FULL[i_dom[i][j]][9];
        ks *= 1.e+12 * 2.  * PI;
        d  *= 1.e-12 * 0.5 / PI;
        dx = QH * d * tau * (particle->kx + 0.5 * dkx) / ks;
        dy = QH * d * tau * (particle->ky + 0.5 * dky) / ks;
        particle->kx += dkx;
        particle->ky += dky;
        particle->x += dx;
        particle->y += dy;
    }

    // check if some particles are out of the device
    mc_particle_coords(particle, &i, &j);


    // Generic boundary conditions for the super-particles
    // ===================================================

    // left edge
    // =========
    // ---Insulator---
    if(particle->x <= 0. && mc_is_boundary_insulator(direction_t.LEFT, j)) {
        particle->x  *= -1.;
        particle->kx *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->x <= 0. && mc_is_boundary_contact(direction_t.LEFT, j)) {
        mc_remove_particle(particle);
        return;
    }

    // right edge
    // ==========
    // ---Insulator---
    if(particle->x >= LX && mc_is_boundary_insulator(direction_t.RIGHT, j)) {
        particle->x = LX - (particle->x - LX);
        particle->kx *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->x >= LX && mc_is_boundary_contact(direction_t.RIGHT, j)) {
        mc_remove_particle(particle);
        return;
    }

    // bottom edge
    // ===========
    // ---Insulator---
    if(particle->y <= 0. && mc_is_boundary_insulator(direction_t.BOTTOM, i)) {
        particle->y  *= -1.;
        particle->ky *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->y <= 0. && mc_is_boundary_contact(direction_t.BOTTOM, i)) {
        mc_remove_particle(particle);
        return;
    }

    // upper edge
    // ==========
    // ---Insulator---
    if(particle->y >= LY && mc_is_boundary_insulator(direction_t.TOP, i)) {
        particle->y = LY - (particle->y - LY);
        particle->ky *= -1.;
        return;
    }
    // ---Schottky or ohmic contact---
    else if(particle->y >= LY && mc_is_boundary_contact(direction_t.TOP, i)) {
        mc_remove_particle(particle);
        return;
    }

}

// ============================================================
