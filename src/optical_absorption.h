/* optical_absorption.h -- This file is part of Archimedes release 1.2.0.
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
// Created on 16 jun.2016, J. Marini
// Last modif. : 16 jun.2016, J. Marini
// ######################################################


real optical_joint_DOS(int material, int conduction_band, int valence_band, real energy) {
    real mc = MSTAR[material][conduction_band],
         mv = MSTAR_VB[material][valence_band];
    real mr = mc * mv / (mc + mv);
    real alpha = alphaK[material][conduction_band];

    real deltaE = DELTAE_VB[material][valence_band] + EG[material];
    real e2 = (mr / mc) * (energy - deltaE);

    real gamma = Q * e2 * (1. + alpha * e2);
    real gamma2 = (1. + 2. * alpha * e2);

    real dos = mr * M * sqrt(2. * mr * M) * sqrt(gamma) * gamma2 / (PI * PI * HBAR * HBAR * HBAR);
    return dos;
}


real absorption_coefficient(int material) {
    return 0.0;
}


real optical_transition_rate(int material, int conduction_band,
                             int valence_band, real photon_energy) {
    real prefactor = PI * HBAR * Q * Q * EG[material]
                   / (3. * EPSR[material] * EPS0 * M * photon_energy);
    real rate = prefactor * optical_joint_DOS(material, conduction_band, valence_band, photon_energy);
    return rate;
}


int electrons_in_cell(int material, int conduction_band, int valence_band, real photon_energy, int i, int j) {
    return 0;
}
