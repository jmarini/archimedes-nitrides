/* ensemblemontecarlo.h -- This file is part of Archimedes release 1.2.0.
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
   but WITHOUT ANparticle->y WARRANTparticle->y; without even the implied warranty of
   MERCHANTABILITparticle->y or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   particle->you should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 25 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Ensemble Monte Carlo method

void EMC(void)
{
    EMC_threaded();

    printf("Actual number of electron super-particles = %d\n", INUM);
    if(INUM > NPMAX) {
        printf("%s: too big actual number of particles\n", progname);
        exit(EXIT_FAILURE);
    }
}

// ============================================================
