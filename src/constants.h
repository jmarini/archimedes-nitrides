/* constants.h -- This file is part of Archimedes release 1.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements both the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It includes some quantum effects by means
   of effective potential method. It is also able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier <sellier@dmi.unict.it>

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


#ifndef ARCHIMEDES_CONSTANTS_H
#define ARCHIMEDES_CONSTANTS_H


// Universal Physical Constants in M.K.S.C. system


static const double KB     = 1.380658e-23;    // Boltzmann constant [J/K]
static const double Q      = 1.60217733e-19;  // electron charge (absolute value) [C]
static const double HBAR   = 1.05457266e-34;  // reduced Planck constant [Js]
static const double EPS0   = 8.854187817e-12; // permittivity of free space [F/m]
static const double M      = 9.1093897e-31;   // electron Mass [kg]
static const double VLIGHT = 2.99792458e8;    // speed of light (m/s)

static const double PI     = 3.141592654;

static const double DE     = 0.002;           // energy step [eV]


// Silicon intrinsic density for room temperature
static const double NI=1.45e16;

// Silicon Schottky contact density (1/m^3)
static const double NGATE=3.9e11;

#endif
