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

#ifndef ARCHIMEDES_MEP_CONSTANTS_H
#define ARCHIMEDES_MEP_CONSTANTS_H


// Universal Physical Constants in M.K.S.C. system

// Silicon low field mobility (m^2/(V*sec))
static const double MIU0=1400.e-4;

// Silicon saturation velocity (m/sec)
static const double VS=1.e5;

// Silicon heavy hole effective mass
static const double mstarhole=0.57;

// Silicon heavy hole low field mobility (m^2/(V*sec))
static const double MIU0hole=0.0471;

// Silicon heavy hole saturation velocity (m/sec)
static const double VShole=1.e5;

#endif
