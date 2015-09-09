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


// ######################################################
// Created on 5 aug.2004, Siracusa, J.M.Sellier
// Last modif. : 17 aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Universal Physical Constants in M.K.S.C. system

// ######################################################
// ######################################################

// Boltzmann constant (Joule/Kelvin)
 const real KB=1.380658e-23;
// Electron charge in absolute value (Coulomb)
 const real Q=1.60217733e-19;
// Reduced Planck constant (Joule*sec)
 const real HBAR=1.05457266e-34;
// Permittivity of free space (F/m)
 const real EPS0=8.854187817e-12;
// Electron Mass (Kg)
 const real M=9.1093897e-31;
// Silicon intrinsic density for room temperature
 const real NI=1.45e16;
// Pi number
 const real PI=3.141592654;
// electron energy step (eV) for the MC method
 const real DE=0.002;
// Silicon low field mobility (m^2/(V*sec))
 real MIU0=1400.e-4;
// Silicon saturation velocity (m/sec)
 real VS=1.e5;
// Silicon heavy hole effective mass
 real mstarhole=0.57;
// Silicon heavy hole low field mobility (m^2/(V*sec))
 real MIU0hole=0.0471;
// Silicon heavy hole saturation velocity (m/sec)
 real VShole=1.e5;
// Silicon Schottky contact density (1/m^3)
 real NGATE=3.9e11;
// ######################################################
// ######################################################
