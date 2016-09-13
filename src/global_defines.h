/* global_defines.h -- This file is part of Archimedes release 1.2.0.
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


#ifndef ARCHIMEDES_DEFINES_H
#define ARCHIMEDES_DEFINES_H

#define real double
#define ON 1
#define OFF 0
#define INVALID -1             // invalid value for most enum-like variables
#define KANE 0                 // conduction band model, kane model
#define PARABOLIC 1            // conduction band model, parabolic approximation
#define FULL 2                 // conduction band model, full band model
#define QEP_BOHM 0             // quantum effective potential, bohm potential
#define QEP_CALIBRATED_BOHM 1  // quantum effective potential, calibrated bohm potential
#define QEP_FULL 2             // quantum effective potential, full effective potential
#define QEP_DENSITY_GRADIENT 3 // quantum effective potential, density gradient
#define MN3 4                  // number of summary values to save per cell
#define DIME 3003              // maximum number of points in energy mesh
#define ITMAX 10000000         // maximum number of monte carlo iterations
#define POISSONITMAX 1500      // maximum number of poisson iterations
#define SMALL 1.e-5            // defines what is a "small" number/delta
#define VMAX 1000000
#define NPMAX 10000000         // maximum number of super-particles
#define MCE 0                  // MCE stands for MC for electrons only
#define MCH 1                  // MCH stands for MC for holes only
#define MCEH 2                 // MCEH stands for MC both for electrons and holes
#define MEPE 3                 // MEPE stands for MEP model for electrons only
#define MEPH 4                 // MEPH stands for MEP model for holes only
#define MEPEH 5                // MEPEH stands for MEP model for electrons and holes
#define GNUPLOTFORMAT 0        // output file in GNUPLOT format
#define MESHFORMAT 1           // output file in Mesh format
#define MAX_VALLEYS 4          // maximum number of valleys

// definition of the material reference table
#define NOAMTIA 17   // Number Of All Material Taken Into Account (excluding SiO2)
// ********************
#define iSIO2 -1     // reference number to Silicon Oxide
#define SILICON 0    // reference number to Silicon
#define GAAS 1       // reference number to GaAs
#define GERMANIUM 2  // reference number to Germanium
#define INSB 3       // reference number to InSb
#define ALSB 4       // reference number to AlSb
#define ALXINXSB 5   // reference number to Al_x In_x Sb
#define ALXIN1XSB 6  // reference number to Al_x In_(1-x) Sb
#define ALAS 7       // reference number to AlAs
#define ALP 8        // reference number to AlP
#define GAP 9        // reference number to GaP
#define GASB 10      // reference number to GaSb
#define INAS 11      // reference number to InAs
#define INP 12       // reference number to InP
#define INXGA1XAS 13 // reference number to In_x Ga_(1-x) As
#define INXAL1XAS 14 // reference number to In_x Al_(1-x) As
#define INXGAXXAS 15 // reference number to In_x Ga_(1-x) As (second zone)
#define GAN 16       // reference number to GaN
// ********************

#define NUMSIO2 2 // maximum number of SiO2 interfaces

#endif
