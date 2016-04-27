/* material_parameters.h -- This file is part of Archimedes release 1.2.0.
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
// Created on 17 apr.2016, J. Marini
// Last modif. : 17 apr.2016, J. Marini
// ######################################################


// =========================
// Band Structure Parameters
// =========================

// Number of valleys
NOVALLEY[SILICON]   = 1;  // X-valley
NOVALLEY[GERMANIUM] = 1;  // G-valley
NOVALLEY[GAAS]      = 2;  // G and L-valleys
NOVALLEY[INSB]      = 1;  // G valley
NOVALLEY[ALSB]      = 1;  // G-valley
NOVALLEY[ALXINXSB]  = 1;  // G-valley
NOVALLEY[ALXIN1XSB] = 1;  // G-valley
NOVALLEY[ALAS]      = 1;  // G-valley
NOVALLEY[ALP]       = 1;  // G-valley
NOVALLEY[GAP]       = 1;  // G-valley
NOVALLEY[GASB]      = 1;  // G-valley
NOVALLEY[INAS]      = 1;  // G-valley
NOVALLEY[INP]       = 1;  // G-valley
NOVALLEY[INXGA1XAS] = 1;  // only G-valley
NOVALLEY[INXAL1XAS] = 1;  // G-valley
NOVALLEY[INXGAXXAS] = 1;  // only G-valley
NOVALLEY[GAN]       = 2;  // G-1, M-L(U), G-3


// Number of equivalent valleys
// Scattering from first index to second index
// Repeated indexes mean scattering between equivalent valleys
memset(ZSCATTER, 0, sizeof(ZSCATTER[0][0][0]) * NOAMTIA * 6 * 6);

ZSCATTER[GAAS][1][1] = 0; // G -> G
ZSCATTER[GAAS][1][2] = 4; // G -> L
ZSCATTER[GAAS][1][3] = 3; // G -> X
ZSCATTER[GAAS][2][1] = 1; // L -> G
ZSCATTER[GAAS][2][2] = 3; // L -> L
ZSCATTER[GAAS][2][3] = 3; // L -> X
ZSCATTER[GAAS][3][1] = 1; // X -> G
ZSCATTER[GAAS][3][2] = 4; // X -> L
ZSCATTER[GAAS][3][3] = 2; // X -> X

ZSCATTER[GAN][1][1] = 0; // G1 -> G1
ZSCATTER[GAN][1][2] = 6; // G1 -> ML
ZSCATTER[GAN][1][3] = 1; // G1 -> G3
ZSCATTER[GAN][2][1] = 1; // ML -> G1
ZSCATTER[GAN][2][2] = 5; // ML -> ML
ZSCATTER[GAN][2][3] = 1; // ML -> G3
ZSCATTER[GAN][3][1] = 1; // G3 -> G1
ZSCATTER[GAN][3][2] = 6; // G3 -> ML
ZSCATTER[GAN][3][3] = 0; // G3 -> G3


// Band minimum energy  - relative to CBM (eV)
// first valley
EMIN[SILICON][1]   = 0.0;    // Sellier, Fischetti, etc.
EMIN[GERMANIUM][1] = 0.173;
EMIN[GAAS][1]      = 0.0;    // Tomizawa
EMIN[INSB][1]      = 0.0;
EMIN[ALSB][1]      = 0.507;
EMIN[ALAS][1]      = 0.767;
EMIN[ALP][1]       = 1.237;
EMIN[GAP][1]       = 0.496;
EMIN[GASB][1]      = 0.0;
EMIN[INAS][1]      = 0.0;
EMIN[INP][1]       = 0.0;
EMIN[GAN][1]       = 0.0;    // G-1
// second valley
EMIN[GAAS][2]      = 0.323;  // L
EMIN[GAN][2]       = 2.1;    // L-M
// third valley
EMIN[GAAS][3]      = 0.48;   // X
EMIN[GAN][3]       = 1.9;    // G-2


// Definition of effective mass for all materials in all valleys
MSTAR[SILICON][1]   = 0.32;    // see Sellier, Tomizawa, etc.
MSTAR[GAAS][1]      = 0.067;   // Gamma-valley -- see Tomizawa
MSTAR[GAAS][2]      = 0.350;   // L-valley     -- see Tomizawa
MSTAR[GAAS][3]      = 0.27;    // X-valley
MSTAR[GERMANIUM][1] = 0.12;    // Gamma valley -- see http://ecee.colorado.edu/~bart/book/effmass.htm#long
MSTAR[INSB][1]      = 0.0135;  // Gamma-valley -- see Ram-Mohan
MSTAR[ALSB][1]      = 0.14;    // Gamma-valley -- See Ram-Mohan
MSTAR[ALAS][1]      = 0.149;   // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
MSTAR[ALP][1]       = 0.22;    // Gamma-valley -- see Ram-Mohan
MSTAR[GAP][1]       = 0.13;    // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
MSTAR[GASB][1]      = 0.039;   // Gamma-valley -- see Ram-Mohan
MSTAR[INAS][1]      = 0.026;   // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
MSTAR[INP][1]       = 0.0795;  // Gamma-valley -- see Ram-Mohan
MSTAR[GAN][1]       = 0.2;     // G-1 -- Foutz, O'Leary, Shur, Eastman
MSTAR[GAN][2]       = 0.4;     // L-M -- Bhapkar & Shur
MSTAR[GAN][3]       = 0.6;     // G-2 -- Bhapkar & Shur


// non-parabolicity coefficients (1/eV)
alphaK[SILICON][1]   = 0.5;    // see Sellier, Tomizawa
alphaK[GERMANIUM][1] = 0.3;  // Gamma valley - Jacoboni Reggiani
alphaK[GAN][1]       = 0.189;      // G-1 -- Foutz, O'Leary, Shur, Eastman
alphaK[GAN][2]       = 0.065;      // L-M -- Bhapkar & Shur
alphaK[GAN][3]       = 0.029;      // G-2 -- Bhapkar & Shur


// Dielectric constant for Silicon Oxide SiO2
// see http://en.wikipedia.org/wiki/Relative_permittivity
EPSRSIO2 = 3.9 * EPS0;

// Dielectric constant for Semiconducting materials
// STATIC
// ======
EPSR[SILICON]   = 11.68;      // see http://en.wikipedia.org/wiki/Relative_permittivity
EPSR[GERMANIUM] = 16.2;       // see http://www.ioffe.ru/SVA/NSM/Semicond/Ge/basic.html
EPSR[GAAS]      = 12.90;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
EPSR[INSB]      = 16.8;       // see http://www.ioffe.ru/SVA/NSM/Semicond/InSb/basic.html
EPSR[ALSB]      = 12.04;      // Fischetti conversations
EPSR[ALAS]      = 12.90-2.84; // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
EPSR[ALP]       = 9.80;       // Fischetti conversations
EPSR[GAP]       = 11.10;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
EPSR[GASB]      = 15.69;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
EPSR[INAS]      = 15.15;      // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
EPSR[INP]       = 12.50;      // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
EPSR[GAN]       = 9.7;        // E. Bellotti & F. Bertazzi
// HIGH FREQUENCY
// ==============
EPF[GAAS] = 10.89;       // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
EPF[INSB] = 15.68;       // Fischetti conversations
EPF[ALSB] = 9.88;        // Fiscehtti conversations
EPF[ALAS] = 10.89-2.73;  // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
EPF[ALP]  = 7.54;        // Fischetti conversations
EPF[GAP]  = 9.11;        // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
EPF[GASB] = 14.44;       // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
EPF[INAS] = 12.3;        // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
EPF[INP]  = 9.61;        // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
EPF[GAN]  = 5.28;        // E. Bellotti & F. Bertazzi

int ii;
for(ii = 0; ii < NOAMTIA; ii++) {
   int i;
   for(i = 0; i < 6; i++){
       HWO[ii][i] = 0.;
       DTK[ii][i] = 0.;
       ZF[ii][i]  = 0.;
   }
}

// Optical phonon scattering energy (eV)
HWO[SILICON][0]   = 0.0120;   // Sellier, Tomizawa
HWO[SILICON][1]   = 0.0185;   // Sellier, Tomizawa
HWO[SILICON][2]   = 0.0190;   // Sellier, Tomizawa
HWO[SILICON][3]   = 0.0474;   // Sellier, Tomizawa
HWO[SILICON][4]   = 0.0612;   // Sellier, Tomizawa
HWO[SILICON][5]   = 0.0590;   // Sellier, Tomizawa
HWO[GERMANIUM][0] = 0.03704;  // Fischetti
HWO[GAAS][0]      = 0.03536;  // Fischetti
HWO[INSB][0]      = 0.02404;  // Fischetti
HWO[ALSB][0]      = 0.0360;   // Fischetti
HWO[ALAS][0]      = 0.05009;  // Fischetti
HWO[ALP][0]       = 0.06211;  // Fischetti
HWO[GAP][0]       = 0.04523;  // Fischetti
HWO[GASB][0]      = 0.02529;  // Fischetti
HWO[INAS][0]      = 0.03008;  // Fischetti
HWO[INP][0]       = 0.04240;  // Fischetti
HWO[GAN][0]       = 0.09212;  // LO -- E. Bellotti & F. Bertazzi
HWO[GAN][1]       = 0.06955;  // TO -- E. Bellotti & F. Bertazzi

// Optical coupling constants (eV/m)
DTK[SILICON][0]   = 0.05e11;  // Jacoboni Reggiani
DTK[SILICON][1]   = 0.08e11;  // Jacoboni Reggiani
DTK[SILICON][2]   = 0.03e11;  // Jacoboni Reggiani
DTK[SILICON][3]   = 0.20e11;  // Jacoboni Reggiani
DTK[SILICON][4]   = 1.14e11;  // Jacoboni Reggiani
DTK[SILICON][5]   = 0.20e11;  // Jacoboni Reggiani
DTK[GERMANIUM][0] = 0.0;      // Jacoboni Reggiani
DTK[GERMANIUM][1] = 0.079e11; // Jacoboni Reggiani
DTK[GERMANIUM][2] = 0.0;      // Jacoboni Reggiani
DTK[GERMANIUM][3] = 0.0;      // Jacoboni Reggiani
DTK[GERMANIUM][4] = 0.95e11;  // Jacoboni Reggiani
DTK[GERMANIUM][5] = 0.0;      // Jacoboni Reggiani
DTK[GAAS][0]      = 1.11e11;  // Sellier, Tomizawa
DTK[INSB][0]      = 0.47e11;  // see ???
DTK[ALSB][0]      = 0.55e11;  // see ???
DTK[ALAS][0]      = 3.0e11;   // see ???
DTK[ALP][0]       = 0.95e11;  // see ???
DTK[GAP][0]       = 5.33e11;  // see ???
DTK[GASB][0]      = 0.94e11;  // see ???
DTK[INAS][0]      = 3.59e11;  // see ???
DTK[INP][0]       = 2.46e11;  // see ???
DTK[GAN][0]       = 1.0e11;   // E. Bellotti & F. Bertazzi
DTK[GAN][1]       = 1.0e11;   // E. Bellotti & F. Bertazzi

// Optical phonon Z-factor
ZF[SILICON][0]   = 1.;  // Sellier
ZF[SILICON][1]   = 1.;  // Sellier
ZF[SILICON][2]   = 4.;  // Sellier
ZF[SILICON][3]   = 4.;  // Sellier
ZF[SILICON][4]   = 1.;  // Sellier
ZF[SILICON][5]   = 4.;  // Sellier
ZF[GERMANIUM][0] = 1.;  // see ???
ZF[GAAS][0]      = 1.;  // Sellier
ZF[INSB][0]      = 1.;  // see ???
ZF[ALSB][0]      = 1.;  // see ???
ZF[ALAS][0]      = 1.;  // see ???
ZF[ALP][0]       = 1.;  // see ???
ZF[GAP][0]       = 1.;  // see ???
ZF[GASB][0]      = 1.;  // see ???
ZF[INAS][0]      = 1.;  // see ???
ZF[INP][0]       = 1.;  // see ???
ZF[GAN][0]       = 1.;  // guess for correct value
// ZF[GAN][1]       = 1.;  // guess for correct value

// Crystal Density (Kg/m^3)
RHO[SILICON]   = 2.33e3;  // Fischetti conversations
RHO[GERMANIUM] = 5.32e3;  // Fischetti conversations
RHO[GAAS]      = 5.36e3;  // Fischetti conversations
RHO[INSB]      = 5.78e3;  // Fischetti conversations
RHO[ALSB]      = 4.26e3;  // Fischetti conversations
RHO[ALAS]      = 3.76e3;  // Fischetti conversations
RHO[ALP]       = 2.40e3;  // Fischetti conversations
RHO[GAP]       = 4.14e3;  // Fischetti conversations
RHO[GASB]      = 5.61e3;  // Fischetti conversations
RHO[INAS]      = 5.67e3;  // Fischetti conversations
RHO[INP]       = 4.81e3;  // Fischetti conversations
RHO[GAN]       = 6.087e3; // E. Bellotti & F. Bertazzi

// Acoustic deformation potential (Joule)
DA[SILICON]   = 9.  * Q;  // Fischetti -- Jacoboni Reggiani
DA[GERMANIUM] = 9.  * Q;  // Fischetti -- Jacoboni Reggiani
DA[GAAS]      = 7.  * Q;  // Fischetti - Gamma valley
DA[INSB]      = 7.  * Q;  // Fischetti
DA[ALSB]      = 4.6 * Q;  // Fischetti
DA[ALAS]      = 9.3 * Q;  // Fischetti
DA[ALP]       = 9.3 * Q;  // Fischetti
DA[GAP]       = 7.4 * Q;  // Fischetti
DA[GASB]      = 9.  * Q;  // Fischetti
DA[INAS]      = 8.2 * Q;  // Fischetti
DA[INP]       = 6.2 * Q;  // Fischetti
DA[GAN]       = 8.3 * Q;  // E. Bellotti & F. Bertazzi

// Longitudinal sound velocity (m/sec)
UL[SILICON]   = 9.18e3;  // Fischetti
UL[GERMANIUM] = 5.4e3;   // Fischetti
UL[GAAS]      = 5.24e3;  // Fischetti
UL[INSB]      = 3.41e3;  // Fischetti
UL[ALSB]      = 4.25e3;  // Fischetti
UL[ALAS]      = 5.65e3;  // Fischetti
UL[ALP]       = 7.41e3;  // Fischetti
UL[GAP]       = 5.85e3;  // Fischetti
UL[GASB]      = 3.97e3;  // Fischetti
UL[INAS]      = 4.28e3;  // Fischetti
UL[INP]       = 5.13e3;  // Fischetti
UL[GAN]       = 6.56e3;  // Foutz, O'Leary, Shur, Eastman


// Lattice constants (m)
LATTCONST[GAAS]      = 565.35e-12;   // CODATA
LATTCONST[SILICON]   = 543.102e-12;  // CODATA
LATTCONST[GERMANIUM] = 564.613e-12;  // CODATA
LATTCONST[ALP]       = 545.10e-12;   // CODATA
LATTCONST[ALAS]      = 565.05e-12;   // CODATA
LATTCONST[ALSB]      = 613.55e-12;   // CODATA
LATTCONST[GAP]       = 545.12e-12;   // CODATA
LATTCONST[GASB]      = 609.59e-12;   // CODATA
LATTCONST[INP]       = 586.87e-12;   // CODATA
LATTCONST[INAS]      = 605.83e-12;   // CODATA
LATTCONST[INSB]      = 647.9e-12;    // CODATA
LATTCONST[GAN]       = 318.9e-12;    // a lattice constant

// electro-mechanical coupling constant
KAV[GAAS] = 0.0252;
KAV[GAN]  = 0.137;

// sp3s* Conduction Band calculated parameters for full band simulations
#include "materials/Silicon.h"
#include "materials/GaAs.h"
#include "materials/Germanium.h"
#include "materials/AlP.h"
#include "materials/AlAs.h"
#include "materials/AlSb.h"
#include "materials/GaP.h"
#include "materials/GaSb.h"
#include "materials/InP.h"
#include "materials/InAs.h"
#include "materials/InSb.h"
#include "materials/GaN.h"
