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

g_materials[SILICON].id   = SILICON;
g_materials[GERMANIUM].id = GERMANIUM;
g_materials[GAAS].id      = GAAS;
g_materials[INSB].id      = INSB;
g_materials[ALSB].id      = ALSB;
g_materials[ALXINXSB].id  = ALXINXSB;
g_materials[ALXIN1XSB].id = ALXIN1XSB;
g_materials[ALAS].id      = ALAS;
g_materials[ALP].id       = ALP;
g_materials[GAP].id       = GAP;
g_materials[GASB].id      = GASB;
g_materials[INAS].id      = INAS;
g_materials[INP].id       = INP;
g_materials[INXGA1XAS].id = INXGA1XAS;
g_materials[INXAL1XAS].id = INXAL1XAS;
g_materials[INXGAXXAS].id = INXGAXXAS;
g_materials[GAN].id       = GAN;


// Number of valleys
g_materials[SILICON].cb.num_valleys   = 1;  // X-valley
g_materials[GERMANIUM].cb.num_valleys = 1;  // G-valley
g_materials[GAAS].cb.num_valleys      = 1;  // G and L-valleys
g_materials[INSB].cb.num_valleys      = 1;  // G valley
g_materials[ALSB].cb.num_valleys      = 1;  // G-valley
g_materials[ALXINXSB].cb.num_valleys  = 1;  // G-valley
g_materials[ALXIN1XSB].cb.num_valleys = 1;  // G-valley
g_materials[ALAS].cb.num_valleys      = 1;  // G-valley
g_materials[ALP].cb.num_valleys       = 1;  // G-valley
g_materials[GAP].cb.num_valleys       = 1;  // G-valley
g_materials[GASB].cb.num_valleys      = 1;  // G-valley
g_materials[INAS].cb.num_valleys      = 1;  // G-valley
g_materials[INP].cb.num_valleys       = 1;  // G-valley
g_materials[INXGA1XAS].cb.num_valleys = 1;  // only G-valley
g_materials[INXAL1XAS].cb.num_valleys = 1;  // G-valley
g_materials[INXGAXXAS].cb.num_valleys = 1;  // only G-valley
g_materials[GAN].cb.num_valleys       = 2;  // G-1, M-L(U), G-3

g_materials[SILICON].vb.num_valleys   = 1;
g_materials[GERMANIUM].vb.num_valleys = 1;
g_materials[GAAS].vb.num_valleys      = 1;
g_materials[INSB].vb.num_valleys      = 1;
g_materials[ALSB].vb.num_valleys      = 1;
g_materials[ALXINXSB].vb.num_valleys  = 1;
g_materials[ALXIN1XSB].vb.num_valleys = 1;
g_materials[ALAS].vb.num_valleys      = 1;
g_materials[ALP].vb.num_valleys       = 1;
g_materials[GAP].vb.num_valleys       = 1;
g_materials[GASB].vb.num_valleys      = 1;
g_materials[INAS].vb.num_valleys      = 1;
g_materials[INP].vb.num_valleys       = 1;
g_materials[INXGA1XAS].vb.num_valleys = 1;
g_materials[INXAL1XAS].vb.num_valleys = 1;
g_materials[INXGAXXAS].vb.num_valleys = 1;
g_materials[GAN].vb.num_valleys       = 3;


g_materials[GAN].affinity = 3.18;


// Number of equivalent valleys
// Scattering from first index to second index
// Repeated indexes mean scattering between equivalent valleys
g_materials[GAAS].zscatter[1][1] = 0; // G -> G
g_materials[GAAS].zscatter[1][2] = 4; // G -> L
g_materials[GAAS].zscatter[1][3] = 3; // G -> X
g_materials[GAAS].zscatter[2][1] = 1; // L -> G
g_materials[GAAS].zscatter[2][2] = 3; // L -> L
g_materials[GAAS].zscatter[2][3] = 3; // L -> X
g_materials[GAAS].zscatter[3][1] = 1; // X -> G
g_materials[GAAS].zscatter[3][2] = 4; // X -> L
g_materials[GAAS].zscatter[3][3] = 2; // X -> X

g_materials[GAN].zscatter[1][1] = 0; // G1 -> G1
g_materials[GAN].zscatter[1][2] = 6; // G1 -> ML
g_materials[GAN].zscatter[1][3] = 1; // G1 -> G3
g_materials[GAN].zscatter[2][1] = 1; // ML -> G1
g_materials[GAN].zscatter[2][2] = 5; // ML -> ML
g_materials[GAN].zscatter[2][3] = 1; // ML -> G3
g_materials[GAN].zscatter[3][1] = 1; // G3 -> G1
g_materials[GAN].zscatter[3][2] = 6; // G3 -> ML
g_materials[GAN].zscatter[3][3] = 0; // G3 -> G3


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
EMIN[GAN][2]       = 1.34;    // L-M
// third valley
EMIN[GAAS][3]      = 0.48;   // X
EMIN[GAN][3]       = 2.14;    // G-2

// first valley
g_materials[SILICON].cb.emin[1]   = 0.0;    // Sellier, Fischetti, etc.
g_materials[GERMANIUM].cb.emin[1] = 0.173;
g_materials[GAAS].cb.emin[1]      = 0.0;    // Tomizawa
g_materials[INSB].cb.emin[1]      = 0.0;
g_materials[ALSB].cb.emin[1]      = 0.507;
g_materials[ALAS].cb.emin[1]      = 0.767;
g_materials[ALP].cb.emin[1]       = 1.237;
g_materials[GAP].cb.emin[1]       = 0.496;
g_materials[GASB].cb.emin[1]      = 0.0;
g_materials[INAS].cb.emin[1]      = 0.0;
g_materials[INP].cb.emin[1]       = 0.0;
g_materials[GAN].cb.emin[1]       = 0.0;    // G-1
// second valley
g_materials[GAAS].cb.emin[2]      = 0.323;  // L
g_materials[GAN].cb.emin[2]       = 1.34;   // L-M (U-3) - see Bulutay et al. Phys.Rev.B v62 p15454 (2000)
// third valley
g_materials[GAAS].cb.emin[3]      = 0.48;   // X
g_materials[GAN].cb.emin[3]       = 2.14;   // G-3 - see Bulutay et al. Phys.Rev.B v62 p15454 (2000)


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

g_materials[SILICON].cb.mstar[1]   = 0.32;    // see Sellier, Tomizawa, etc.
g_materials[GAAS].cb.mstar[1]      = 0.067;   // Gamma-valley -- see Tomizawa
g_materials[GAAS].cb.mstar[2]      = 0.350;   // L-valley     -- see Tomizawa
g_materials[GAAS].cb.mstar[3]      = 0.27;    // X-valley
g_materials[GERMANIUM].cb.mstar[1] = 0.12;    // Gamma valley -- see http://ecee.colorado.edu/~bart/book/effmass.htm#long
g_materials[INSB].cb.mstar[1]      = 0.0135;  // Gamma-valley -- see Ram-Mohan
g_materials[ALSB].cb.mstar[1]      = 0.14;    // Gamma-valley -- See Ram-Mohan
g_materials[ALAS].cb.mstar[1]      = 0.149;   // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
g_materials[ALP].cb.mstar[1]       = 0.22;    // Gamma-valley -- see Ram-Mohan
g_materials[GAP].cb.mstar[1]       = 0.13;    // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
g_materials[GASB].cb.mstar[1]      = 0.039;   // Gamma-valley -- see Ram-Mohan
g_materials[INAS].cb.mstar[1]      = 0.026;   // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
g_materials[INP].cb.mstar[1]       = 0.0795;  // Gamma-valley -- see Ram-Mohan
g_materials[GAN].cb.mstar[1]       = 0.2;     // G-1 -- Foutz, O'Leary, Shur, Eastman
g_materials[GAN].cb.mstar[2]       = 0.4;     // L-M -- Bhapkar & Shur
g_materials[GAN].cb.mstar[3]       = 0.6;     // G-2 -- Bhapkar & Shur


// non-parabolicity coefficients (1/eV)
alphaK[SILICON][1]   = 0.5;    // see Sellier, Tomizawa
alphaK[GERMANIUM][1] = 0.3;  // Gamma valley - Jacoboni Reggiani
alphaK[GAN][1]       = 0.189;      // G-1 -- Foutz, O'Leary, Shur, Eastman
alphaK[GAN][2]       = 0.065;      // L-M -- Bhapkar & Shur
alphaK[GAN][3]       = 0.029;      // G-2 -- Bhapkar & Shur

g_materials[SILICON].cb.alpha[1]   = 0.5;    // see Sellier, Tomizawa
g_materials[GERMANIUM].cb.alpha[1] = 0.3;  // Gamma valley - Jacoboni Reggiani
g_materials[GAN].cb.alpha[1]       = 0.189;      // G-1 -- Foutz, O'Leary, Shur, Eastman
g_materials[GAN].cb.alpha[2]       = 0.065;      // L-M -- Bhapkar & Shur
g_materials[GAN].cb.alpha[3]       = 0.029;      // G-2 -- Bhapkar & Shur


MSTAR_VB[GAN][0] = 1.4; // heavy hole
MSTAR_VB[GAN][1] = 0.3; // light hole
MSTAR_VB[GAN][2] = 0.6; // split-off

g_materials[GAN].vb.mstar[0] = 1.4; // heavy hole
g_materials[GAN].vb.mstar[1] = 0.3; // light hole
g_materials[GAN].vb.mstar[2] = 0.6; // split-off

DELTAE_VB[GAN][0] = 0.0;
DELTAE_VB[GAN][1] = 0.008; // difference between heavy hole and light hole bands
DELTAE_VB[GAN][2] = 0.04;  // difference between heavy hole and split-off bands

g_materials[GAN].vb.emin[0] = 0.0;
g_materials[GAN].vb.emin[1] = 0.008; // difference between heavy hole and light hole bands
g_materials[GAN].vb.emin[2] = 0.04;  // difference between heavy hole and split-off bands

g_materials[GAN].vb.alpha[0] = 0.0;
g_materials[GAN].vb.alpha[1] = 0.0;
g_materials[GAN].vb.alpha[2] = 0.0;


// precomputed constants
for(int m = 0; m < NOAMTIA; ++m) {
    for(int v = 1; v <= g_materials[m].cb.num_valleys; ++v) {
        g_materials[m].cb.smh[v] = sqrt(2. * g_materials[m].cb.mstar[v] * M * Q) / HBAR;
        g_materials[m].cb.hhm[v] = HBAR * HBAR / (2. * g_materials[m].cb.mstar[v] * M  * Q);
        g_materials[m].cb.hm[v]  = HBAR / (g_materials[m].cb.mstar[v] * M);
    }
}

// correction factor - multiplied by the calculated absorption constant
// used to correct for the fact that we are not taking into account all
// possible transistions which leads to a lower absorption rate
g_materials[SILICON].abs_correction   = 1.;
g_materials[GERMANIUM].abs_correction = 1.;
g_materials[GAAS].abs_correction      = 1.;
g_materials[INSB].abs_correction      = 1.;
g_materials[ALSB].abs_correction      = 1.;
g_materials[ALXINXSB].abs_correction  = 1.;
g_materials[ALXIN1XSB].abs_correction = 1.;
g_materials[ALAS].abs_correction      = 1.;
g_materials[ALP].abs_correction       = 1.;
g_materials[GAP].abs_correction       = 1.;
g_materials[GASB].abs_correction      = 1.;
g_materials[INAS].abs_correction      = 1.;
g_materials[INP].abs_correction       = 1.;
g_materials[INXGA1XAS].abs_correction = 1.;
g_materials[INXAL1XAS].abs_correction = 1.;
g_materials[INXGAXXAS].abs_correction = 1.;
g_materials[GAN].abs_correction       = 6.;


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

g_materials[SILICON].eps_static   = 11.68;      // see http://en.wikipedia.org/wiki/Relative_permittivity
g_materials[GERMANIUM].eps_static = 16.2;       // see http://www.ioffe.ru/SVA/NSM/Semicond/Ge/basic.html
g_materials[GAAS].eps_static      = 12.90;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
g_materials[INSB].eps_static      = 16.8;       // see http://www.ioffe.ru/SVA/NSM/Semicond/InSb/basic.html
g_materials[ALSB].eps_static      = 12.04;      // Fischetti conversations
g_materials[ALAS].eps_static      = 12.90-2.84; // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
g_materials[ALP].eps_static       = 9.80;       // Fischetti conversations
g_materials[GAP].eps_static       = 11.10;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
g_materials[GASB].eps_static      = 15.69;      // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
g_materials[INAS].eps_static      = 15.15;      // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
g_materials[INP].eps_static       = 12.50;      // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
g_materials[GAN].eps_static       = 9.7;        // E. Bellotti & F. Bertazzi

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

g_materials[GAAS].eps_hf = 10.89;       // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
g_materials[INSB].eps_hf = 15.68;       // Fischetti conversations
g_materials[ALSB].eps_hf = 9.88;        // Fiscehtti conversations
g_materials[ALAS].eps_hf = 10.89-2.73;  // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
g_materials[ALP].eps_hf  = 7.54;        // Fischetti conversations
g_materials[GAP].eps_hf  = 9.11;        // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
g_materials[GASB].eps_hf = 14.44;       // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
g_materials[INAS].eps_hf = 12.3;        // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
g_materials[INP].eps_hf  = 9.61;        // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
g_materials[GAN].eps_hf  = 5.28;        // E. Bellotti & F. Bertazzi


int m = 0,
    i = 0;
for(m = 0; m < NOAMTIA; m++) {
   for(i = 0; i < 6; i++){
       HWO[m][i] = 0.;
       DTK[m][i] = 0.;
       ZF[m][i]  = 0.;
       g_materials[m].hwo[i] = 0.;
       g_materials[m].dtk[i] = 0.;
       g_materials[m].zf[i]  = 0;
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

g_materials[SILICON].hwo[0]   = 0.0120;   // Sellier, Tomizawa
g_materials[SILICON].hwo[1]   = 0.0185;   // Sellier, Tomizawa
g_materials[SILICON].hwo[2]   = 0.0190;   // Sellier, Tomizawa
g_materials[SILICON].hwo[3]   = 0.0474;   // Sellier, Tomizawa
g_materials[SILICON].hwo[4]   = 0.0612;   // Sellier, Tomizawa
g_materials[SILICON].hwo[5]   = 0.0590;   // Sellier, Tomizawa
g_materials[GERMANIUM].hwo[0] = 0.03704;  // Fischetti
g_materials[GAAS].hwo[0]      = 0.03536;  // Fischetti
g_materials[INSB].hwo[0]      = 0.02404;  // Fischetti
g_materials[ALSB].hwo[0]      = 0.0360;   // Fischetti
g_materials[ALAS].hwo[0]      = 0.05009;  // Fischetti
g_materials[ALP].hwo[0]       = 0.06211;  // Fischetti
g_materials[GAP].hwo[0]       = 0.04523;  // Fischetti
g_materials[GASB].hwo[0]      = 0.02529;  // Fischetti
g_materials[INAS].hwo[0]      = 0.03008;  // Fischetti
g_materials[INP].hwo[0]       = 0.04240;  // Fischetti
g_materials[GAN].hwo[0]       = 0.09212;  // LO -- E. Bellotti & F. Bertazzi
g_materials[GAN].hwo[1]       = 0.06955;  // TO -- E. Bellotti & F. Bertazzi

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

g_materials[SILICON].dtk[0]   = 0.05e11;  // Jacoboni Reggiani
g_materials[SILICON].dtk[1]   = 0.08e11;  // Jacoboni Reggiani
g_materials[SILICON].dtk[2]   = 0.03e11;  // Jacoboni Reggiani
g_materials[SILICON].dtk[3]   = 0.20e11;  // Jacoboni Reggiani
g_materials[SILICON].dtk[4]   = 1.14e11;  // Jacoboni Reggiani
g_materials[SILICON].dtk[5]   = 0.20e11;  // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[0] = 0.0;      // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[1] = 0.079e11; // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[2] = 0.0;      // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[3] = 0.0;      // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[4] = 0.95e11;  // Jacoboni Reggiani
g_materials[GERMANIUM].dtk[5] = 0.0;      // Jacoboni Reggiani
g_materials[GAAS].dtk[0]      = 1.11e11;  // Sellier, Tomizawa
g_materials[INSB].dtk[0]      = 0.47e11;  // see ???
g_materials[ALSB].dtk[0]      = 0.55e11;  // see ???
g_materials[ALAS].dtk[0]      = 3.0e11;   // see ???
g_materials[ALP].dtk[0]       = 0.95e11;  // see ???
g_materials[GAP].dtk[0]       = 5.33e11;  // see ???
g_materials[GASB].dtk[0]      = 0.94e11;  // see ???
g_materials[INAS].dtk[0]      = 3.59e11;  // see ???
g_materials[INP].dtk[0]       = 2.46e11;  // see ???
g_materials[GAN].dtk[0]       = 1.0e11;   // E. Bellotti & F. Bertazzi
g_materials[GAN].dtk[1]       = 1.0e11;   // E. Bellotti & F. Bertazzi

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
ZF[GAN][1]       = 1.;  // guess for correct value

g_materials[SILICON].zf[0]   = 1.;  // Sellier
g_materials[SILICON].zf[1]   = 1.;  // Sellier
g_materials[SILICON].zf[2]   = 4.;  // Sellier
g_materials[SILICON].zf[3]   = 4.;  // Sellier
g_materials[SILICON].zf[4]   = 1.;  // Sellier
g_materials[SILICON].zf[5]   = 4.;  // Sellier
g_materials[GERMANIUM].zf[0] = 1.;  // see ???
g_materials[GAAS].zf[0]      = 1.;  // Sellier
g_materials[INSB].zf[0]      = 1.;  // see ???
g_materials[ALSB].zf[0]      = 1.;  // see ???
g_materials[ALAS].zf[0]      = 1.;  // see ???
g_materials[ALP].zf[0]       = 1.;  // see ???
g_materials[GAP].zf[0]       = 1.;  // see ???
g_materials[GASB].zf[0]      = 1.;  // see ???
g_materials[INAS].zf[0]      = 1.;  // see ???
g_materials[INP].zf[0]       = 1.;  // see ???
g_materials[GAN].zf[0]       = 1.;  // guess for correct value
g_materials[GAN].zf[1]       = 1.;  // guess for correct value

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

g_materials[SILICON].rho   = 2.33e3;  // Fischetti conversations
g_materials[GERMANIUM].rho = 5.32e3;  // Fischetti conversations
g_materials[GAAS].rho      = 5.36e3;  // Fischetti conversations
g_materials[INSB].rho      = 5.78e3;  // Fischetti conversations
g_materials[ALSB].rho      = 4.26e3;  // Fischetti conversations
g_materials[ALAS].rho      = 3.76e3;  // Fischetti conversations
g_materials[ALP].rho       = 2.40e3;  // Fischetti conversations
g_materials[GAP].rho       = 4.14e3;  // Fischetti conversations
g_materials[GASB].rho      = 5.61e3;  // Fischetti conversations
g_materials[INAS].rho      = 5.67e3;  // Fischetti conversations
g_materials[INP].rho       = 4.81e3;  // Fischetti conversations
g_materials[GAN].rho       = 6.087e3; // E. Bellotti & F. Bertazzi

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

g_materials[SILICON].da   = 9.  * Q;  // Fischetti -- Jacoboni Reggiani
g_materials[GERMANIUM].da = 9.  * Q;  // Fischetti -- Jacoboni Reggiani
g_materials[GAAS].da      = 7.  * Q;  // Fischetti - Gamma valley
g_materials[INSB].da      = 7.  * Q;  // Fischetti
g_materials[ALSB].da      = 4.6 * Q;  // Fischetti
g_materials[ALAS].da      = 9.3 * Q;  // Fischetti
g_materials[ALP].da       = 9.3 * Q;  // Fischetti
g_materials[GAP].da       = 7.4 * Q;  // Fischetti
g_materials[GASB].da      = 9.  * Q;  // Fischetti
g_materials[INAS].da      = 8.2 * Q;  // Fischetti
g_materials[INP].da       = 6.2 * Q;  // Fischetti
g_materials[GAN].da       = 8.3 * Q;  // E. Bellotti & F. Bertazzi

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

g_materials[SILICON].ul   = 9.18e3;  // Fischetti
g_materials[GERMANIUM].ul = 5.4e3;   // Fischetti
g_materials[GAAS].ul      = 5.24e3;  // Fischetti
g_materials[INSB].ul      = 3.41e3;  // Fischetti
g_materials[ALSB].ul      = 4.25e3;  // Fischetti
g_materials[ALAS].ul      = 5.65e3;  // Fischetti
g_materials[ALP].ul       = 7.41e3;  // Fischetti
g_materials[GAP].ul       = 5.85e3;  // Fischetti
g_materials[GASB].ul      = 3.97e3;  // Fischetti
g_materials[INAS].ul      = 4.28e3;  // Fischetti
g_materials[INP].ul       = 5.13e3;  // Fischetti
g_materials[GAN].ul       = 6.56e3;  // Foutz, O'Leary, Shur, Eastman


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

g_materials[GAAS].lattice_const      = 565.35e-12;   // CODATA
g_materials[SILICON].lattice_const   = 543.102e-12;  // CODATA
g_materials[GERMANIUM].lattice_const = 564.613e-12;  // CODATA
g_materials[ALP].lattice_const       = 545.10e-12;   // CODATA
g_materials[ALAS].lattice_const      = 565.05e-12;   // CODATA
g_materials[ALSB].lattice_const      = 613.55e-12;   // CODATA
g_materials[GAP].lattice_const       = 545.12e-12;   // CODATA
g_materials[GASB].lattice_const      = 609.59e-12;   // CODATA
g_materials[INP].lattice_const       = 586.87e-12;   // CODATA
g_materials[INAS].lattice_const      = 605.83e-12;   // CODATA
g_materials[INSB].lattice_const      = 647.9e-12;    // CODATA
g_materials[GAN].lattice_const       = 318.9e-12;    // a lattice constant

// electro-mechanical coupling constant
KAV[GAAS] = 0.0252;
KAV[GAN]  = 0.137;

g_materials[GAAS].kav = 0.0252;
g_materials[GAN].kav  = 0.137;

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

// III-V Semiconductor materials energy gap (depending on the lattice temperature)
    printf("\n");
    EG[SILICON]=1.21-3.333e-4*g_config->lattice_temp;
    printf("EG[SILICON]      = %g\n",EG[SILICON]);
    EG[GERMANIUM]=0.747-3.587e-4*g_config->lattice_temp;
    printf("EG[GERMANIUM]    = %g\n",EG[GERMANIUM]);
    EG[GAAS]=1.54-4.036e-4*g_config->lattice_temp;
    printf("EG[GAAS]         = %g\n",EG[GAAS]);
    EG[INSB]=0.2446-2.153e-4*g_config->lattice_temp;
    printf("EG[INSB]         = %g\n",EG[INSB]);
    EG[ALSB]=1.696-2.20e-4*g_config->lattice_temp;
    printf("EG[ALSB]         = %g\n",EG[ALSB]);
    EG[ALAS]=2.314-3.0e-4*g_config->lattice_temp;
    printf("EG[ALAS]         = %g\n",EG[ALAS]);
    EG[ALP]=2.51-3.333e-4*g_config->lattice_temp;
    printf("EG[ALP]         = %g\n",EG[ALP]);
    EG[GAP]=2.35-2.667e-4*g_config->lattice_temp;
    printf("EG[GAP]         = %g\n",EG[GAP]);
    EG[GASB]=0.81-3.667e-4*g_config->lattice_temp;
    printf("EG[GASB]         = %g\n",EG[GASB]);
    EG[INAS]=0.434-2.601e-4*g_config->lattice_temp;
    printf("EG[INAS]         = %g\n",EG[INAS]);
    EG[INP]=1.445-3.296e-4*g_config->lattice_temp;
    printf("EG[INP]         = %g\n",EG[INP]);
    EG[GAN]=3.47 - 7.7e-4 * g_config->lattice_temp * g_config->lattice_temp / (g_config->lattice_temp + 600.);
    printf("EG[GAN]         = %g\n", EG[GAN]);
    printf("\n");

    g_materials[SILICON].Eg=1.21-3.333e-4*g_config->lattice_temp;
    g_materials[GERMANIUM].Eg=0.747-3.587e-4*g_config->lattice_temp;
    g_materials[GAAS].Eg=1.54-4.036e-4*g_config->lattice_temp;
    g_materials[INSB].Eg=0.2446-2.153e-4*g_config->lattice_temp;
    g_materials[ALSB].Eg=1.696-2.20e-4*g_config->lattice_temp;
    g_materials[ALAS].Eg=2.314-3.0e-4*g_config->lattice_temp;
    g_materials[ALP].Eg=2.51-3.333e-4*g_config->lattice_temp;
    g_materials[GAP].Eg=2.35-2.667e-4*g_config->lattice_temp;
    g_materials[GASB].Eg=0.81-3.667e-4*g_config->lattice_temp;
    g_materials[INAS].Eg=0.434-2.601e-4*g_config->lattice_temp;
    g_materials[INP].Eg=1.445-3.296e-4*g_config->lattice_temp;
    g_materials[GAN].Eg=3.47 - 7.7e-4 * g_config->lattice_temp * g_config->lattice_temp / (g_config->lattice_temp + 600.);

    if(g_config->conduction_band == KANE ||
       g_config->conduction_band == PARABOLIC ||
       g_config->conduction_band == FULL) {
        // USED WHATEVER IS THE CONDUCTION BAND FOR THE INITIAL PSUEDO WAVE VECTOR
        // OF THE PSEUDO PARTICLES
        // all the following non-parabolicity coefficients depend on lattice temperature
        // non-parabolicity coefficient for GaAs in the GAMMA-valley
        alphaK[GAAS][1]=pow(1.-MSTAR[GAAS][1],2.)/(EG[GAAS]+EMIN[GAAS][1]);//expected value = 0.611
        printf("alphaK_gamma[GaAs] = %g\n",alphaK[GAAS][1]);
        // non-parabolicity coefficient for GaAs in the L-valley
        alphaK[GAAS][2]=pow(1.-MSTAR[GAAS][2],2.)/(EG[GAAS]+EMIN[GAAS][2]);//expected value = 0.242;
        printf("alphaK_L[GaAs]     = %g\n",alphaK[GAAS][2]);
        // non-parabolicity coefficient for InSb in the GAMMA-valley
        alphaK[INSB][1]=pow(1.-MSTAR[INSB][1],2.)/(EG[INSB]+EMIN[INSB][1]);//5.59;
        printf("alphaK_gamma[InSb] = %g\n",alphaK[INSB][1]);
        // non-parabolicity coefficient for AlSb in the GAMMA-valley
        alphaK[ALSB][1]=pow(1.-MSTAR[ALSB][1],2.)/(EG[ALSB]+EMIN[ALSB][1]);//0.321;
        printf("alphaK_gamma[AlSb] = %g\n",alphaK[ALSB][1]);
        // non-parabolicity coefficient for AlAs in the GAMMA-valley
        alphaK[ALAS][1]=pow(1.-MSTAR[ALAS][1],2.)/(EG[ALAS]+EMIN[ALAS][1]);
        printf("alphaK_gamma[AlAs] = %g\n",alphaK[ALAS][1]);
        // non-parabolicity coefficient for AlP in the GAMMA-valley
        alphaK[ALP][1]=pow(1.-MSTAR[ALP][1],2.)/(EG[ALP]+EMIN[ALP][1]);
        printf("alphaK_gamma[AlP] = %g\n",alphaK[ALP][1]);
        // non-parabolicity coefficient for GaP in the GAMMA-valley
        alphaK[GAP][1]=pow(1.-MSTAR[GAP][1],2.)/(EG[GAP]+EMIN[GAP][1]);
        printf("alphaK_gamma[GaP] = %g\n",alphaK[GAP][1]);
        // non-parabolicity coefficient for GaSb in the GAMMA-valley
        alphaK[GASB][1]=pow(1.-MSTAR[GASB][1],2.)/(EG[GASB]+EMIN[GASB][1]);
        printf("alphaK_gamma[GaSb] = %g\n",alphaK[GASB][1]);
        // non-parabolicity coefficient for InAs in the GAMMA-valley
        alphaK[INAS][1]=pow(1.-MSTAR[INAS][1],2.)/(EG[INAS]+EMIN[INAS][1]);
        printf("alphaK_gamma[InAs] = %g\n",alphaK[INAS][1]);
        // non-parabolicity coefficient for InP in the GAMMA-valley
        alphaK[INP][1]=pow(1.-MSTAR[INP][1],2.)/(EG[INP]+EMIN[INP][1]);
        printf("alphaK_gamma[InP] = %g\n",alphaK[INP][1]);
        alphaK[GAN][1] = pow(1. - MSTAR[GAN][1], 2.) / (EG[GAN] + EMIN[GAN][1]); // expected value = 0.189
        printf("alphaK_gamma1[GAN] = %g\n", alphaK[GAN][1]);


        g_materials[GAAS].cb.alpha[1] = pow(1.-g_materials[GAAS].cb.mstar[1],2.)/(g_materials[GAAS].Eg+g_materials[GAAS].cb.emin[1]);//expected value = 0.611
        g_materials[GAAS].cb.alpha[2] = pow(1.-g_materials[GAAS].cb.mstar[2],2.)/(g_materials[GAAS].Eg+g_materials[GAAS].cb.emin[2]);//expected value = 0.242;
        g_materials[INSB].cb.alpha[1] = pow(1.-g_materials[INSB].cb.mstar[1],2.)/(g_materials[INSB].Eg+g_materials[INSB].cb.emin[1]);//5.59;
        g_materials[ALSB].cb.alpha[1] = pow(1.-g_materials[ALSB].cb.mstar[1],2.)/(g_materials[ALSB].Eg+g_materials[ALSB].cb.emin[1]);//0.321;
        g_materials[ALAS].cb.alpha[1] = pow(1.-g_materials[ALAS].cb.mstar[1],2.)/(g_materials[ALAS].Eg+g_materials[ALAS].cb.emin[1]);
        g_materials[ALP].cb.alpha[1]  = pow(1.-g_materials[ALP].cb.mstar[1],2.)/(g_materials[ALP].Eg+g_materials[ALP].cb.emin[1]);
        g_materials[GAP].cb.alpha[1]  = pow(1.-g_materials[GAP].cb.mstar[1],2.)/(g_materials[GAP].Eg+g_materials[GAP].cb.emin[1]);
        g_materials[GASB].cb.alpha[1] = pow(1.-g_materials[GASB].cb.mstar[1],2.)/(g_materials[GASB].Eg+g_materials[GASB].cb.emin[1]);
        g_materials[INAS].cb.alpha[1] = pow(1.-g_materials[INAS].cb.mstar[1],2.)/(g_materials[INAS].Eg+g_materials[INAS].cb.emin[1]);
        g_materials[INP].cb.alpha[1]  = pow(1.-g_materials[INP].cb.mstar[1],2.)/(g_materials[INP].Eg+g_materials[INP].cb.emin[1]);
        g_materials[GAN].cb.alpha[1]  = pow(1. - g_materials[GAN].cb.mstar[1], 2.) / (g_materials[GAN].Eg + g_materials[GAN].cb.emin[1]); // expected value = 0.189

    }

    // Semiconductor compounds
    // ***
    // Relative dielectric constant for semiconductor compounds
    EPSR[ALXINXSB]=XVAL[ALXINXSB]*EPSR[ALSB]+XVAL[ALXINXSB]*EPSR[INSB];
    EPSR[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPSR[ALSB]+(1.-XVAL[ALXIN1XSB])*EPSR[INSB];
    EPSR[INXGA1XAS]=XVAL[INXGA1XAS]*EPSR[INAS]+(1.-XVAL[INXGA1XAS])*EPSR[GAAS];
    EPSR[INXAL1XAS]=XVAL[INXAL1XAS]*EPSR[INAS]+(1.-XVAL[INXAL1XAS])*EPSR[ALAS];
    EPSR[INXGAXXAS]=XVAL[INXGAXXAS]*EPSR[INAS]+XVAL[INXGAXXAS]*EPSR[GAAS];
    // semiconductor compounds high frequency dieletric constant
    EPF[ALXINXSB]=XVAL[ALXINXSB]*(EPF[ALSB]+EPF[INSB]);
    EPF[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPF[ALSB]+(1.-XVAL[ALXIN1XSB])*EPF[INSB];
    EPF[INXGA1XAS]=XVAL[INXGA1XAS]*EPF[INAS]+(1.-XVAL[INXGA1XAS])*EPF[GAAS];
    EPF[INXAL1XAS]=XVAL[INXAL1XAS]*EPF[INAS]+(1.-XVAL[INXAL1XAS])*EPF[ALAS];
    EPF[INXGAXXAS]=XVAL[INXGAXXAS]*EPF[INAS]+XVAL[INXGAXXAS]*EPF[GAAS];
    // semiconductor compounds optical phonon scattering energy (eV)
    HWO[ALXINXSB][0]=XVAL[ALXINXSB]*(HWO[ALSB][0]+HWO[INSB][0]);
    HWO[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*HWO[ALSB][0]+(1.-XVAL[ALXIN1XSB])*HWO[INSB][0];
    HWO[INXGA1XAS][0]=XVAL[INXGA1XAS]*HWO[INAS][0]+(1.-XVAL[INXGA1XAS])*HWO[GAAS][0];
    HWO[INXAL1XAS][0]=XVAL[INXAL1XAS]*HWO[INAS][0]+(1.-XVAL[INXAL1XAS])*HWO[ALAS][0];
    HWO[INXGAXXAS][0]=XVAL[INXGAXXAS]*HWO[INAS][0]+XVAL[INXGAXXAS]*HWO[GAAS][0];
    // semiconductor compounds optical coupling constants (eV/m)
    DTK[ALXINXSB][0]=XVAL[ALXINXSB]*(DTK[ALSB][0]+DTK[INSB][0]);
    DTK[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*DTK[ALSB][0]+(1.-XVAL[ALXIN1XSB])*DTK[INSB][0];
    DTK[INXGA1XAS][0]=XVAL[INXGA1XAS]*DTK[INAS][0]+(1.-XVAL[INXGA1XAS])*DTK[GAAS][0];
    DTK[INXAL1XAS][0]=XVAL[INXAL1XAS]*DTK[INAS][0]+(1.-XVAL[INXAL1XAS])*DTK[ALAS][0];
    DTK[INXGAXXAS][0]=XVAL[INXGAXXAS]*DTK[INAS][0]+XVAL[INXGAXXAS]*DTK[GAAS][0];
    // semiconductor compounds optical phonon Z-factor
    ZF[ALXINXSB][0]=XVAL[ALXINXSB]*(ZF[ALSB][0]+ZF[INSB][0]);
    ZF[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*ZF[ALSB][0]+(1.-XVAL[ALXIN1XSB])*ZF[INSB][0];
    ZF[INXGA1XAS][0]=XVAL[INXGA1XAS]*ZF[INAS][0]+(1.-XVAL[INXGA1XAS])*ZF[GAAS][0];
    ZF[INXAL1XAS][0]=XVAL[INXAL1XAS]*ZF[INAS][0]+(1.-XVAL[INXAL1XAS])*ZF[ALAS][0];
    ZF[INXGAXXAS][0]=XVAL[INXGAXXAS]*ZF[INAS][0]+XVAL[INXGAXXAS]*ZF[GAAS][0];
    // semiconductor compounds Crystal Density (Kg/m^3)
    RHO[ALXINXSB]=XVAL[ALXINXSB]*(RHO[ALSB]+RHO[INSB]);
    RHO[ALXIN1XSB]=XVAL[ALXIN1XSB]*RHO[ALSB]+(1.-XVAL[ALXIN1XSB])*RHO[INSB];
    RHO[INXGA1XAS]=XVAL[INXGA1XAS]*RHO[INAS]+(1.-XVAL[INXGA1XAS])*RHO[GAAS];
    RHO[INXAL1XAS]=XVAL[INXAL1XAS]*RHO[INAS]+(1.-XVAL[INXAL1XAS])*RHO[ALAS];
    RHO[INXGAXXAS]=XVAL[INXGAXXAS]*RHO[INAS]+XVAL[INXGAXXAS]*RHO[GAAS];
    // semiconductor compounds acoustic deformation potential (Joule)
    DA[ALXINXSB]=XVAL[ALXINXSB]*(DA[ALSB]+DA[INSB]);
    DA[ALXIN1XSB]=XVAL[ALXIN1XSB]*DA[ALSB]+(1.-XVAL[ALXIN1XSB])*DA[INSB];
    DA[INXGA1XAS]=XVAL[INXGA1XAS]*DA[INAS]+(1.-XVAL[INXGA1XAS])*DA[GAAS];
    DA[INXAL1XAS]=XVAL[INXAL1XAS]*DA[INAS]+(1.-XVAL[INXAL1XAS])*DA[ALAS];
    DA[INXGAXXAS]=XVAL[INXGAXXAS]*DA[INAS]+XVAL[INXGAXXAS]*DA[GAAS];
    // semiconductor compounds longitudinal sound velocity (m/sec)
    UL[ALXINXSB]=XVAL[ALXINXSB]*(UL[ALSB]+UL[INSB]);
    UL[ALXIN1XSB]=XVAL[ALXIN1XSB]*UL[ALSB]+(1.-XVAL[ALXIN1XSB])*UL[INSB];
    UL[INXGA1XAS]=XVAL[INXGA1XAS]*UL[INAS]+(1.-XVAL[INXGA1XAS])*UL[GAAS];
    UL[INXAL1XAS]=XVAL[INXAL1XAS]*UL[INAS]+(1.-XVAL[INXAL1XAS])*UL[ALAS];
    UL[INXGAXXAS]=XVAL[INXGAXXAS]*UL[INAS]+XVAL[INXGAXXAS]*UL[GAAS];
    // semiconductor compounds energy gap
    EG[ALXINXSB]=XVAL[ALXINXSB]*(EG[ALSB]+EG[INSB]);
    EG[ALXIN1XSB]=XVAL[ALXIN1XSB]*EG[ALSB]+(1.-XVAL[ALXIN1XSB])*EG[INSB];
    EG[INXGA1XAS]=XVAL[INXGA1XAS]*EG[INAS]+(1.-XVAL[INXGA1XAS])*EG[GAAS];
    EG[INXAL1XAS]=XVAL[INXAL1XAS]*EG[INAS]+(1.-XVAL[INXAL1XAS])*EG[ALAS];
    EG[INXGAXXAS]=XVAL[INXGAXXAS]*EG[INAS]+XVAL[INXGAXXAS]*EG[GAAS];
    // semiconductor compounds energy minimum of GAMMA-valley
    EMIN[ALXINXSB][1]=XVAL[ALXINXSB]*(EMIN[ALSB][1]+EMIN[INSB][1]);
    EMIN[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*EMIN[ALSB][1]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][1];
    EMIN[INXGA1XAS][1]=XVAL[INXGA1XAS]*EMIN[INAS][1]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][1];
    EMIN[INXAL1XAS][1]=XVAL[INXAL1XAS]*EMIN[INAS][1]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][1];
    EMIN[INXGAXXAS][1]=XVAL[INXGAXXAS]*EMIN[INAS][1]+XVAL[INXGAXXAS]*EMIN[GAAS][1];
    // semiconductor compounds energy minimum 0f L-valley
    EMIN[ALXINXSB][2]=XVAL[ALXINXSB]*(EMIN[ALSB][2]+EMIN[INSB][2]);
    EMIN[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*EMIN[ALSB][2]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][2];
    EMIN[INXGA1XAS][2]=XVAL[INXGA1XAS]*EMIN[INAS][2]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][2];
    EMIN[INXAL1XAS][2]=XVAL[INXAL1XAS]*EMIN[INAS][2]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][2];
    EMIN[INXGAXXAS][2]=XVAL[INXGAXXAS]*EMIN[INAS][2]+XVAL[INXGAXXAS]*EMIN[GAAS][2];
    // GAMMA-valley effective mass
    MSTAR[ALXINXSB][1]=XVAL[ALXINXSB]*(MSTAR[ALSB][1]+MSTAR[INSB][1]);
    MSTAR[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*MSTAR[ALSB][1]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][1];
    MSTAR[INXGA1XAS][1]=XVAL[INXGA1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][1];
    MSTAR[INXAL1XAS][1]=XVAL[INXAL1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][1];
    MSTAR[INXGAXXAS][1]=XVAL[INXGAXXAS]*MSTAR[INAS][1]+XVAL[INXGAXXAS]*MSTAR[GAAS][1];
    // L-valley effective mass
    MSTAR[ALXINXSB][2]=XVAL[ALXINXSB]*(MSTAR[ALSB][2]+MSTAR[INSB][2]);
    MSTAR[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*MSTAR[ALSB][2]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][2];
    MSTAR[INXGA1XAS][2]=XVAL[INXGA1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][2];
    MSTAR[INXAL1XAS][2]=XVAL[INXAL1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][2];
    MSTAR[INXGAXXAS][2]=XVAL[INXGAXXAS]*MSTAR[INAS][2]+XVAL[INXGAXXAS]*MSTAR[GAAS][2];
    // non-parabolicity coefficient for semiconductor compounds in the GAMMA-valley
    alphaK[ALXINXSB][1]=XVAL[ALXINXSB]*(alphaK[ALSB][1]+alphaK[INSB][1]);
    alphaK[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*alphaK[ALSB][1]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][1];
    alphaK[INXGA1XAS][1]=XVAL[INXGA1XAS]*alphaK[INAS][1]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][1];
    alphaK[INXAL1XAS][1]=XVAL[INXAL1XAS]*alphaK[INAS][1]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][1];
    alphaK[INXGAXXAS][1]=XVAL[INXGAXXAS]*alphaK[INAS][1]+XVAL[INXGAXXAS]*alphaK[GAAS][1];
    // non-parabolicity coefficient for Al_x In_(1-x) Sb in the L-valley
    alphaK[ALXINXSB][2]=XVAL[ALXINXSB]*(alphaK[ALSB][2]+alphaK[INSB][2]);
    alphaK[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*alphaK[ALSB][2]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][2];
    alphaK[INXGA1XAS][2]=XVAL[INXGA1XAS]*alphaK[INAS][2]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][2];
    alphaK[INXAL1XAS][2]=XVAL[INXAL1XAS]*alphaK[INAS][2]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][2];
    alphaK[INXGAXXAS][2]=XVAL[INXGAXXAS]*alphaK[INAS][2]+XVAL[INXGAXXAS]*alphaK[GAAS][2];
