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
g_materials[SILICON].cb.alpha[1]   = 0.5;    // see Sellier, Tomizawa
g_materials[GERMANIUM].cb.alpha[1] = 0.3;    // Gamma valley - Jacoboni Reggiani
g_materials[GAN].cb.alpha[1]       = 0.189;  // G-1 -- Foutz, O'Leary, Shur, Eastman
g_materials[GAN].cb.alpha[2]       = 0.065;  // L-M -- Bhapkar & Shur
g_materials[GAN].cb.alpha[3]       = 0.029;  // G-2 -- Bhapkar & Shur


g_materials[GAN].vb.mstar[0] = 1.4; // heavy hole
g_materials[GAN].vb.mstar[1] = 0.3; // light hole
g_materials[GAN].vb.mstar[2] = 0.6; // split-off

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
       g_materials[m].hwo[i] = 0.;
       g_materials[m].dtk[i] = 0.;
       g_materials[m].zf[i]  = 0;
   }
}

// Optical phonon scattering energy (eV)
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
    g_materials[SILICON].Eg=1.21-3.333e-4*g_config->lattice_temp;
    printf("EG[SILICON]      = %g\n",g_materials[SILICON].Eg);
    g_materials[GERMANIUM].Eg=0.747-3.587e-4*g_config->lattice_temp;
    printf("EG[GERMANIUM]    = %g\n",g_materials[GERMANIUM].Eg);
    g_materials[GAAS].Eg=1.54-4.036e-4*g_config->lattice_temp;
    printf("EG[GAAS]         = %g\n",g_materials[GAAS].Eg);
    g_materials[INSB].Eg=0.2446-2.153e-4*g_config->lattice_temp;
    printf("EG[INSB]         = %g\n",g_materials[INSB].Eg);
    g_materials[ALSB].Eg=1.696-2.20e-4*g_config->lattice_temp;
    printf("EG[ALSB]         = %g\n",g_materials[ALSB].Eg);
    g_materials[ALAS].Eg=2.314-3.0e-4*g_config->lattice_temp;
    printf("EG[ALAS]         = %g\n",g_materials[ALAS].Eg);
    g_materials[ALP].Eg=2.51-3.333e-4*g_config->lattice_temp;
    printf("EG[ALP]         = %g\n",g_materials[ALP].Eg);
    g_materials[GAP].Eg=2.35-2.667e-4*g_config->lattice_temp;
    printf("EG[GAP]         = %g\n",g_materials[GAP].Eg);
    g_materials[GASB].Eg=0.81-3.667e-4*g_config->lattice_temp;
    printf("EG[GASB]         = %g\n",g_materials[GASB].Eg);
    g_materials[INAS].Eg=0.434-2.601e-4*g_config->lattice_temp;
    printf("EG[INAS]         = %g\n",g_materials[INAS].Eg);
    g_materials[INP].Eg=1.445-3.296e-4*g_config->lattice_temp;
    printf("EG[INP]         = %g\n",g_materials[INP].Eg);
    g_materials[GAN].Eg=3.47 - 7.7e-4 * g_config->lattice_temp * g_config->lattice_temp / (g_config->lattice_temp + 600.);
    printf("EG[GAN]         = %g\n", g_materials[GAN].Eg);
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
        // USED WHATEVER IS THcE CONDUCTION BAND FOR THE INITIAL PSUEDO WAVE VECTOR
        // OF THE PSEUDO PARTICLES
        // all the following non-parabolicity coefficients depend on lattice temperature
        // non-parabolicity coefficient for GaAs in the GAMMA-valley
        g_materials[GAAS].cb.alpha[1]=pow(1.-g_materials[GAAS].cb.mstar[1],2.)/(g_materials[GAAS].Eg+g_materials[GAAS].cb.emin[1]);//expected value = 0.611
        printf("alphaK_gamma[GaAs] = %g\n",g_materials[GAAS].cb.alpha[1]);
        // non-parabolicity coefficient for GaAs in the L-valley
        g_materials[GAAS].cb.alpha[2]=pow(1.-g_materials[GAAS].cb.mstar[2],2.)/(g_materials[GAAS].Eg+g_materials[GAAS].cb.emin[2]);//expected value = 0.242;
        printf("alphaK_L[GaAs]     = %g\n",g_materials[GAAS].cb.alpha[2]);
        // non-parabolicity coefficient for InSb in the GAMMA-valley
        g_materials[INSB].cb.alpha[1]=pow(1.-g_materials[INSB].cb.mstar[1],2.)/(g_materials[INSB].Eg+g_materials[INSB].cb.emin[1]);//5.59;
        printf("alphaK_gamma[InSb] = %g\n",g_materials[INSB].cb.alpha[1]);
        // non-parabolicity coefficient for AlSb in the GAMMA-valley
        g_materials[ALSB].cb.alpha[1]=pow(1.-g_materials[ALSB].cb.mstar[1],2.)/(g_materials[ALSB].Eg+g_materials[ALSB].cb.emin[1]);//0.321;
        printf("alphaK_gamma[AlSb] = %g\n",g_materials[ALSB].cb.alpha[1]);
        // non-parabolicity coefficient for AlAs in the GAMMA-valley
        g_materials[ALAS].cb.alpha[1]=pow(1.-g_materials[ALAS].cb.mstar[1],2.)/(g_materials[ALAS].Eg+g_materials[ALAS].cb.emin[1]);
        printf("alphaK_gamma[AlAs] = %g\n",g_materials[ALAS].cb.alpha[1]);
        // non-parabolicity coefficient for AlP in the GAMMA-valley
        g_materials[ALP].cb.alpha[1]=pow(1.-g_materials[ALP].cb.mstar[1],2.)/(g_materials[ALP].Eg+g_materials[ALP].cb.emin[1]);
        printf("alphaK_gamma[AlP] = %g\n",g_materials[ALP].cb.alpha[1]);
        // non-parabolicity coefficient for GaP in the GAMMA-valley
        g_materials[GAP].cb.alpha[1]=pow(1.-g_materials[GAP].cb.mstar[1],2.)/(g_materials[GAP].Eg+g_materials[GAP].cb.emin[1]);
        printf("alphaK_gamma[GaP] = %g\n",g_materials[GAP].cb.alpha[1]);
        // non-parabolicity coefficient for GaSb in the GAMMA-valley
        g_materials[GASB].cb.alpha[1]=pow(1.-g_materials[GASB].cb.mstar[1],2.)/(g_materials[GASB].Eg+g_materials[GASB].cb.emin[1]);
        printf("alphaK_gamma[GaSb] = %g\n",g_materials[GASB].cb.alpha[1]);
        // non-parabolicity coefficient for InAs in the GAMMA-valley
        g_materials[INAS].cb.alpha[1]=pow(1.-g_materials[INAS].cb.mstar[1],2.)/(g_materials[INAS].Eg+g_materials[INAS].cb.emin[1]);
        printf("alphaK_gamma[InAs] = %g\n",g_materials[INAS].cb.alpha[1]);
        // non-parabolicity coefficient for InP in the GAMMA-valley
        g_materials[INP].cb.alpha[1]=pow(1.-g_materials[INP].cb.mstar[1],2.)/(g_materials[INP].Eg+g_materials[INP].cb.emin[1]);
        printf("alphaK_gamma[InP] = %g\n",g_materials[INP].cb.alpha[1]);
        g_materials[GAN].cb.alpha[1] = pow(1. - g_materials[GAN].cb.mstar[1], 2.) / (g_materials[GAN].Eg + g_materials[GAN].cb.emin[1]); // expected value = 0.189
        printf("alphaK_gamma1[GAN] = %g\n", g_materials[GAN].cb.alpha[1]);


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
    g_materials[ALXINXSB].eps_static  = XVAL[ALXINXSB]  * g_materials[ALSB].eps_static + XVAL[ALXINXSB]  * g_materials[INSB].eps_static;
    g_materials[ALXIN1XSB].eps_static = XVAL[ALXIN1XSB] * g_materials[ALSB].eps_static + XVAL[ALXIN1XSB] * g_materials[INSB].eps_static;
    g_materials[INXGA1XAS].eps_static = XVAL[INXGA1XAS] * g_materials[INAS].eps_static + XVAL[INXGA1XAS] * g_materials[GAAS].eps_static;
    g_materials[INXAL1XAS].eps_static = XVAL[INXAL1XAS] * g_materials[INAS].eps_static + XVAL[INXAL1XAS] * g_materials[ALAS].eps_static;
    g_materials[INXGAXXAS].eps_static = XVAL[INXGAXXAS] * g_materials[INAS].eps_static + XVAL[INXGAXXAS] * g_materials[GAAS].eps_static;

    // semiconductor compounds high frequency dieletric constant
    g_materials[ALXINXSB].eps_hf  = XVAL[ALXINXSB]  * g_materials[ALSB].eps_hf + XVAL[ALXINXSB]  * g_materials[INSB].eps_hf;
    g_materials[ALXIN1XSB].eps_hf = XVAL[ALXIN1XSB] * g_materials[ALSB].eps_hf + XVAL[ALXIN1XSB] * g_materials[INSB].eps_hf;
    g_materials[INXGA1XAS].eps_hf = XVAL[INXGA1XAS] * g_materials[INAS].eps_hf + XVAL[INXGA1XAS] * g_materials[GAAS].eps_hf;
    g_materials[INXAL1XAS].eps_hf = XVAL[INXAL1XAS] * g_materials[INAS].eps_hf + XVAL[INXAL1XAS] * g_materials[ALAS].eps_hf;
    g_materials[INXGAXXAS].eps_hf = XVAL[INXGAXXAS] * g_materials[INAS].eps_hf + XVAL[INXGAXXAS] * g_materials[GAAS].eps_hf;

    // semiconductor compounds optical phonon scattering energy (eV)
    g_materials[ALXINXSB].hwo[0]=XVAL[ALXINXSB]*(g_materials[ALSB].hwo[0]+g_materials[INSB].hwo[0]);
    g_materials[ALXIN1XSB].hwo[0]=XVAL[ALXIN1XSB]*g_materials[ALSB].hwo[0]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].hwo[0];
    g_materials[INXGA1XAS].hwo[0]=XVAL[INXGA1XAS]*g_materials[INAS].hwo[0]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].hwo[0];
    g_materials[INXAL1XAS].hwo[0]=XVAL[INXAL1XAS]*g_materials[INAS].hwo[0]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].hwo[0];
    g_materials[INXGAXXAS].hwo[0]=XVAL[INXGAXXAS]*g_materials[INAS].hwo[0]+XVAL[INXGAXXAS]*g_materials[GAAS].hwo[0];
    // semiconductor compounds optical coupling constants (eV/m)
    g_materials[ALXINXSB].dtk[0]=XVAL[ALXINXSB]*(g_materials[ALSB].dtk[0]+g_materials[INSB].dtk[0]);
    g_materials[ALXIN1XSB].dtk[0]=XVAL[ALXIN1XSB]*g_materials[ALSB].dtk[0]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].dtk[0];
    g_materials[INXGA1XAS].dtk[0]=XVAL[INXGA1XAS]*g_materials[INAS].dtk[0]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].dtk[0];
    g_materials[INXAL1XAS].dtk[0]=XVAL[INXAL1XAS]*g_materials[INAS].dtk[0]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].dtk[0];
    g_materials[INXGAXXAS].dtk[0]=XVAL[INXGAXXAS]*g_materials[INAS].dtk[0]+XVAL[INXGAXXAS]*g_materials[GAAS].dtk[0];
    // semiconductor compounds optical phonon Z-factor
    g_materials[ALXINXSB].zf[0]=XVAL[ALXINXSB]*(g_materials[ALSB].zf[0]+g_materials[INSB].zf[0]);
    g_materials[ALXIN1XSB].zf[0]=XVAL[ALXIN1XSB]*g_materials[ALSB].zf[0]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].zf[0];
    g_materials[INXGA1XAS].zf[0]=XVAL[INXGA1XAS]*g_materials[INAS].zf[0]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].zf[0];
    g_materials[INXAL1XAS].zf[0]=XVAL[INXAL1XAS]*g_materials[INAS].zf[0]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].zf[0];
    g_materials[INXGAXXAS].zf[0]=XVAL[INXGAXXAS]*g_materials[INAS].zf[0]+XVAL[INXGAXXAS]*g_materials[GAAS].zf[0];
    // semiconductor compounds Crystal Density (Kg/m^3)
    g_materials[ALXINXSB].rho=XVAL[ALXINXSB]*(g_materials[ALSB].rho+g_materials[INSB].rho);
    g_materials[ALXIN1XSB].rho=XVAL[ALXIN1XSB]*g_materials[ALSB].rho+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].rho;
    g_materials[INXGA1XAS].rho=XVAL[INXGA1XAS]*g_materials[INAS].rho+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].rho;
    g_materials[INXAL1XAS].rho=XVAL[INXAL1XAS]*g_materials[INAS].rho+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].rho;
    g_materials[INXGAXXAS].rho=XVAL[INXGAXXAS]*g_materials[INAS].rho+XVAL[INXGAXXAS]*g_materials[GAAS].rho;
    // semiconductor compounds acoustic deformation potential (Joule)
    g_materials[ALXINXSB].da=XVAL[ALXINXSB]*(g_materials[ALSB].da+g_materials[INSB].da);
    g_materials[ALXIN1XSB].da=XVAL[ALXIN1XSB]*g_materials[ALSB].da+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].da;
    g_materials[INXGA1XAS].da=XVAL[INXGA1XAS]*g_materials[INAS].da+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].da;
    g_materials[INXAL1XAS].da=XVAL[INXAL1XAS]*g_materials[INAS].da+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].da;
    g_materials[INXGAXXAS].da=XVAL[INXGAXXAS]*g_materials[INAS].da+XVAL[INXGAXXAS]*g_materials[GAAS].da;
    // semiconductor compounds longitudinal sound velocity (m/sec)
    g_materials[ALXINXSB].ul=XVAL[ALXINXSB]*(g_materials[ALSB].ul+g_materials[INSB].ul);
    g_materials[ALXIN1XSB].ul=XVAL[ALXIN1XSB]*g_materials[ALSB].ul+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].ul;
    g_materials[INXGA1XAS].ul=XVAL[INXGA1XAS]*g_materials[INAS].ul+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].ul;
    g_materials[INXAL1XAS].ul=XVAL[INXAL1XAS]*g_materials[INAS].ul+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].ul;
    g_materials[INXGAXXAS].ul=XVAL[INXGAXXAS]*g_materials[INAS].ul+XVAL[INXGAXXAS]*g_materials[GAAS].ul;
    // semiconductor compounds energy gap
    g_materials[ALXINXSB].Eg=XVAL[ALXINXSB]*(g_materials[ALSB].Eg+g_materials[INSB].Eg);
    g_materials[ALXIN1XSB].Eg=XVAL[ALXIN1XSB]*g_materials[ALSB].Eg+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].Eg;
    g_materials[INXGA1XAS].Eg=XVAL[INXGA1XAS]*g_materials[INAS].Eg+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].Eg;
    g_materials[INXAL1XAS].Eg=XVAL[INXAL1XAS]*g_materials[INAS].Eg+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].Eg;
    g_materials[INXGAXXAS].Eg=XVAL[INXGAXXAS]*g_materials[INAS].Eg+XVAL[INXGAXXAS]*g_materials[GAAS].Eg;
    // semiconductor compounds energy minimum of GAMMA-valley
    g_materials[ALXINXSB].cb.emin[1]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.emin[1]+g_materials[INSB].cb.emin[1]);
    g_materials[ALXIN1XSB].cb.emin[1]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.emin[1]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.emin[1];
    g_materials[INXGA1XAS].cb.emin[1]=XVAL[INXGA1XAS]*g_materials[INAS].cb.emin[1]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.emin[1];
    g_materials[INXAL1XAS].cb.emin[1]=XVAL[INXAL1XAS]*g_materials[INAS].cb.emin[1]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.emin[1];
    g_materials[INXGAXXAS].cb.emin[1]=XVAL[INXGAXXAS]*g_materials[INAS].cb.emin[1]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.emin[1];
    // semiconductor compounds energy minimum 0f L-valley
    g_materials[ALXINXSB].cb.emin[2]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.emin[2]+g_materials[INSB].cb.emin[2]);
    g_materials[ALXIN1XSB].cb.emin[2]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.emin[2]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.emin[2];
    g_materials[INXGA1XAS].cb.emin[2]=XVAL[INXGA1XAS]*g_materials[INAS].cb.emin[2]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.emin[2];
    g_materials[INXAL1XAS].cb.emin[2]=XVAL[INXAL1XAS]*g_materials[INAS].cb.emin[2]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.emin[2];
    g_materials[INXGAXXAS].cb.emin[2]=XVAL[INXGAXXAS]*g_materials[INAS].cb.emin[2]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.emin[2];
    // GAMMA-valley effective mass
    g_materials[ALXINXSB].cb.mstar[1]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.mstar[1]+g_materials[INSB].cb.mstar[1]);
    g_materials[ALXIN1XSB].cb.mstar[1]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.mstar[1]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.mstar[1];
    g_materials[INXGA1XAS].cb.mstar[1]=XVAL[INXGA1XAS]*g_materials[INAS].cb.mstar[1]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.mstar[1];
    g_materials[INXAL1XAS].cb.mstar[1]=XVAL[INXAL1XAS]*g_materials[INAS].cb.mstar[1]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.mstar[1];
    g_materials[INXGAXXAS].cb.mstar[1]=XVAL[INXGAXXAS]*g_materials[INAS].cb.mstar[1]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.mstar[1];
    // L-valley effective mass
    g_materials[ALXINXSB].cb.mstar[2]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.mstar[2]+g_materials[INSB].cb.mstar[2]);
    g_materials[ALXIN1XSB].cb.mstar[2]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.mstar[2]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.mstar[2];
    g_materials[INXGA1XAS].cb.mstar[2]=XVAL[INXGA1XAS]*g_materials[INAS].cb.mstar[2]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.mstar[2];
    g_materials[INXAL1XAS].cb.mstar[2]=XVAL[INXAL1XAS]*g_materials[INAS].cb.mstar[2]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.mstar[2];
    g_materials[INXGAXXAS].cb.mstar[2]=XVAL[INXGAXXAS]*g_materials[INAS].cb.mstar[2]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.mstar[2];
    // non-parabolicity coefficient for semiconductor compounds in the GAMMA-valley
    g_materials[ALXINXSB].cb.alpha[1]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.alpha[1]+g_materials[INSB].cb.alpha[1]);
    g_materials[ALXIN1XSB].cb.alpha[1]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.alpha[1]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.alpha[1];
    g_materials[INXGA1XAS].cb.alpha[1]=XVAL[INXGA1XAS]*g_materials[INAS].cb.alpha[1]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.alpha[1];
    g_materials[INXAL1XAS].cb.alpha[1]=XVAL[INXAL1XAS]*g_materials[INAS].cb.alpha[1]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.alpha[1];
    g_materials[INXGAXXAS].cb.alpha[1]=XVAL[INXGAXXAS]*g_materials[INAS].cb.alpha[1]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.alpha[1];
    // non-parabolicity coefficient for Al_x In_(1-x) Sb in the L-valley
    g_materials[ALXINXSB].cb.alpha[2]=XVAL[ALXINXSB]*(g_materials[ALSB].cb.alpha[2]+g_materials[INSB].cb.alpha[2]);
    g_materials[ALXIN1XSB].cb.alpha[2]=XVAL[ALXIN1XSB]*g_materials[ALSB].cb.alpha[2]+(1.-XVAL[ALXIN1XSB])*g_materials[INSB].cb.alpha[2];
    g_materials[INXGA1XAS].cb.alpha[2]=XVAL[INXGA1XAS]*g_materials[INAS].cb.alpha[2]+(1.-XVAL[INXGA1XAS])*g_materials[GAAS].cb.alpha[2];
    g_materials[INXAL1XAS].cb.alpha[2]=XVAL[INXAL1XAS]*g_materials[INAS].cb.alpha[2]+(1.-XVAL[INXAL1XAS])*g_materials[ALAS].cb.alpha[2];
    g_materials[INXGAXXAS].cb.alpha[2]=XVAL[INXGAXXAS]*g_materials[INAS].cb.alpha[2]+XVAL[INXGAXXAS]*g_materials[GAAS].cb.alpha[2];
