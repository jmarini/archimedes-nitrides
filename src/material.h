/* particle.h -- This file is part of Archimedes release 1.2.0.
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


#ifndef ARCHIMEDES_MATERIAL_H
#define ARCHIMEDES_MATERIAL_H

#include "global_defines.h"


typedef struct {
    int num_valleys;

    double emin[MAX_VALLEYS];
    double mstar[MAX_VALLEYS];
    double alpha[MAX_VALLEYS];

    // precomputed constants
    double smh[MAX_VALLEYS];  // sqrt(2 * m* * m_e * q) / hbar
    double hhm[MAX_VALLEYS];  // hbar^2 / (2 * m* * m_e * q)
    double hm[MAX_VALLEYS];   // hbar / (m* * m_e)
} Band_Info;


typedef struct {
    int id;

    double Eg;             // band gap [eV]
    double affinity;       // electron affinity [eV]
    Band_Info cb;
    Band_Info vb;

    int  zscatter[MAX_VALLEYS][MAX_VALLEYS];

    double abs_correction; // absorption coefficient correction
    double eps_static;     // static dielectric constant (relative)
    double eps_hf;         // high frequency dielectric constant (relative)
    double hwo[6];         // optical phonon energy [eV]
    double dtk[6];         // optical coupling constant [eV/m]
    int    zf[6];          // optical phonon z-factor
    double da;             // acoustic deformation potential [J]
    double ul;             // longitudinal sound velocity [m/s]
    double rho;            // crystal density [kg/m^3]
    double lattice_const;  // crystal lattice constant [m]
    double kav;            // electro-mechanical coupling constant

} Material;


extern Material g_materials[NOAMTIA];


Material material_node(int i, int j);

char* mc_material_name(Material *material);


#endif
